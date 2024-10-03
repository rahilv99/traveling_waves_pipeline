%% 1. WAVE EXTRACTION
% inputs: raw data ouput from attn processing pipeline (wave_info struct)
% ouputs: all_pres and all_recog structs
% fields of output include processed data for plotting (refer to end of code)
clear; close all;

addpath(genpath('/Volumes/Rahil_FRNU/ZaghloulCodebase'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/macro_tws_pal'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/circstat-matlab-master'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/(7) Analysis/functions'))

eeg_path = fullfile('/Volumes','FRNU/','data/', 'eeg/');

load("/Volumes/Rahil_FRNU/Scripts/soc.mat")
load("/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct_full.mat")

cycle_min = 1; % this seems like  reasonable minimum, but could make it longer if we're interested in more sustained waves. perhaps more sustained waves serve a different purpose in attention?
nperm = 100; % change if necessary
z_thresh = 2; % change if necessary
el_dist= 25; 
delete(gcp('nocreate'));
parpool('local',10)

p=1; r=1; f=1;

for i = 1:length(soc)
    sess_i = soc(i);

    curr = s(sess_i); %sess_i = 311 

    subj = curr.pat_ID;
    sess = curr.session;

    try
        mono_elec = readtable(fullfile(eeg_path, subj,'/tal/atlas/atlas_monopolar_simple.csv'));

        folder_name = fullfile('/Volumes/Rahil_FRNU/Data', subj,sess);
        %%
        load(fullfile(folder_name, 'allwave.mat')) % wave_info
        load(fullfile(folder_name, 'clusters.mat'))  % sess_clust
        load(fullfile(folder_name, 'pow.mat')) % sess_power_struct
        load(fullfile(folder_name,'sess_struct_pc.mat')) % sess_struct

        out = sprintf("Loaded files for patient %s, %s", subj, sess);
        disp(out)
        %%
        % load + clean events
        load(fullfile(eeg_path, subj,'/behavioral/attentionTask/', sess,'/events.mat'));

        default_root_EEG_dir = '/Volumes/Shares/FRNU/data/eeg';
        eeg_file = regexprep(unique({events.eegfile}), default_root_EEG_dir,eeg_path);
        eeg_file = eeg_file( ~cellfun(@isempty,eeg_file)  );%remove empty eegfiles references
        eeg_file = regexprep(eeg_file,'noreref','processedBP');
        events = setField(events,'eegfile',eeg_file{1});

        waveletFreqs = logspace(log10(2),log10(32),200); % 200 wavelets from 2-32 Hz, log space
        resamp=1000;
        w  = 6;
        
        % step to match clusters between encoding, cue, and retrieval

        clust_sets = match_clusters_across_task(sess_clust);
        for cl = 1:size(clust_sets,1)

            try
            for n = 1:3
                % check if there is a cluster in the encoding, cue or retrieval
                % columns
                if clust_sets(cl,n)>0

                    cln =  clust_sets(cl,n);

                    if strcmp(wave_info(cln).task_phase, 'presentation')
                        clust_el = sess_clust(cln).pres_channame;
                        correct = [sess_struct.evs_pres.responseCorrect];
                        cued = [strcmp({sess_struct.evs_pres.asteriskStr},'*before')];
                        n_cued = [strcmp({sess_struct.evs_pres.asteriskStr},'*none')];
                        p_ep_clust = zeros(length([sess_struct.evs_pres.RT]), size(sess_struct.raw_pres,3)); % events x [time points + 1 s buffer on either side]
                        raw_eeg = sess_struct.raw_pres;
                    end
                    if strcmp(wave_info(cln).task_phase, 'recog')
                        clust_el = sess_clust(cln).recog_channame;
                        correct = [sess_struct.evs_pres.responseCorrect];
                        cued = [strcmp({sess_struct.evs_pres.asteriskStr},'*before')];
                        n_cued = [strcmp({sess_struct.evs_pres.asteriskStr},'*none')];
                        p_ep_clust = zeros(length([sess_struct.evs_recog.RT]), size(sess_struct.raw_recog,3));
                        raw_eeg = sess_struct.raw_recog;
                    end
                    %% USE EEGOFFSETS TO OBTAIN FREE RECALL WORDS -- NEED TO IMPLEMENT TO GET FREE RECALL DATA. for now, continue statement passes to next cluster
                    if strcmp(wave_info(cln).task_phase, 'free recall')
                        continue;
                        clust_el = sess_clust(cln).fr_channame;
                        correct = [sess_struct.evs_fr.FRannotateCounts(2)];
                        p_ep_clust = zeros(length([sess_struct.evs_fr.RT]), size(sess_struct.raw_fr,3));
                        raw_eeg = sess_struct.raw_fr;
                    end
%% Section 1: Prevalence Script
                    % Desikan mapping (group by lobes)
                    tal = wave_info(cln).cltal;
    
                    tw_reg = mono_elec(ismember(mono_elec.chanName, clust_el),{'chanName','label_desikan'});
                    tw_reg.y = tal(ismember(tw_reg.chanName, clust_el),2);

                    lobes = cellfun(@(x,y) aparc2lobe_atlptlmtl(x,y),tw_reg{:,2}, num2cell(tw_reg{:,3}),'UniformOutput',false); 
                    lobes = cell2table(lobes,"VariableNames","lobes");
                    lobes = groupcounts(lobes,'lobes');
                    
                    clust_el_nums = ismember(sess_struct.chan_names, clust_el);
                    num_elec = .5*sum(clust_el_nums);
                    freq_ind = [waveletFreqs > wave_info(cln).mfreq - 2 & waveletFreqs < wave_info(cln).mfreq + 2 ]; % frequency range around peak frequency, index of wavelets
                    dur_th = cycle_min*(1000/wave_info(cln).mfreq); % duration threshold using min num of cycles of oscillation - used in p_epi_withminoscduration funciton below
                    %% refined wave_logic 

                    %pow_pgd_logic = zeros(length(wave_info(cl).rt), size(wave_info(cl).CLcorrelation,2));

                    parfor ev = 1:size(p_ep_clust,1)
                        unionvec=calcPepisode_palwaves(squeeze(raw_eeg(clust_el_nums,ev,:)),events(1), waveletFreqs, w); % cluster elecs x wavelets x time
                        unionvec_freqcollapse = squeeze(any(unionvec(:,freq_ind,:),2)); %
                        p_ep_clust(ev,:) = squeeze(sum(unionvec_freqcollapse,1)>num_elec);
                    end
                    pgd_logic_ev = squeeze(nanmean(wave_info(cln).CLcorrelation,1));
                    if strcmp(wave_info(cln).task_phase, 'presentation')
                        p_ep_clust = p_ep_clust(:,1001:end-2000);
                        task_phase_str = 'Presentation';

                    end
                    if strcmp(wave_info(cln).task_phase, 'recog')
                        p_ep_clust = p_ep_clust(:,1001:end-2000);
                        task_phase_str = 'Recognition';
                    end
                    if strcmp(wave_info(cln).task_phase, 'fr')
                        p_ep_clust = p_ep_clust(:,1001:end-1500); 
                        task_phase_str = 'Free Recall';
                    end

                    [p_ep_withmindur] = p_epi_withminoscduration(p_ep_clust, dur_th); %this function zeros out any parts of p episode matrix that are less than minumum duration set
                    %shuffling pgd function
                    %[new_pgd_logic] = pgd_shuffled(wave_info(cln), nperm, z_thresh, el_dist);

                    load(fullfile(folder_name,'shuffled',cln)); % contains variable 'new_pgd_logic' calculated using biowulf
                    pow_pgd_logic = new_pgd_logic  & p_ep_clust;
                    dip_th= 5; 
                    pow_pgd_logic_nodip = wave_logic_allow_dips(pow_pgd_logic, dip_th); % this is just like the function above except it fills little dips (0s) to 1s

                    cued = (cued == 1 & correct == 1);
                    n_cued = (n_cued == 1 & correct == 1);

% figures originally in script (IGNORE)
                    % for each cluster, this plots for encoding, cue, and retrieval
                    % % of trials with oscillations
                    % avg PGD
                    % % of trials with waves (osc + pgd)                    
%                     figure(cl)
%                     set(gcf,'color','w');
% 
%                     % 1. % of trials with an oscillation
%                     subplot(3,3,n)
%                     plot(100*(sum(p_ep_clust(cued,:),1)/size(p_ep_clust(cued,:),1)), 'g'); hold on;
%                     plot(100*(sum(p_ep_clust(n_cued,:),1)/size(p_ep_clust(n_cued,:),1)), 'b');
% 
%                     ylabel('% trials with oscillation')
%                     legend(['cued, ', num2str(length(cued))], ['no cue, ', num2str(length(n_cued))])
%                     legend boxoff
%                     % Fix this to reflect pres, recog, and fr related events
%                     if n == 1
%                         rectangle('Position', [2000 0 1000 100], 'FaceColor', [0 0 0.0 0.05])
%                     end
%                     if n == 2
%                         rectangle('Position', [1000 0 rt 100], 'FaceColor', [0 0 0.0 0.05])
%                     end
%                     if n == 3
%                         plot([5000 5000], [0 100], 'k', 'LineWidth', 2)
%                     end
%                     ylim([100*min(nanmean(p_ep_clust,2)) 100*max(nanmean(p_ep_clust,2))])
% 
%                     title(['OSCILLATIONS during ', task_phase_str])
                    % 2. avg wave "strength" (CLcorrelation)
%                     subplot(3,3,n+3)
%                     shadedErrorBar([],squeeze(nanmean(nanmean(wave_info(cln).CLcorrelation(:,:,cued),1),3)),std(squeeze(nanmean(wave_info(cln).CLcorrelation(:,:,cued),1))')/sqrt(sum(cued)),{'g','markerfacecolor',[0.4660 0.6740 0.1880]},1);
%                     hold on;
%                     shadedErrorBar([],squeeze(nanmean(nanmean(wave_info(cln).CLcorrelation(:,:,n_cued),1),3)),std(squeeze(nanmean(wave_info(cln).CLcorrelation(:,:,n_cued),1))')/sqrt(sum(n_cued)),{'b','markerfacecolor',[0 0.4470 0.7410]},1);
%                     ylabel('average PGD across channels')
%                     if n == 1
%                         rectangle('Position', [2000 0 1000 1], 'FaceColor', [0 0 0.0 0.05])
%                     end
%                     if n == 2
%                         rectangle('Position', [1000 0 rt 1], 'FaceColor', [0 0 0.0 0.05])
%                     end
%                     if n == 3
%                         plot([5000 5000], [0 1], 'k', 'LineWidth', 2)
%                     end
%                     ylim([min(nanmean(nanmean(wave_info(cln).CLcorrelation,1),3))-.1 max(nanmean(nanmean(wave_info(cln).CLcorrelation,1),3))+.1])
%                     title(['PGD during ', task_phase_str])
                    % 3. % of trials with wave
%                     subplot(3,3,n+6)
%                     plot(100*(sum(pow_pgd_logic(cued,:),1)/size(pow_pgd_logic(cued,:),1)), 'g'); hold on;
%                     plot(100*(sum(pow_pgd_logic(n_cued,:),1)/size(pow_pgd_logic(n_cued,:),1)), 'b');
%                     ylabel('% trials PGD and oscillations')
%                     if n == 1
%                         rectangle('Position', [2000 0 1000 100], 'FaceColor', [0 0 0.0 0.05])
%                     end
%                     if n == 2
%                         rectangle('Position', [1000 0 rt 100], 'FaceColor', [0 0 0.0 0.05])
%                     end
%                     if n == 3
%                         plot([5000 5000], [0 100], 'k', 'LineWidth', 2)
%                     end
%                     ylim([100*min(nanmean(pow_pgd_logic,2)) 100*max(nanmean(pow_pgd_logic,2))])
% 
%                     title(['WAVES during ', task_phase_str])
% 
%                     sgtitle([subj, ': ', num2str(wave_info(cln).mfreq), ' Hz cluster in ', strjoin(unique(tw_reg.lobe_desikan))])
%                     set(gcf,'position',[0,0,1900,1500]);
% 
%                     out_dir = fullfile('/Volumes/frnu/Rahil/Data/', subj, sess, 'vis');
%                     if ~exist(out_dir, 'dir') % make folder if doesn't exist already
%                         mkdir(out_dir);
%                     end
% 
%                     save(fullfile(folder_name, "allwave.mat"), 'wave_info', '-v7.3');
                
%% Section 2: Direction Script
                    wave_dir = squeeze(circ_mean(wave_info(cln).ang,[],1))';
                    wave_logic = pow_pgd_logic_nodip;

                    trial_pref_dir = nan(1,size(wave_dir,1));
                    for t = 1:size(wave_dir,1)
                        if sum(wave_logic(t,:)) > 0
                            trial_dirs = wave_dir(t,wave_logic(t,:));
                            converted_trial_dirs = convert_angle_to_rad_in_yzplane(wave_info(cln).cltal, trial_dirs); 
                            [thetahat kappa] = circ_vmpar(converted_trial_dirs);
                            trial_pref_dir(t) = thetahat;
                        end
                    end
                    
                   % Percent preferred direction
                   yz_dirs = convert_angle_to_rad_in_yzplane(wave_info(cln).cltal, wave_dir(wave_logic)); 
                    [pref_dir kappa] = circ_vmpar(yz_dirs);
                    if pref_dir < 0
                        pref_dir = pref_dir + 2*pi;
                    end
                    a1 = pref_dir - pi/2;
                    a2 = pref_dir + pi/2;
                    if a1 < 0
                        a1 = a1 + 2*pi;
                    end
                    if a2 < 0
                        a2 = a2 + 2*pi;
                    end
                    if a2 > 2*pi
                        a2 = a2 - 2*pi;
                    end
                    percent_pref_dir_corr = [];
                    percent_pref_dir_incorr = [];

                    for tt = 1:size(wave_dir,2)
                        dirs_t = wave_dir(wave_logic(:,tt)' & cued,tt);
                        dirs_t_uncued = wave_dir(wave_logic(:,tt)' & n_cued,tt);
                        if ~isempty(dirs_t)
                            dirs_t = convert_angle_to_rad_in_yzplane(wave_info(cln).cltal, dirs_t); 
                        end
                        if ~isempty(dirs_t_uncued)
                            dirs_t_uncued = convert_angle_to_rad_in_yzplane(wave_info(cln).cltal, dirs_t_uncued); 
                        end
                        if pref_dir < pi/2
                            percent_pref_dir_corr(tt) = sum(dirs_t>a1 | dirs_t < a2)/length(dirs_t);
                            percent_pref_dir_incorr(tt) = sum(dirs_t_uncued>a1 | dirs_t_uncued < a2)/length(dirs_t_uncued);
                        elseif pref_dir > 3*pi/2
                            percent_pref_dir_corr(tt) = sum(dirs_t>a1 | dirs_t < a2)/length(dirs_t);
                            percent_pref_dir_incorr(tt) = sum(dirs_t_uncued>a1 | dirs_t_uncued < a2)/length(dirs_t_uncued);
                        else
                            percent_pref_dir_corr(tt) = sum(dirs_t>a1 & dirs_t < a2)/length(dirs_t);
                            percent_pref_dir_incorr(tt) = sum(dirs_t_uncued>a1 & dirs_t_uncued < a2)/length(dirs_t_uncued);
                        end
                    end
                
                    dc_corr = [];
                    dc_incorr = [];
                    for tt = 1:size(wave_dir,2)
                        dirs_t = wave_dir(wave_logic(:,tt)' & cued,tt);
                        if ~isempty(dirs_t)
                            dirs_t = convert_angle_to_rad_in_yzplane(wave_info(cln).cltal,dirs_t);
                        end
                        dc_corr(tt) = rbar(dirs_t);

                        dirs_t = wave_dir(wave_logic(:,tt)' & n_cued,tt);
                        if ~isempty(dirs_t)
                            dirs_t = convert_angle_to_rad_in_yzplane(wave_info(cln).cltal,dirs_t);
                        end
                        dc_incorr(tt) = rbar(dirs_t);
                    end

                    % DC across cluster, averaged across trials over time
                    dcpl_corr = [];
                    dcpl_incorr = [];
                    wave_dir_el = wave_info(cln).ang;

                    for tt = 1:size(wave_dir_el,2)
                        dirs_t = wave_dir_el(:,tt,wave_logic(:,tt)' & cued);

                        for ii = 1:size(dirs_t,1)
                            if ~isempty(dirs_t(ii,:,:))
                                dirs_t(ii,:) = convert_angle_to_rad_in_yzplane(wave_info(cln).cltal,dirs_t(ii,:,:));
                            end
                        end
                        dcpl_corr(tt) = squeeze(nanmean(rbar(dirs_t))); 

                        dirs_t = wave_dir_el(:,tt,wave_logic(:,tt)' & n_cued);
                        dcpl_incorr(tt) = squeeze(nanmean(rbar(dirs_t)));
                    end

                    phv_corr = [];
                    phv_incorr = [];
                    for tt = 1:size(wave_dir,2)
                        phv_corr(tt) = wave_info(cln).mfreq/(10*squeeze(nanmean(nanmean(wave_info(cln).spatial_frequency(:,tt,wave_logic(:,tt)' & cued)))));
                        phv_incorr(tt) = wave_info(cln).mfreq/(10*squeeze(nanmean(nanmean(wave_info(cln).spatial_frequency(:,tt,wave_logic(:,tt)' & n_cued)))));
                    end


 %% Section 3: True Angle
                    pres_cued_words = find(cued);
                    pres_n_cued_words = find(n_cued);
                    
                    recog_RT = {sess_struct.evs_recog.RT};
                    
                    cued_words = {sess_struct.evs_pres(cued).encodedWord};
                    n_cued_words = {sess_struct.evs_pres(n_cued).encodedWord};

                    cued_cell = cell([length(cued_words),2]);
                    n_cued_cell = cell([length(n_cued_words),2]);



                    for word = 1:length(cued_words)
                        % index of correctly recalled cued word
                        pres_word = pres_cued_words(word);
                        recog_word = find([strcmp({sess_struct.evs_recog.targetWord},cued_words(word))]);
                        RT = cell2mat(recog_RT(recog_word));

                        % Check if wave occurs in pres and recog
                        if  RT < 4000
                            lb = max(RT,1000)+1;

                            % get angles during cue, pres and recog for current word
                            if strcmp(task_phase_str, 'Presentation')
                                dir = squeeze(circ_mean(wave_info(cln).ang(:, 1:3000, pres_word),[],1))';
                                cued_cell(word,2) =  {wave_logic(pres_word,:)};
                            end
                            if strcmp(task_phase_str, 'Recognition')
                                dir = squeeze(circ_mean(wave_info(cln).ang(:, lb:(1000+RT), recog_word),[],1))';
                                cued_cell(word,2) =  {wave_logic(recog_word,:)};
                            end
                            cued_ang = convert_angle_to_rad_in_yzplane(wave_info(clust_sets(cl,1)).cltal, dir)';
                            cued_cell(word,1) = {cued_ang};
                        end
                    end

                    for word = 1:length(n_cued_words)
                        % index of correctly recalled cued word
                        pres_word = pres_cued_words(word);
                        recog_word = find([strcmp({sess_struct.evs_recog.targetWord},cued_words(word))]);
                        RT = cell2mat(recog_RT(recog_word));

                        % Check if wave occurs in pres and recog
                        if  RT < 4000
                            lb = max(RT,1000)+1;

                            % get angles during cue, pres and recog for current word
                            if strcmp(task_phase_str, 'Presentation')
                                dir = squeeze(circ_mean(wave_info(cln).ang(:, 1:3000, pres_word),[],1))';
                                n_cued_cell(word,2) =  {wave_logic(pres_word,:)};
                            end
                            if strcmp(task_phase_str, 'Recognition')
                                dir = squeeze(circ_mean(wave_info(cln).ang(:, lb:(1000+RT), recog_word),[],1))';
                                n_cued_cell(word,2) =  {wave_logic(recog_word,:)};
                            end
                            cued_ang = convert_angle_to_rad_in_yzplane(wave_info(clust_sets(cl,1)).cltal, dir)';
                            n_cued_cell(word,1) = {cued_ang};
                        end
                    end

%% Section 4: Wave Prevalence Type
                    cp = cued_cell(:,2);
                    nonEmptyIndex = ~cellfun(@isempty, cp);
                    waves = cp(nonEmptyIndex);
                
                
                    % Sliding window
                    window_size = 500;
                    step_size = 20;
                    prevalence_len = nan([((3000-window_size)/step_size +1),1]);
                    prevalence_num = nan([((3000-window_size)/step_size +1),1]);
                    prevalence_t = 0:step_size:3000-window_size;
                
                    for w = 1:length(prevalence_t)
                        t_i = prevalence_t(w);
                        t_f = t_i + window_size;
                
                        evs = cellfun(@(x) (x(1,t_i+1:t_f)), waves, 'UniformOutput', false);
                        evs = cellfun(@bwconncomp, evs, 'UniformOutput', false);
                        % CC is a struct containing NumObjects (number of waves) and
                        % PixelIdxList (cell array with linear indicies of each wave)
                    
                        epoch_info = struct();
                    
                        for tr = 1:length(evs)
                    
                            CC = evs{tr,1};
                            for wav = 1:length(CC)
                                len = cellfun(@length,CC.PixelIdxList,'UniformOutput', false);
                                wave_idx = cell2mat(len(cell2mat(len) > 5));
                    
                                % populate struct for current epoch
                                epoch_info(tr).len = mean(wave_idx);
                                epoch_info(tr).num = length(wave_idx);
                            end
                        end
                    
                        avg_len = nanmean([epoch_info.len]);
                        avg_num = mean([epoch_info.num]);
                    
                        prevalence_len(w) = avg_len;
                        prevalence_num(w) = avg_num;
                    end
                    
                    temp_struct.cued_prev_len = prevalence_len;
                    temp_struct.cued_prev_num = prevalence_num;
                    temp_struct.prev_t = prevalence_t;
                
                    np = n_cued_cell(:,2);
                    nonEmptyIndex = ~cellfun(@isempty, np);
                    waves = np(nonEmptyIndex);
                
                    prevalence_len = nan([((3000-window_size)/step_size +1),1]);
                    prevalence_num = nan([((3000-window_size)/step_size +1),1]);
                
                    for w = 1:length(prevalence_t)
                        t_i = prevalence_t(w);
                        t_f = t_i + window_size;
                
                        evs = cellfun(@(x) (x(1,t_i+1:t_f)), waves, 'UniformOutput', false);
                        evs = cellfun(@bwconncomp, evs, 'UniformOutput', false);
                        % CC is a struct containing NumObjects (number of waves) and
                        % PixelIdxList (cell array with linear indicies of each wave)
                    
                        epoch_info = struct();
                    
                        for tr = 1:length(evs)
                    
                            CC = evs{tr,1};
                            for wav = 1:length(CC)
                                len = cellfun(@length,CC.PixelIdxList,'UniformOutput', false);
                                wave_idx = cell2mat(len(cell2mat(len) > 5));
                    
                                % populate struct for current epoch
                                epoch_info(tr).len = mean(wave_idx);
                                epoch_info(tr).num = length(wave_idx);
                            end
                        end
                    
                        avg_len = nanmean([epoch_info.len]);
                        avg_num = mean([epoch_info.num]);
                    
                        prevalence_len(w) = avg_len;
                        prevalence_num(w) = avg_num;
                    end
                    
                    temp_struct.n_cued_prev_len = prevalence_len;
                    temp_struct.n_cued_prev_num = prevalence_num;


%% Results of processing raw data from wave_info
                    % fill all_pres, all_recog, or all_fr
                    temp_struct.sub = wave_info(cln).sub;
                    temp_struct.sess = wave_info(cln).sess;
                    temp_struct.cltal = wave_info(cln).cltal;
                    temp_struct.mfreq = wave_info(cln).mfreq;
                    temp_struct.CLcorrelation = wave_info(cln).CLcorrelation;
                    temp_struct.ang_raw = wave_info(cln).ang;

                    % 1st section
                    temp_struct.wave_logic = pow_pgd_logic;
                    temp_struct.p_episode = p_ep_clust;
                    temp_struct.rt = mean([sess_struct.evs_pres.RT]);
                    temp_struct.cued = cued;
                    temp_struct.n_cued = n_cued;
                    temp_struct.loc = lobes.lobes;
                    temp_struct.loc_percent = lobes.Percent;
                    temp_struct.rt = {-999};
                    temp_struct.elec_names = string(clust_el(ismember(clust_el,mono_elec.chanName)));
                    

                    % 2nd section
                    temp_struct.trial_dir_cued = trial_pref_dir(~isnan(trial_pref_dir) & cued);
                    temp_struct.trial_dir_n_cued = trial_pref_dir(~isnan(trial_pref_dir)  & n_cued);
                    temp_struct.percent_pref_dir_cued = 100*percent_pref_dir_corr;
                    temp_struct.percent_pref_dir_n_cued = 100*percent_pref_dir_incorr;
                    temp_struct.dc_cued = dc_corr;
                    temp_struct.dc_n_cued = dc_incorr;
                    temp_struct.dc_inclust_cued = dcpl_corr;
                    temp_struct.dc_inclust_n_cued = dcpl_incorr;
                    temp_struct.phv_corr = phv_corr;
                    temp_struct.phv_incorr = phv_incorr;

                    %3rd section (corrected angle)
                    % This gives us access to the full epoch if we ever want to do something other than average (below)
                    temp_struct.cued_ang = cued_cell;
                    temp_struct.n_cued_ang = n_cued_cell;
                    
                    if strcmp(task_phase_str, 'Presentation')
                        % If presentation epoch, get only angle during encoding.
                        waves = cellfun(@(x, idx) (x(2001:3000)), cued_cell(:,1), cued_cell(:,2), 'UniformOutput', false);
                        waves = cellfun(@(x, idx) (x(idx(2001:3000))), cued_cell(:,1), cued_cell(:,2), 'UniformOutput', false);
                        avg = cellfun(@(x) (circ_mean(x)),waves,'UniformOutput', false);
                        nonEmptyIndex = ~cellfun(@isempty, avg);
                        temp_struct.cued_avg = cell2mat(avg(nonEmptyIndex))';

                        waves = cellfun(@(x, idx) (x(2001:3000)), cued_cell(:,1), cued_cell(:,2), 'UniformOutput', false);
                        waves = cellfun(@(x, idx) (x(idx(2001:3000))), n_cued_cell(:,1), n_cued_cell(:,2), 'UniformOutput', false);
                        avg = cellfun(@(x) (circ_mean(x)),waves,'UniformOutput', false);
                        nonEmptyIndex = ~cellfun(@isempty, avg);
                        temp_struct.n_cued_avg = cell2mat(avg(nonEmptyIndex))';
                    else
                        waves = cellfun(@(x, idx) (x(idx)), cued_cell(:,1), cued_cell(:,2), 'UniformOutput', false);
                        avg = cellfun(@(x) (circ_mean(x)),waves,'UniformOutput', false);
                        nonEmptyIndex = ~cellfun(@isempty, avg);
                        temp_struct.cued_avg = cell2mat(avg(nonEmptyIndex))';

                        waves = cellfun(@(x, idx) (x(idx)), n_cued_cell(:,1), n_cued_cell(:,2), 'UniformOutput', false);
                        avg = cellfun(@(x) (circ_mean(x)),waves,'UniformOutput', false);
                        nonEmptyIndex = ~cellfun(@isempty, avg);
                        temp_struct.n_cued_avg = cell2mat(avg(nonEmptyIndex))';
                    end
                    temp_struct.dif = circ_mean(temp_struct.cued_avg) - circ_mean(temp_struct.n_cued_avg);

                    % Angle Plot
                    % Obtain wave angles by ang*wave_info
                    C = cued_cell;
                    waves = cellfun(@(x, idx) (x.*idx), C(:,1), C(:,2), 'UniformOutput', false);
                    nonEmptyIndex = ~cellfun(@isempty, waves);
                    waves = cell2mat(waves(nonEmptyIndex));
                    waves(waves==0)= NaN;
                    temp_struct.cued_signal = circ_mean_nan(waves',2);
                    % Obtain wave angles by ang*wave_info
                    C = n_cued_cell;
                    waves = cellfun(@(x, idx) (x.*idx), C(:,1), C(:,2), 'UniformOutput', false);
                    nonEmptyIndex = ~cellfun(@isempty, waves);
                    waves = cell2mat(waves(nonEmptyIndex));
                    waves(waves==0)= NaN;
                    temp_struct.n_cued_signal = circ_mean_nan(waves',2);

                    
                    if strcmp(task_phase_str, 'Presentation')
                        all_pres(p) = temp_struct;
                        p = p+1;
                    end
                    if strcmp(task_phase_str, 'Recognition')
                        all_recog(r) = temp_struct;
                        r = r+1;
                    end

                    if strcmp(task_phase_str, 'Free Recall')
                        all_fr(f) = temp_struct;
                        f = f+1;
                    end

                    clear cued; clear n_cued; clear pow_pgd_logic; clear p_ep_clust;
                end
            end
            catch
                fprintf("Error on n %d in cl %d\n", n, cl);
            end
        end
    catch e
        fprintf("Error on patient %s, %s\n", subj, sess);
        fprintf(1,'The identifier was: %s\n',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
if exist(all_pres,'var')
    save('/Volumes/Rahil_FRNU/Scripts/all_pres.mat','all_pres', '-v7.3');
end
if exist(all_recog,'var')
    save('/Volumes/Rahil_FRNU/Scripts/all_recog.mat','all_recog', '-v7.3');
end
if exist(all_fr,'var')
    save('/Volumes/Rahil_FRNU/Scripts/all_fr.mat','all_fr', '-v7.3');
end

delete(gcp('nocreate'));