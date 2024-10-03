%% 5a. Preprocessing - Get data to feed into circular linear regression model (local)
% Adapated from pal_phase_byclust from honghui's code to calculate phase waves in sternberg wm task and Uma's previous code to calculate phase waves in TH
% Rahil Verma 6/18/2024 rahil.verma@duke.edu

%% setting up + initializing

% set path
clear; close all;

%addpath(genpath('/Users/mohanur/Documents/MATLAB/eeg_toolbox_rhino'));
%addpath(genpath('/Users/mohanur/Documents/MATLAB/iEEGDataBase'));
%addpath(genpath('/Users/mohanur/Documents/MATLAB/circstat-matlab-master'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/ZaghloulCodebase'));
%addpath(genpath('Z:\people\Rahil\ZaghloulCodebase'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/circstat-matlab-master'));
%addpath(genpath('Z:\people\Rahil\Scripts\circstat-matlab-master'))

projfolder='/Volumes/Rahil_FRNU/Scripts';
%projfolder = 'Z:\people\Rahil\Scripts';

% loading sessions structures (created separately)
load("/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct2.mat")
s=s2;
%load('Z:\people\Rahil\Scripts\attn_sessions_struct.mat');

%addAttachedFiles(gcp,['circ_lin_regress_2D.m'])

%% Biowolf text instructions
% set biowulf directory (make sure you have CMND + K’ed already) to:
% smb://hpcdrive.nih.gov/data
%biowulfSaveDir = 'Y:\Rahil\Scripts\Code\processing';
biowulfSaveDir = '/Volumes/frnu-1/Rahil/Processing';
if ~exist(biowulfSaveDir, 'dir') % make folder if doesn't exist already
    mkdir(biowulfSaveDir);
end
cd(biowulfSaveDir)

%% Fill out important information here:
txtFile.name = 'vermar5_attn_pipeline_processing.swarm18-43_p2.txt';
txtFile.jobName = 'vermar5_attn_pipeline_processing';
txtFile.numTasks = 2; % should always be 2 * number of parallel workers (plus a few more to not cause delays or crashes)
txtFile.GB = 6; % some of the 'norm' partitions have a maximum memory of 247 GB. Keep this number under 247 for smaller wait times.
% if you need to go above 350 GB, use the 'largemem' partition
txtFile.scratch = '--gres=lscratch:4'; % it's a good idea to ask for the same amount of /lscratch/ space as you need for the function itself
txtFile.function = 'attn_fitting_biowolf_v1'; % the name of the function
txtFile.time = '01:45:00'; % anything that'll take more than 10 days should use the 'unlimited' partition. I'm pretty sure the format is: DAYS:HOURS:MINUTES:SECONDS

txtFile.introText = 'matlab -nodisplay -r';
txtFile.outroText = 'exit;" ';
txtFile.dir = '/data/FRNU/Rahil/Processing'; % Biowulf location where you store your script

fileID = fopen(txtFile.name, 'w');

% write first few required lines
fprintf(fileID, '#SWARM --job-name=%s\n', txtFile.jobName);
fprintf(fileID, '#SWARM --module=matlab\n');
fprintf(fileID, '#SWARM -t %d -g %d --time %s %s\n', ...
    txtFile.numTasks, txtFile.GB, txtFile.time, txtFile.scratch);

fprintf(fileID, '\n \n'); % line break

%% loop through sessions + load  power structs

waveletFreqs = logspace(log10(2),log10(32),200); % 200 wavelets from 2-32 Hz, log space
resampleFreq= 1000;
el_dist = 25;

for sess_i = 169:length(s)
    patID = s(sess_i).pat_ID;
    sessID = s(sess_i).session;
    % Mac path to data
    folder_name = fullfile('/Volumes/frnu-1/Rahil/Data', patID,sessID);
    %folder_name = fullfile('Y:\Rahil\Data', allFileNames{sess_i});
    % Biowolf path to data
    data_path = fullfile('/data/FRNU/Rahil/Data/', patID, sessID);


    fprintf("Processing session %d \n", sess_i);

    try
        % load session cluster data
        load(fullfile(folder_name, "clusters.mat"));
        load(fullfile(folder_name, "sess_struct_pc.mat"));
    catch
        fprintf('Files not found for %s, %s \n',patID,sessID);
        continue;
    end

    ct = 1;

    %% loop through pres clusters
    try
        pres_cl = sess_clust(find(~cellfun(@isempty,{sess_clust.pres_freq})));
        pres_cl = pres_cl(~cellfun(@(C) isnan(C), {pres_cl.pres_freq}));
        for cl = 1:length(pres_cl)
            mf = pres_cl(cl).pres_freq;
            tal = pres_cl(cl).pres_tal;
            elecs = pres_cl(cl).pres_elecnum;
            channames = pres_cl(cl).pres_channame;
            pres_clust_eeg = sess_struct.raw_pres(ismember(sess_struct.chan_names, pres_cl(cl).pres_channame), :,:);
            pres_phases =  [];
            pres_amp = [];
            for el = 1:length(elecs)
                eegFilt = buttfilt(squeeze(pres_clust_eeg(el,:,:)),[mf*.85 mf/.85],resampleFreq,'bandpass',2); %             narrowband filtering for phase
                %             eegFilt = buttfilt(squeeze(enc_clust_eeg(el,:,:)),[2 32],resampleFreq,'bandpass',2); %generalized phase
                h = hilbert(eegFilt')'; %hilbert only goes across columns, unfortunately...
                h = h(:,1001:end-1000); %remove the buffering
                pres_phases(el,:,:) = angle(h); % chan x event x time
                pres_amp(el,:,:) = abs(h);
            end
    
            [ch,eb,tm]=size(pres_phases);
    
            ref_phase = NaN(eb,tm);
            for e=1:eb
                for t=1:tm
                    ref_phase(e,t)=circ_mean(pres_phases(:,e,t));
                end
            end
            individual_mean_phase = NaN(ch,tm);
            relative_phase = [];
            for j=1:ch
                dif=squeeze(pres_phases(j,:,:))-ref_phase;
                for v=1:tm
                    individual_mean_phase(j,v)=circ_mean(dif(:,v));
                end
                relative_phase(:,:,j)=dif;
            end
    
            tal_cent = [[tal(:,1) - mean(tal(:,1))],[tal(:,2) - mean(tal(:,2))],[tal(:,3) - mean(tal(:,3))]];
            [coeff,score,~] = pca(tal_cent); % PCA to get the electrode plane
            score=score(:,1:2); % skip the third dimension, which is the least variance explained.
    
            %% Delegate circular linear regression model fitting to biowolf
            for el = 1:length(elecs)
                if sum(pres_cl(cl).pres_dist_elec(el,:) < el_dist) > 4
                    %for all tt in v and all i in eb
                    % Output to biowolf swarm file
                    functText = sprintf("%s('%d', '%d', '%s'); exit;", ...
                        txtFile.function, ct, el, data_path);
                    fprintf(fileID, '%s %scd(%s%s%s); %s%s\n', ...
                        txtFile.introText, '"', "'", txtFile.dir,  "'", functText, '"');
                end
            end
            % TO BE FILLED IN BY BIOWOLF FUNCTION
            % wave_info(ct).angle=ang;
            % wave_info(ct).CLcorrelation=cl_corr;
            % wave_info(ct).spatial_frequency=spatial_freq;
            % wave_info(ct).CLcorrelation=cl_corr;
            wave_info(ct).dims = size(pres_phases);
            wave_info(ct).relative_phase = relative_phase;
            wave_info(ct).score = score;
            wave_info(ct).sub = patID;
            wave_info(ct).sess = sessID;
            wave_info(ct).task_phase = 'presentation'; 
            %wave_info(ct).rt=[sess_struct.evs_pres.RT];
            wave_info(ct).cltal=tal;
            wave_info(ct).mfreq=mf;
            wave_info(ct).dist = pres_cl(cl).pres_dist_elec;
            wave_info(ct).elecs = el;
            ct = ct + 1; 
            clear relative_phase;
        end
    catch
        disp("No presentation data")
    end

    %% loop through recog clsuters
    try
        recog_cl = sess_clust(find(~cellfun(@isempty,{sess_clust.recog_freq})));
        recog_cl = recog_cl(~cellfun(@(C) isnan(C), {recog_cl.recog_freq}));
        for cl = 1:length(recog_cl)
            mf = recog_cl(cl).recog_freq;
            tal = recog_cl(cl).recog_tal;
            elecs = recog_cl(cl).recog_elecnum;
            channames = recog_cl(cl).recog_channame;
            recog_clust_eeg = sess_struct.raw_recog(ismember(sess_struct.chan_names, recog_cl(cl).recog_channame), :,:); %starts 7 seconds before cue onset, ends 7 seconds after
            recog_phases =  [];
            recog_amp = [];
            for el = 1:length(elecs)
                eegFilt = buttfilt(squeeze(recog_clust_eeg(el,:,:)),[mf*.85 mf/.85],resampleFreq,'bandpass',2);
                %             eegFilt = buttfilt(squeeze(cue_clust_eeg(el,:,:)),[2 32],resampleFreq,'bandpass',2); %generalized phase
    
                h = hilbert(eegFilt')'; %hilbert only goes across columns, unfortunately...
                h = h(:,1001:end-1000); %remove the buffering and extra time so it starts 1 second before cue and ends 2 seconds after cue off screen
                recog_phases(el,:,:) = angle(h); % chan x event x time
                recog_amp(el,:,:) = abs(h);
    
            end
    
            [ch,eb,tm]=size(recog_phases);
            ref_phase = NaN(eb,tm);
            for e=1:eb
                for t=1:tm
                    ref_phase(e,t)=circ_mean(recog_phases(:,e,t));
                end
            end
            individual_mean_phase = NaN(ch,tm);
            relative_phase = [];
            for j=1:ch
                dif=squeeze(recog_phases(j,:,:))-ref_phase;
                for v=1:tm
                    individual_mean_phase(j,v)=circ_mean(dif(:,v));
                end
                relative_phase(:,:,j)=dif;
            end
    
            tal_cent = [[tal(:,1) - mean(tal(:,1))],[tal(:,2) - mean(tal(:,2))],[tal(:,3) - mean(tal(:,3))]];
            [coeff,score,~] = pca(tal_cent); % PCA to get the electrode plane
            score=score(:,1:2); % skip the third dimention, which is the least variance explained.
           
             %% Delegate circular linear regression model fitting to biowolf
            for el = 1:length(elecs)
                if sum(recog_cl(cl).recog_dist_elec(el,:) < el_dist) > 4
                    %for all tt in v and all i in eb
                    % Output to biowolf swarm file
                    functText = sprintf("%s('%d', '%d', '%s'); exit;", ...
                        txtFile.function, ct, el, data_path);
                    fprintf(fileID, '%s %scd(%s%s%s); %s%s\n', ...
                        txtFile.introText, '"', "'", txtFile.dir,  "'", functText, '"');
                end
            end
            % TO BE FILLED IN BY BIOWOLF FUNCTION
            % wave_info(ct).angle=ang;
            % wave_info(ct).CLcorrelation=cl_corr;
            % wave_info(ct).spatial_frequency=spatial_freq;
            % wave_info(ct).CLcorrelation=cl_corr;
            wave_info(ct).dims = size(recog_phases);
            wave_info(ct).relative_phase = relative_phase;
            wave_info(ct).score = score;
            wave_info(ct).sub = patID;
            wave_info(ct).sess = sessID;
            wave_info(ct).task_phase = 'recog'; 
            %wave_info(ct).rt=[sess_struct.evs_recog.RT];
            wave_info(ct).cltal=tal;
            wave_info(ct).mfreq=mf;
            wave_info(ct).dist = recog_cl(cl).recog_dist_elec;
            wave_info(ct).elecs = el;
            ct = ct + 1; 
            clear relative_phase;
        end
    catch
        disp("No recognition clusters");
    end

    %% loop through retrieval clusters
    try
        fr_cl = sess_clust(find(~cellfun(@isempty,{sess_clust.fr_freq})));
        fr_cl = fr_cl(~cellfun(@(C) isnan(C), {fr_cl.fr_freq}));
        for cl = 1:length(fr_cl)
            mf = fr_cl(cl).fr_freq;
            tal = fr_cl(cl).fr_tal;
            elecs = fr_cl(cl).fr_elecnum;
            channames = fr_cl(cl).fr_channame;
            fr_clust_eeg = sess_struct.raw_fr(ismember(sess_struct.chan_names, fr_cl(cl).fr_channame), :,:);
            fr_phases =  [];
            fr_amp = [];
            for el = 1:length(elecs)
                eegFilt = buttfilt(squeeze(fr_clust_eeg(el,:,:)),[mf*.85 mf/.85],resampleFreq,'bandpass',2);
                %             eegFilt = buttfilt(squeeze(ret_clust_eeg(el,:,:)),[2 32],resampleFreq,'bandpass',2); %generalized phase
    
                h = hilbert(eegFilt')'; %hilbert only goes across columns, unfortunately...
                h = h(:,1001:end-1500); %remove the buffering and extra time so it starts 1 second before cue and ends 2 seconds after cue off screen
                fr_phases(el,:,:) = angle(h); % chan x event x time
                fr_amp(el,:,:) = abs(h);
    
            end
    
            [ch,eb,tm]=size(fr_phases);
            ref_phase = NaN(eb,tm);
            for e=1:eb
                for t=1:tm
                    ref_phase(e,t)=circ_mean(fr_phases(:,e,t));
                end
            end
            individual_mean_phase = NaN(ch,tm);
            relative_phase = [];
            for j=1:ch
                dif=squeeze(fr_phases(j,:,:))-ref_phase;
                for v=1:tm
                    individual_mean_phase(j,v)=circ_mean(dif(:,v));
                end
                relative_phase(:,:,j)=dif;
            end
    
            tal_cent = [[tal(:,1) - mean(tal(:,1))],[tal(:,2) - mean(tal(:,2))],[tal(:,3) - mean(tal(:,3))]];
            [coeff,score,~] = pca(tal_cent); % PCA to get the electrode plane
            score=score(:,1:2); % skip the third dimention, which is the least variance explained.
    
             % Delegate circular linear regression model fitting to biowolf
            for el = 1:length(elecs)
                if sum(fr_cl(cl).fr_dist_elec(el,:) < el_dist) > 4
                    %for all tt in v and all i in eb
                    % Output to biowolf swarm file
                    functText = sprintf("%s('%d', '%d', '%s'); exit;", ...
                        txtFile.function, ct, el, data_path);
                    fprintf(fileID, '%s %scd(%s%s%s); %s%s\n', ...
                        txtFile.introText, '"', "'", txtFile.dir,  "'", functText, '"');
                end
            end
            %TO BE FILLED IN BY BIOWOLF FUNCTION
            %wave_info(ct).angle=ang;
            %wave_info(ct).CLcorrelation=cl_corr;
            %wave_info(ct).spatial_frequency=spatial_freq;
            %wave_info(ct).CLcorrelation=cl_corr;
            wave_info(ct).dims = size(fr_phases);
            wave_info(ct).relative_phase = relative_phase;
            wave_info(ct).score = score;
            wave_info(ct).sub = patID;
            wave_info(ct).sess = sessID;
            wave_info(ct).task_phase = 'free recall'; 
            wave_info(ct).rt=[sess_struct.evs_fr.RT];
            wave_info(ct).cltal=tal;
            wave_info(ct).mfreq=mf;
            wave_info(ct).dist = fr_cl(cl).fr_dist_elec;
            wave_info(ct).elecs = el;
            ct = ct + 1; 
            clear relative_phase;
        end
    catch
        disp("No free recall data")
    end
%% Save output
    if exist('wave_info','var')
        out_dir = fullfile(folder_name, 'Processing');
        if ~exist(out_dir, 'dir') % make folder if doesn't exist already
            mkdir(out_dir);
        end
        save(fullfile(out_dir,"allwave.mat"), 'wave_info', '-v7.3');
    end
    clear fr_cl; clear pres_cl; clear fr_cl; clear wave_info;
end
