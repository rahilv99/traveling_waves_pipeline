%% 4. Find all trials with clusters in specified reigion and frequency, plot desired figure (type)

clear; close all;
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/(7) Analysis/functions'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/circstat-matlab-master'));

load("/Volumes/Rahil_FRNU/Scripts/all_pres.mat");

%% Controls for plotting
freq_high = 12;
freq_low = 8;
normalized = true;
% types: wave_prevalance | direction_clust | percent_pref_dir | phase_velocity
type = 'wave_prevalance';
ev_type = 'pres';

subj = 'all';
sess = 'all_';
save = true;
% SPLIT INTO SEPERATE STRUCTS BASED ON LOCATION
switch ev_type
    case 'pres'
        evs = all_pres;
        t = 5000;
    case 'recog'
        evs = all_recog;
        t= 6000;
end

a=0;
b=0;
c=0;
d=0;
e=0;
f=0;
valid = false;

s = split(sess, '_');
id = append(subj,' ',s(1),' ',s(2));

for i = 1:numel(evs)
    if evs(i).mfreq < freq_high && evs(i).mfreq > freq_low % && any(strcmp(evs(i).sub, subj) & strcmp(evs(i).sess,sess))
        %if (sum(evs(i).cltal(:,1) > 0)/length(evs(i).cltal(:,1))) > .5 % if over half of electrodes are on right side
            valid = true;
            significance = evs(i).loc_percent > 40;
            if any((strcmp(evs(i).loc, 'anteriortemporal') & significance))
                atl(a+1) = evs(i);
                a=a+1;
            end
            if any((strcmp(evs(i).loc, 'frontal') & significance))
                fl(c+1) = evs(i);
                c=c+1;
            end
    %         if any((strcmp(evs(i).loc, 'parietal') & significance))
    %             pl(e+1) = evs(i);
    %             e=e+1;
    %         end
    %         if any((strcmp(evs(i).loc, 'posteriortemporal') & significance))
    %             ptl(f+1) = evs(i);
    %             f=f+1;
    %         end
        %end
    end
end

if valid
    figure();
    
    %% Anterior Temporal Lobe
    if exist('atl','var')
        % Initialize a cell array to store the values of 'loc'
        
        cued = nan(size(atl,2),t);
        uncued = nan(size(atl,2),t);
        sz = 0;
        
        
        % Extract 'loc' field values using a loop
        for i = 1:numel(atl)
            switch type
                case 'wave_prevalance'
                    cued(i,:) = nanmean(atl(i).wave_logic(atl(i).cued,:),1);
                    uncued(i,:) = nanmean(atl(i).wave_logic(atl(i).n_cued,:),1);
                case 'direction_clust'
                    cued(i,:) = atl(i).dc_inclust_cued;
                    uncued(i,:) = atl(i).dc_inclust_n_cued;
                case 'percent_pref_dir'
                    cued(i,:) = atl(i).percent_pref_dir_cued;
                    uncued(i,:) = atl(i).percent_pref_dir_n_cued;
                case 'phase_velocity'
                    cued(i,:) = atl(i).phv_cued;
                    uncued(i,:) = atl(i).phv_n_cued;
            end
            sz = sz+1;
        end
        subplot(1,2,2);
        cued = cued';
        uncued = uncued';
        %sliding window avg
        window_size = 20;
        step_size = 1;


        cued_smooth = nan([((size(cued,1)-window_size)/step_size),size(cued,2)]);
        uncued_smooth = nan([((size(uncued,1)-window_size)/step_size),size(cued,2)]);

        for i = 1:(size(cued,1)-window_size)/step_size
            t_i = i*step_size;
            t_f = t_i + window_size;

            for tr = 1:size(cued,2)
                samp = cued(t_i:t_f,tr);
                cued_smooth(i,tr) = nanmean(samp);
            end
            for tr = 1:size(uncued,2)
                samp = uncued(t_i:t_f,tr);
                uncued_smooth(i,tr) = nanmean(samp);
            end
        end
        
        if sz>1
            if normalized
                cued_smooth = normalize(cued_smooth,1);
                uncued_smooth = normalize(uncued_smooth,1);
            end
            Lineprops.col={'g'};
            mseb(1+window_size/2:step_size:(t-window_size/2), squeeze(nanmean(cued_smooth')), squeeze(nanstd(cued_smooth'))./sqrt(sz), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(1+window_size/2:step_size:(t-window_size/2), squeeze(nanmean(uncued_smooth')), squeeze(nanstd(uncued_smooth'))./sqrt(sz), Lineprops); hold on;
            top = max(nanmean(cued_smooth'));
            bottom = min(nanmean(uncued_smooth'));
        else
            plot(1:t, cued, 'g'); hold on;
            plot(1:t, uncued, 'r');
            top = max(max(cued_smooth));
            bottom = min(min(uncued_smooth));
        end

        [h_1, pval] = stat_multcomp(cued_smooth - uncued_smooth);

        h_1=double(h_1');
        h_1(h_1==0) = NaN;
        plot(1+window_size/2:step_size:(t-window_size/2), h_1*(bottom-0.5), 'color','m', 'LineWidth', 10);


        switch ev_type
            case 'pres'
                bands = [1000 1500; 2000 3000];                                                      
                xp = [bands fliplr(bands)];                                                 
                yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';              
                for k = 1:size(bands,1)                                               
                    patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.1, 'EdgeColor','black', 'LineWidth', 2)
                end

            case 'recog'
                rectangle('Position', [1000 0 1000 1], 'FaceColor', [0 0 0.0 0.05])
        end

        title("Anterior Temporal Lobe")
        switch type
            case 'wave_prevalance'
                ylabel("normalized number of trials with waves")
            case 'direction_clust'
                ylabel("rbar")
            case 'percent_pref_dir'
                ylabel("normalized number of waves in preferred direction")
            case 'phase_velocity'
                ylabel("normalized phase velocity")
        end
        xlabel("ms")
        legend(['cued, ', num2str(sz)], ['no cue, ', num2str(sz)],"multiple comparisons test")
    end
    
    %% Frontal Lobe
    if exist('fl','var')
        
        % Initialize a cell array to store the values of 'loc'
        cued = nan(size(fl,2),t);
        uncued = nan(size(fl,2),t);
        sz = 0;
        
        % Extract 'loc' field values using a loop
        for i = 1:numel(fl)
            switch type
                case 'wave_prevalance'
                    cued(i,:) = nanmean(fl(i).wave_logic(fl(i).cued,:),1);
                    uncued(i,:) = nanmean(fl(i).wave_logic(fl(i).n_cued,:),1);
                case 'direction_clust'
                    cued(i,:) = fl(i).dc_inclust_cued;
                    uncued(i,:) = fl(i).dc_inclust_n_cued;
                case 'percent_pref_dir'
                    cued(i,:) = fl(i).percent_pref_dir_cued;
                    uncued(i,:) = fl(i).percent_pref_dir_n_cued;
                case 'phase_velocity'
                    cued(i,:) = atl(i).phv_cued;
                    uncued(i,:) = atl(i).phv_n_cued;
            end
            sz = sz+1;
        end
        subplot(1,2,1);
        cued = cued';
        uncued = uncued';
        
        %sliding window avg
        window_size = 20;
        step_size = 1;

        cued_smooth = nan([((size(cued,1)-window_size)/step_size),size(cued,2)]);
        uncued_smooth = nan([((size(uncued,1)-window_size)/step_size),size(cued,2)]);

        for i = 1:(size(cued,1)-window_size)/step_size
            t_i = i*step_size;
            t_f = t_i + window_size;

            for tr = 1:size(cued,2)
                samp = cued(t_i:t_f,tr);
                cued_smooth(i,tr) = nanmean(samp);
            end
            for tr = 1:size(uncued,2)
                samp = uncued(t_i:t_f,tr);
                uncued_smooth(i,tr) = nanmean(samp);
            end
        end
        
        if sz>1
            if normalized
                cued_smooth = normalize(cued_smooth,1);
                uncued_smooth = normalize(uncued_smooth,1);
            end
            Lineprops.col={'g'};
            mseb(1+window_size/2:step_size:(t-window_size/2), squeeze(nanmean(cued_smooth')), squeeze(nanstd(cued_smooth'))./sqrt(sz), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(1+window_size/2:step_size:(t-window_size/2), squeeze(nanmean(uncued_smooth')), squeeze(nanstd(uncued_smooth'))./sqrt(sz), Lineprops); hold on;
            top = max(nanmean(cued_smooth'));
            bottom = min(nanmean(uncued_smooth'));
        else
            plot(1:t, cued, 'g'); hold on;
            plot(1:t, uncued, 'r');
            top = max(max(cued_smooth));
            bottom = min(min(uncued_smooth));
        end

        [h_1, pval] = stat_multcomp(cued_smooth - uncued_smooth);

        h_1=double(h_1');
        h_1(h_1==0) = NaN;
        plot(1+window_size/2:step_size:(t-window_size/2), h_1*(bottom-0.5), 'color','m', 'LineWidth', 10);

        switch ev_type
            case 'pres'
                bands = [1000 1500; 2000 3000];                                                      
                xp = [bands fliplr(bands)];                                                 
                yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';              
                for k = 1:size(bands,1)                                               
                    patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.1, 'EdgeColor','black', 'LineWidth', 2)
                end
            case 'recog'
                rectangle('Position', [1000 0 1000 1], 'FaceColor', [0 0 0.0 0.05])
        end
        
        title("Frontal Lobe")
        switch type
            case 'wave_prevalance'
                ylabel("normalized number of trials with waves")
            case 'direction_clust'
                ylabel("rbar")
            case 'percent_pref_dir'
                ylabel("normalized number of waves in preferred direction")
            case 'phase_velocity'
                ylabel("normalized phase velocity")
        end
        xlabel("ms")
        legend(['cued, ', num2str(sz)], ['no cue, ', num2str(sz)], "multiple comparisons test")
    end
  
    %% Parietal Lobe
    if exist('pl','var')
        % Initialize a cell array to store the values of 'loc'
        cued = nan(size(pl,2),t);
        uncued = nan(size(pl,2),t);
        sz = 0;
        
        % Extract 'loc' field values using a loop
        for i = 1:numel(pl)
            switch type
                case 'wave_prevalance'
                    cued(i,:) = nanmean(pl(i).wave_logic(pl(i).cued,:),1);
                    uncued(i,:) = nanmean(pl(i).wave_logic(pl(i).n_cued,:),1);
                case 'direction_clust'
                    cued(i,:) = pl(i).dc_inclust_cued;
                    uncued(i,:) = pl(i).dc_inclust_n_cued;
                case 'percent_pref_dir'
                    cued(i,:) = pl(i).percent_pref_dir_cued;
                    uncued(i,:) = pl(i).percent_pref_dir_n_cued;
                case 'phase_velocity'
                    cued(i,:) = atl(i).phv_cued;
                    uncued(i,:) = atl(i).phv_n_cued;
            end
            sz = sz+1;
        end
        subplot(2,3,5);
        
        cued = cued';
        uncued = uncued';
        
        if sz>1
            if normalized
                cued = normalize(cued,1);
                uncued = normalize(uncued,1);
            end
            Lineprops.col={'g'};
            mseb(1:t, squeeze(nanmean(cued')), squeeze(nanstd(cued'))./sqrt(sz), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(1:t, squeeze(nanmean(uncued')), squeeze(nanstd(uncued'))./sqrt(sz), Lineprops);
            top = max(nanmean(cued'));
            bottom = min(nanmean(uncued'));
        else
            plot(1:t, cued, 'g'); hold on;
            plot(1:t, uncued, 'r');
            top = max(max(cued));
            bottom = min(min(uncued));
        end

        [h_1, pval] = stat_multcomp(cued - uncued);

        h_1=double(h_1');
        h_1(h_1==0) = NaN;
        plot(1:t, h_1*(bottom-0.5), 'color','m', 'LineWidth', 4);
        
        switch ev_type
            case 'pres'
                rectangle('Position', [2000 bottom-0.5 1000 top-bottom+0.5], 'FaceColor', [0 0 0.0 0.05])
                xline(1000, 'LineWidth',1.5)
                xline(1500, 'LineWidth',1.5)
            case 'recog'
                rectangle('Position', [1000 0 1000 1], 'FaceColor', [0 0 0.0 0.05])
        end
        
        title("Parietal Lobe")
        switch type
            case 'wave_prevalance'
                ylabel("normalized number of trials with waves")
            case 'direction_clust'
                ylabel("rbar")
            case 'percent_pref_dir'
                ylabel("normalized number of waves in preferred direction")
        end
        xlabel("ms")
        legend(['cued, ', num2str(sz)], ['no cue, ', num2str(sz)])
    end
    
    %% Posterior Temporal Lobe
    if exist('ptl','var')
        % Initialize a cell array to store the values of 'loc'
        cued = nan(size(ptl,2),t);
        uncued = nan(size(ptl,2),t);
        sz = 0;
        
        % Extract 'loc' field values using a loop
        for i = 1:numel(ptl)
            switch type
                case 'wave_prevalance'
                    cued(i,:) = nanmean(ptl(i).wave_logic(ptl(i).cued,:),1);
                    uncued(i,:) = nanmean(ptl(i).wave_logic(ptl(i).n_cued,:),1);
                case 'direction_clust'
                    cued(i,:) = ptl(i).dc_inclust_cued;
                    uncued(i,:) = ptl(i).dc_inclust_n_cued;
                case 'percent_pref_dir'
                    cued(i,:) = ptl(i).percent_pref_dir_cued;
                    uncued(i,:) = ptl(i).percent_pref_dir_n_cued;
                case 'phase_velocity'
                    cued(i,:) = atl(i).phv_cued;
                    uncued(i,:) = atl(i).phv_n_cued;
            end
            sz = sz+1;
        end
        subplot(2,3,6);
        cued = cued';
        uncued = uncued';
        
        if sz>1
            if normalized
                cued = normalize(cued,1);
                uncued = normalize(uncued,1);
            end
            Lineprops.col={'g'};
            mseb(1:t, squeeze(nanmean(cued')), squeeze(nanstd(cued'))./sqrt(sz), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(1:t, squeeze(nanmean(uncued')), squeeze(nanstd(uncued'))./sqrt(sz), Lineprops);
            top = max(nanmean(cued'));
            bottom = min(nanmean(uncued'));
        else
            plot(1:t, cued, 'g'); hold on;
            plot(1:t, uncued, 'r');
            top = max(max(cued));
            bottom = min(min(uncued));
        end

        [h_1, pval] = stat_multcomp(cued - uncued);

        h_1=double(h_1');
        h_1(h_1==0) = NaN;
        plot(1:t, h_1*(bottom-0.5), 'color','m', 'LineWidth', 10);
        
        switch ev_type
            case 'pres'
                rectangle('Position', [2000 bottom-0.5 1000 top-bottom+0.5], 'FaceColor', [0 0 0.0 0.05])
                xline(1000, 'LineWidth',1.5)
                xline(1500, 'LineWidth',1.5)
            case 'recog'
                rectangle('Position', [1000 0 1000 1], 'FaceColor', [0 0 0.0 0.05])
        end
        
        title("Posterior Temporal Lobe")
        switch type
            case 'wave_prevalance'
                ylabel("normalized number of trials with waves")
            case 'direction_clust'
                ylabel("rbar")
            case 'percent_pref_dir'
                ylabel("normalized number of waves in preferred direction")
        end
        xlabel("ms")
        legend(['cued, ', num2str(sz)], ['no cue, ', num2str(sz)])
    end
    
    %% Save figure
    set(gcf,'position',[0,0,3000,1000]);
    
    switch type
        case 'wave_prevalance'
            sgtitle(['Average Wave Prevalance Across Sessions with Clusters Between ', num2str(freq_low), ' and ', num2str(freq_high), ' Hz']);
            fontsize(gcf, 28, 'points');
            if save
                saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Prevelance/', append(num2str(freq_low), "_to_",num2str(freq_high), "/wave_prevalance.png")));
            end
        case 'direction_clust'
            sgtitle(['Average Directional Consistency Across Sessions with Clusters Between ', num2str(freq_low), ' and ', num2str(freq_high), ' Hz']);
            fontsize(gcf, 28, 'points');
            if save
                saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high), "/direction_clust.png")));
            end
        case 'percent_pref_dir'
            sgtitle(['Average Percent Preferred Direction Across Sessions with Clusters Between ', num2str(freq_low), ' and ', num2str(freq_high), ' Hz']);
            fontsize(gcf, 28, 'points');
            if save
                saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high), "/percent_pref_dir.png")));
            end
        case 'phase_velocity'
            sgtitle(['Average Phase Velocity Across Sessions with Clusters Between ', num2str(freq_low), ' and ', num2str(freq_high), ' Hz']);
            fontsize(gcf, 28, 'points');
            if save
                saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high), "/phase_velocity.png")));
            end
    end
    
else
    disp("Requested data not found. Try different parameters.");
end

%% Multiple Comparisions Shuffling
