%% 4. Find all trials with clusters in specified reigion and frequency, plot desired figure (type)

clear; close all;
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/(7) Analysis/functions'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/circstat-matlab-master'));

load("/Volumes/Rahil_FRNU/Scripts/all_pres.mat");


%% Controls for plotting
freq_high = 12;
freq_low = 8;
% types: direction_compare | pres_compare | recog_compare | prevalence | angle_plot
type = 'prevalence';

%% SPLIT INTO SEPERATE STRUCTS BASED ON LOCATION

evs = all_pres;

a=0;
b=0;
c=0;
d=0;

for i = 1:numel(evs)

    if evs(i).mfreq < freq_high && evs(i).mfreq > freq_low
            significance = evs(i).loc_percent > 40;
            temp_struct = all_pres(index);
            temp_struct.cltal = evs(i).cltal;
            %if (sum(evs(i).cltal(:,1) > 0)/length(evs(i).cltal(:,1))) > .5 % if over half of electrodes are on right side

                if any((strcmp(evs(i).loc, 'anteriortemporal') & significance))
                    atl(a+1) = temp_struct;
                    a=a+1;
                end
                if any((strcmp(evs(i).loc, 'frontal') & significance))
                    fl(b+1) = temp_struct;
                    b=b+1;
                end
                %             if any((strcmp(evs(i).loc, 'parietal') & significance))
                %                 pl(c+1) = temp_struct;
                %                 c=c+1;
                %             end
                %             if any((strcmp(evs(i).loc, 'posteriortemporal') & significance))
                %                 ptl(d+1) = temp_struct;
                %                 d=d+1;
                %             end
           %end
    end
end
f=figure();
%% Anterior Temporal Lobe
if exist("atl", 'var')
    
    cued = [atl.cued_avg];
    n_cued = [atl.n_cued_avg];
    dif = [atl.dif];

    switch type
        case 'direction_compare'
            % If pres, compares angle while word is on
            % If recog, compares angle from word onset to RT
            % If fr, undecided (I think whole epoch)
            subplot(2,2,1);
            polarhistogram(cued_pres,30,'FaceColor','green','FaceAlpha',.3); hold on;
            polarhistogram(n_cued_pres,30, 'FaceColor','red','FaceAlpha',.3);
            legend(['pres cued, ', num2str(length(cued_pres))], ['pres not cued, ', num2str(length(n_cued_pres))])
            legend boxoff
            title('Cued vs Not Cued Within Trial ATL')

            subplot(2,2,2);
            polarhistogram(dif,45,'FaceColor','blue','FaceAlpha',.3);
            legend(['dif, ', num2str(length(cued_dif))])
            legend boxoff
            title('[Cued - Not Cued] Within Session Direction ATL')

        case 'prevalence'
            cued_num = [atl.cued_prev_num];
            uncued_num = [atl.n_cued_prev_num];

            cued_len = [atl.cued_prev_len];
            uncued_len = [atl.n_cued_prev_len];

            cued_num = normalize(cued_num,1);
            uncued_num = normalize(uncued_num,1);

            cued_len = normalize(cued_len,1);
            uncued_len = normalize(uncued_len,1);

            t_series = atl.prev_t;
            t_series = t_series +250;

            subplot(221);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_num')), squeeze(nanstd(cued_num'))./sqrt(length(atl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_num')), squeeze(nanstd(uncued_num'))./sqrt(length(atl)), Lineprops); hold on;

            [h_1, pval] = stat_multcomp(cued_num - uncued_num);

            h_1=double(h_1');
            h_1(h_1==0) = NaN;
            plot(t_series, h_1*0.9, 'color','m', 'LineWidth', 10);

            legend(['cued, ', num2str(length(atl))], ['no cue, ', num2str(length(atl))], 'Multiple Comparisons Test', 'Location', 'southeast', 'AutoUpdate','off')
            title('Number of Waves in ATL')
            xlabel('Time (ms)');
            ylabel('Normalized Number of Waves');
            
            bands = [1000 1500; 2000 3000];
            xp = [bands fliplr(bands)];
            yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';
            for k = 1:size(bands,1)
                patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.1, 'EdgeColor','black', 'LineWidth', 2)
            end

            subplot(222);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_len')), squeeze(nanstd(cued_len'))./sqrt(length(atl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_len')), squeeze(nanstd(uncued_len'))./sqrt(length(atl)), Lineprops); hold on;

            [h_1, pval] = stat_multcomp(cued_len - uncued_len);

            h_1=double(h_1');
            h_1(h_1==0) = NaN;
            plot(t_series, h_1*0.9, 'color','m', 'LineWidth', 10);
            top = max(nanmean(cued_len'));
            bottom = min(nanmean(uncued_len'));

            legend(['cued, ', num2str(length(atl))], ['no cue, ', num2str(length(atl))], 'Multiple Comparisons Test','Location', 'southwest', 'AutoUpdate','off')
            title('Average Duration of Waves in ATL');
            xlabel('Time (ms)');
            ylabel('Normalized Duration of Waves');
            
            bands = [1000 1500; 2000 3000];
            xp = [bands fliplr(bands)];
            yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';
            for k = 1:size(bands,1)
                patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.1, 'EdgeColor','black', 'LineWidth', 2)
            end

        case 'angle_plot'

            cued = [atl.cued_signal];
            cued = nanmean(cued,2);
            n_cued = [atl.n_cued_signal];
            n_cued = nanmean(n_cued,2);

            % define time series
            t = (1:3000)/1000;
            t=t';
            % works to plot multiple sessions, one session, or avg
            t_cued= repmat(t,1,size(cued,2));
            t_n_cued = repmat(t,1,size(n_cued,2));
            
            % data for cue off and word presentation markers
            circle = (1:3000)/(1000/pi);
            t_cue_on = ones(1,3000);
            t_cue_off = ones(1,3000)*1.5;
            t_word_on = ones(1,3000)*2;
            
            subplot(1,2,1);

            polarplot(cued,t_cued(:,1),'g.', 'LineWidth', 2); hold on;
            polarplot(n_cued,t_n_cued(:,1),'r.', 'LineWidth', 2); hold on;
            polarplot(circle,t_cue_off,'b--','LineWidth', 2); hold on;
            polarplot(circle,t_word_on,'m--','LineWidth', 2);  hold on;
            polarplot(circle,t_cue_on,'c--', 'LineWidth', 2); 
            legend('cue', 'baseline (not cued)','cue turns on','cue turns off', 'word presentation')
            legend boxoff
            title('Average Angle Across Trials in 0 to 3 Seconds During Presentation in ATL')
    end
end

%% Frontal Lobe
if exist("fl", 'var')
    cued = [fl.cued_avg];
    n_cued = [fl.n_cued_avg];
    dif = [fl.dif];

    switch type
        case 'direction_compare'
            % If pres, compares angle while word is on
            % If recog, compares angle from word onset to RT
            % If fr, undecided (I think whole epoch)
            subplot(2,2,1);
            polarhistogram(cued_pres,30,'FaceColor','green','FaceAlpha',.3); hold on;
            polarhistogram(n_cued_pres,30, 'FaceColor','red','FaceAlpha',.3);
            legend(['pres cued, ', num2str(length(cued_pres))], ['pres not cued, ', num2str(length(n_cued_pres))])
            legend boxoff
            title('Cued vs Not Cued Within Trial FL')

            subplot(2,2,2);
            polarhistogram(dif,45,'FaceColor','blue','FaceAlpha',.3);
            legend(['dif, ', num2str(length(cued_dif))])
            legend boxoff
            title('[Cued - Not Cued] Within Session Direction FL')

        case 'prevalence'
            cued_num = [fl.cued_prev_num];
            uncued_num = [fl.n_cued_prev_num];

            cued_len = [fl.cued_prev_len];
            uncued_len = [fl.n_cued_prev_len];

            cued_num = normalize(cued_num,1);
            uncued_num = normalize(uncued_num,1);

            cued_len = normalize(cued_len,1);
            uncued_len = normalize(uncued_len,1);

            t_series = fl.prev_t;
            t_series = t_series +250;

            subplot(223);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_num')), squeeze(nanstd(cued_num'))./sqrt(length(fl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_num')), squeeze(nanstd(uncued_num'))./sqrt(length(fl)), Lineprops); hold on;

            [h_1, ~] = stat_multcomp(cued_num - uncued_num);

            h_1=double(h_1');
            h_1(h_1==0) = NaN;
            plot(t_series, h_1*1.5, 'color','m', 'LineWidth', 10);
            top = max(nanmean(cued_num'));
            bottom = min(nanmean(uncued_num'));

            legend(['cued, ', num2str(length(fl))], ['no cue, ', num2str(length(fl))], 'Multiple Comparisons Test', 'Location', 'northwest', 'AutoUpdate', 'off')
            title('Number of Waves in FL')
            xlabel('Time (ms)');
            ylabel('Normalized Number of Waves');
            
            bands = [1000 1500; 2000 3000];
            xp = [bands fliplr(bands)];
            yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';
            for k = 1:size(bands,1)
                patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.1, 'EdgeColor','black', 'LineWidth', 2)
            end


            subplot(224);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_len')), squeeze(nanstd(cued_len'))./sqrt(length(fl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_len')), squeeze(nanstd(uncued_len'))./sqrt(length(fl)), Lineprops); hold on;

            [h_1, pval] = stat_multcomp(cued_len - uncued_len);

            h_1=double(h_1');
            h_1(h_1==0) = NaN;
            plot(t_series, h_1*2, 'color','m', 'LineWidth', 10);
            top = max(nanmean(cued_len'));
            bottom = min(nanmean(uncued_len'));

            legend(['cued, ', num2str(length(fl))], ['no cue, ', num2str(length(fl))], 'Multiple Comparisons Test' ,'Location', 'northwest', 'AutoUpdate', 'off')
            title('Average Duration of Waves in FL');
            xlabel('Time (ms)');
            ylabel('Normalized Duration of Waves');
            
            bands = [1000 1500; 2000 3000];
            xp = [bands fliplr(bands)];
            yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';
            for k = 1:size(bands,1)
                patch(xp(k,:), yp(k,:), [1 1 1]*0.25, 'FaceAlpha',0.1, 'EdgeColor','black', 'LineWidth', 2)
            end

        case 'angle_plot'

            cued = [fl.cued_signal];
            cued = nanmean(cued,2);
            n_cued = [fl.n_cued_signal];
            n_cued = nanmean(n_cued,2);

            % define time series
            t = (1:3000)/1000;
            t=t';
            % works to plot multiple sessions, one session, or avg
            t_cued= repmat(t,1,size(cued,2));
            t_n_cued = repmat(t,1,size(n_cued,2));
            
            % data for cue off and word presentation markers
            circle = (1:3000)/(1000/pi);
            t_cue_on = ones(1,3000);
            t_cue_off = ones(1,3000)*1.5;
            t_word_on = ones(1,3000)*2;
            
            subplot(1,2,2);

            polarplot(cued,t_cued(:,1),'g.', 'LineWidth', 2); hold on;
            polarplot(n_cued,t_n_cued(:,1),'r.', 'LineWidth', 2); hold on;
            polarplot(circle,t_cue_off,'b--','LineWidth', 2); hold on;
            polarplot(circle,t_word_on,'m--','LineWidth', 2);  hold on;
            polarplot(circle,t_cue_on,'c--', 'LineWidth', 2); 
            legend('cue', 'baseline (not cued)','cue turns on','cue turns off', 'word presentation')
            legend boxoff
            title('0 to 3 Seconds During Presentation in FL')

    end
end

%% Parietal Lobe
if exist("pl", 'var')
    cued = [pl.cued_avg];
    n_cued = [pl.n_cued_avg];
    dif = [pl.dif];

    switch type
        case 'direction_compare'
            % If pres, compares angle while word is on
            % If recog, compares angle from word onset to RT
            % If fr, undecided (I think whole epoch)
            subplot(2,2,1);
            polarhistogram(cued_pres,30,'FaceColor','green','FaceAlpha',.3); hold on;
            polarhistogram(n_cued_pres,30, 'FaceColor','red','FaceAlpha',.3);
            legend(['pres cued, ', num2str(length(cued_pres))], ['pres not cued, ', num2str(length(n_cued_pres))])
            legend boxoff
            title('Cued vs Not Cued Within Trial PL')

            subplot(2,2,2);
            polarhistogram(dif,45,'FaceColor','blue','FaceAlpha',.3);
            legend(['dif, ', num2str(length(cued_dif))])
            legend boxoff
            title('[Cued - Not Cued] Within Session Direction PL')

        case 'prevalence'
            cued_num = [pl.cued_prev_num];
            uncued_num = [pl.n_cued_prev_num];

            cued_len = [pl.cued_prev_len];
            uncued_len = [pl.n_cued_prev_len];

            cued_num = normalize(cued_num,1);
            uncued_num = normalize(uncued_num,1);

            cued_len = normalize(cued_len,1);
            uncued_len = normalize(uncued_len,1);

            t_series = pl.prev_t;

            subplot(3,4,9);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_num')), squeeze(nanstd(cued_num'))./sqrt(length(pl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_num')), squeeze(nanstd(uncued_num'))./sqrt(length(pl)), Lineprops);
            top = max(nanmean(cued_num'));
            bottom = min(nanmean(uncued_num'));

            legend(['cued, ', num2str(length(pl))], ['no cue, ', num2str(length(pl))])
            title('Number of Waves in PL (sliding window)')
            xlabel('Time (ms)');
            ylabel('Normalized Number of Waves Averaged Across Trials');
            rectangle('Position', [2000 bottom-0.5 1000 top-bottom+0.5], 'FaceColor', [0 0 0.0 0.05])

            subplot(3,4,10);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_len')), squeeze(nanstd(cued_len'))./sqrt(length(pl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_len')), squeeze(nanstd(uncued_len'))./sqrt(length(pl)), Lineprops);
            top = max(nanmean(cued_len'));
            bottom = min(nanmean(uncued_len'));

            legend(['cued, ', num2str(length(pl))], ['no cue, ', num2str(length(pl))])
            title('Average Duration of Waves in PL (sliding window)');
            xlabel('Time (ms)');
            ylabel('Normalized Duration of Waves Averaged Across Trials');
            rectangle('Position', [2000 bottom-0.5 1000 top-bottom+0.5], 'FaceColor', [0 0 0.0 0.05])
        case 'angle_plot'

            cued = [pl.cued_signal];
            cued = nanmean(cued,2);
            n_cued = [pl.n_cued_signal];
            n_cued = nanmean(n_cued,2);

            % define time series
            t = (1:3000)/1000;
            t=t';
            % works to plot multiple sessions, one session, or avg
            t_cued= repmat(t,1,size(cued,2));
            t_n_cued = repmat(t,1,size(n_cued,2));
            
            % data for cue off and word presentation markers
            circle = (1:3000)/(1000/pi);
            t_cue_on = ones(1,3000);
            t_cue_off = ones(1,3000)*1.5;
            t_word_on = ones(1,3000)*2;
            
            subplot(1,1,1);

            polarplot(cued,t_cued(:,1),'g.', 'LineWidth', 2); hold on;
            polarplot(n_cued,t_n_cued(:,1),'r.', 'LineWidth', 2); hold on;
            polarplot(circle,t_cue_off,'b--','LineWidth', 2); hold on;
            polarplot(circle,t_word_on,'m--','LineWidth', 2);  hold on;
            polarplot(circle,t_cue_on,'c--', 'LineWidth', 2); 
            legend('cue', 'baseline (not cued)','cue turns on','cue turns off', 'word presentation')
            legend boxoff
            title('0 to 3 Seconds During Presentation in PL')
    end
end
%% Posterior Temporal Lobe
if exist("ptl", 'var')
    cued = [ptl.cued_avg];
    n_cued = [ptl.n_cued_avg];
    dif = [ptl.dif];

    switch type
        case 'direction_compare'
            % If pres, compares angle while word is on
            % If recog, compares angle from word onset to RT
            % If fr, undecided (I think whole epoch)
            subplot(2,2,1);
            polarhistogram(cued_pres,30,'FaceColor','green','FaceAlpha',.3); hold on;
            polarhistogram(n_cued_pres,30, 'FaceColor','red','FaceAlpha',.3);
            legend(['pres cued, ', num2str(length(cued_pres))], ['pres not cued, ', num2str(length(n_cued_pres))])
            legend boxoff
            title('Cued vs Not Cued Within Trial PTL')

            subplot(2,2,2);
            polarhistogram(dif,45,'FaceColor','blue','FaceAlpha',.3);
            legend(['dif, ', num2str(length(cued_dif))])
            legend boxoff
            title('[Cued - Not Cued] Within Session Direction PTL')
        case 'prevalence'
            cued_num = [ptl.cued_prev_num];
            uncued_num = [ptl.n_cued_prev_num];

            cued_len = [ptl.cued_prev_len];
            uncued_len = [ptl.n_cued_prev_len];

            cued_num = normalize(cued_num,1);
            uncued_num = normalize(uncued_num,1);

            cued_len = normalize(cued_len,1);
            uncued_len = normalize(uncued_len,1);

            t_series = ptl.prev_t;

            subplot(3,4,11);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_num')), squeeze(nanstd(cued_num'))./sqrt(length(ptl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_num')), squeeze(nanstd(uncued_num'))./sqrt(length(ptl)), Lineprops);
            top = max(nanmean(cued_num'));
            bottom = min(nanmean(uncued_num'));

            legend(['cued, ', num2str(length(ptl))], ['no cue, ', num2str(length(ptl))])
            title('Number of Waves in PTL (sliding window)')
            xlabel('Time (ms)');
            ylabel('Normalized Number of Waves Averaged Across Trials');
            rectangle('Position', [2000 bottom-0.5 1000 top-bottom+0.5], 'FaceColor', [0 0 0.0 0.05])

            subplot(3,4,12);
            Lineprops.col={'g'};
            mseb(t_series, squeeze(nanmean(cued_len')), squeeze(nanstd(cued_len'))./sqrt(length(ptl)), Lineprops); hold on;
            Lineprops.col={'r'};
            mseb(t_series, squeeze(nanmean(uncued_len')), squeeze(nanstd(uncued_len'))./sqrt(length(ptl)), Lineprops);
            top = max(nanmean(cued_len'));
            bottom = min(nanmean(uncued_len'));

            legend(['cued, ', num2str(length(ptl))], ['no cue, ', num2str(length(ptl))])
            title('Average Duration of Waves in PTL');
            xlabel('Time (ms)');
            ylabel('Normalized Duration of Waves Averaged Across Trials');
            rectangle('Position', [2000 bottom-0.5 1000 top-bottom+0.5], 'FaceColor', [0 0 0.0 0.05])
        case 'angle_plot'

            cued = [ptl.cued_signal];
            cued = nanmean(cued,2);
            n_cued = [ptl.n_cued_signal];
            n_cued = nanmean(n_cued,2);

            % define time series
            t = (1:3000)/1000;
            t=t';
            % works to plot multiple sessions, one session, or avg
            t_cued= repmat(t,1,size(cued,2));
            t_n_cued = repmat(t,1,size(n_cued,2));
            
            % data for cue off and word presentation markers
            circle = (1:3000)/(1000/pi);
            t_cue_on = ones(1,3000);
            t_cue_off = ones(1,3000)*1.5;
            t_word_on = ones(1,3000)*2;
            
            subplot(1,1,1);

            polarplot(cued,t_cued(:,1),'g.', 'LineWidth', 2); hold on;
            polarplot(n_cued,t_n_cued(:,1),'r.', 'LineWidth', 2); hold on;
            polarplot(circle,t_cue_off,'b--','LineWidth', 2); hold on;
            polarplot(circle,t_word_on,'m--','LineWidth', 2);  hold on;
            polarplot(circle,t_cue_on,'c--', 'LineWidth', 2); 
            legend('cue', 'baseline (not cued)','cue turns on','cue turns off', 'word presentation')
            legend boxoff
            title('0 to 3 Seconds During Presentation in PTL')
    end
end

%% Save figure
set(f,'position',[0,0,3000,1000]);


switch type
    case 'direction_compare'
        sgtitle('Cued vs Not Cued Average Direction by Trial in Presentation and Recognition');
        fontsize(f, 26, 'points');
        saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high),'/direction_compare.png')));
    case 'pres_compare'
        sgtitle('Cued vs Not Cued Average Presentation Direction by Trial');
        fontsize(f, 26, 'points');
        f.WindowState = 'maximized';
        set(f,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        print -dpdf -vector epsFig
        %saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high),'/pres_direction_compare.pdf')));
    case 'recog_compare'
        sgtitle('Cued vs Not Cued Average Recognition Direction by Trial');
        fontsize(f, 26, 'points');
        saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high),'/recog_direction_compare.png')));
    case 'consistency'
        sgtitle('Consistency of Electrodes Within Cluster by Trial');
        fontsize(f, 26, 'points');
        saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high),'/correlation.png')));
    case 'cue_locked'
        sgtitle('Cued vs Not Cued Average Direction by Trial Cue Locked');
        fontsize(f, 26, 'points');
        saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high),'/cue_locked_direction_compare.png')));
    case 'prevalence'
        sgtitle(['Trends in Wave Prevalence Averaged Across Clusters Between ', num2str(freq_low), ' and ', num2str(freq_high), ' Hz']);
        fontsize(f, 26, 'points');
        saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Prevelance/', append(num2str(freq_low), "_to_",num2str(freq_high),'/prevalence_trends.png')));
    case 'angle_plot'
        %sgtitle(['Angle from 0 to 3 seconds in Clusters Between ', num2str(freq_low), ' and ', num2str(freq_high), ' Hz']);
        fontsize(f, 26, 'points');
        saveas(gcf, fullfile('/Volumes/Rahil_FRNU/Figures/Direction/', append(num2str(freq_low), "_to_",num2str(freq_high),'/angle_plot.png')));
end

%% Multiple Comparisions Shuffling

%
% [h_1, pval] = stat_multcomp(cued - uncued);
%
% h_1=double(h_1');
% h_1(h_1==0) = NaN;
% plot([1:5000], h_1*.18, 'color','r', 'LineWidth', 2);

