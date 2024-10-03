%% Plot angle of all electrodes in ATL and FL at max wave prevalence timing (can add parietal) 

clear; close all;
% add paths
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/ZaghloulCodebase'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/circstat-matlab-master'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/(7) Analysis/functions'));


load("/Volumes/Rahil_FRNU/Scripts/all_pres.mat"); % use wave prevalence from wave_logic (in loc_filter), ang
evs = all_pres;
root = '/Volumes/FRNU-1/data/eeg';

% Controls for plotting
freq_high = 12;
freq_low = 8;

a=0;
for i = 1:numel(evs)
    if evs(i).mfreq < freq_high && evs(i).mfreq > freq_low
        significance = evs(i).loc_percent > 40;
        if any((strcmp(evs(i).loc, 'anteriortemporal') & significance)) || any((strcmp(evs(i).loc, 'frontal') & significance))
            clusters(a+1) = evs(i);
            a=a+1;
        end
    end
end

subjs = {clusters.sub};

% Cross-subject
bds = braindata2.xLoadSubj(subjs, root);
bda = braindata2();
bda.loadAverage();

% Get data
angs_cued = [];
times_cued = [];
angs_n_cued = [];
times_n_cued = [];
xyzs = [];
chanNames = [];
for i=1:length(bds)
    % get tal and names for electrodes in cluster
    xyz = [array2table(clusters(i).elec_names', 'VariableNames', {'chanName'}),array2table(clusters(i).cltal, 'VariableNames', {'x','y','z'})];

    % get angles for all elecs in xyz at max wave prevalence 
    prevalence = nanmean(clusters(i).wave_logic(clusters(i).cued,:),1);
    [~,t_cued] = max(prevalence(1000:end));
    t_cued=t_cued+1000;

    prevalence = nanmean(clusters(i).wave_logic(clusters(i).n_cued,:),1);
    [~,t_n_cued] = max(prevalence(1000:end));
    t_n_cued=t_n_cued+1000;
    
    ang_cued = nan(height(xyz),1);
    ang_n_cued = nan(height(xyz),1);
    for elec = 1:height(xyz)
        % This might throw an error - look at the array being fed into convert_ang and the return - dims might be wrong
        ang_cued(elec) = circ_mean_nan(convert_angle_to_rad_in_yzplane(clusters(i).cltal, squeeze(clusters(i).ang(elec,t_cued,clusters(i).cued))')',1);
        ang_n_cued(elec) = circ_mean_nan(convert_angle_to_rad_in_yzplane(clusters(i).cltal, squeeze(clusters(i).ang(elec,t_cued,clusters(i).n_cued))')',1);
    end
    % populate data arrays
    angs_cued{i} = ang_cued;
    angs_n_cued{i} = ang_n_cued;
    times_cued{i} = repmat(t_cued,height(xyz.chanName),1);
    times_n_cued{i} = repmat(t_n_cued,height(xyz.chanName),1);
    xyzs{i} = xyz;
    chanNames{i} = xyz.chanName;
end

%% Cued Timing
% Average data
[t_cued_avg, m, data_mat] = braindata2.xLeadData2ROI(bds, times_cued, chanNames, 'minSubj',1);

% number of clusters that contributed to the data at each roic
cl_num = nan(size(data_mat,1),1);
for roic = 1:size(data_mat,1)
    cl_num(roic) = nnz(~isnan(data_mat(roic,:)));
end
inds = find(cl_num == 0);
cl_num(inds) = nan;

% Combine hemispheres
nROI = length(m) / 2;
ml   = m(1:nROI);
mr   = m(1+nROI:end);

% Get index of rh start
rhi  = sum(ml) + 1;

% Get ROI vertices on fsaverage
VPerROI_lh = bda.roic2roi(bda.roi.ROIC_lh(ml), 'lh');
VPerROI_rh = bda.roic2roi(bda.roi.ROIC_rh(mr), 'rh');
VPerROI    = [VPerROI_lh; VPerROI_rh];

figure()
% Plot
s1=subplot(121);
bp = bda.ezplot([],s1);
bp.plotRegionsData(VPerROI, t_cued_avg(m,:), 'rh_begin',rhi);
%bp.plotRegionsData(VPerROI, cl_num(m,:), 'rh_begin',rhi);
bp.view('right');
%colorbar();
title("Cued");

%% Not cued Timing
% Average data
[t_n_cued_avg, m, data_mat] = braindata2.xLeadData2ROI(bds, times_n_cued, chanNames, 'minSubj',1);

% number of clusters that contributed to the data at each roic
cl_num = nan(size(data_mat,1),1);
for roic = 1:size(data_mat,1)
    cl_num(roic) = nnz(~isnan(data_mat(roic,:)));
end
inds = find(cl_num == 0);
cl_num(inds) = nan;

% Combine hemispheres
nROI = length(m) / 2;
ml   = m(1:nROI);
mr   = m(1+nROI:end);

% Get index of rh start
rhi  = sum(ml) + 1;

% Get ROI vertices on fsaverage
VPerROI_lh = bda.roic2roi(bda.roi.ROIC_lh(ml), 'lh');
VPerROI_rh = bda.roic2roi(bda.roi.ROIC_rh(mr), 'rh');
VPerROI    = [VPerROI_lh; VPerROI_rh];

% Plot
s2 = subplot(122);
bp = bda.ezplot([],s2);
bp.plotRegionsData(VPerROI, t_n_cued_avg(m,:), 'rh_begin',rhi);
%bp.plotRegionsData(VPerROI, cl_num(m,:), 'rh_begin',rhi);
bp.view('right');

c = colorbar('Ticks',[1000,2000,3000,4000,5000],'TickLabels',{'1000 (cue on)','2000 (word on)','3000 (word off)','4000','5000'}, 'AxisLocation','out');
c.Label.String = 'Presentation Epoch (ms)';

%c = colorbar();
%c.Label.String = 'Number of Clusters';

title("Not Cued");


sgtitle("Timing of Max Theta Wave Prevalence Trends in Frontal Lobe and Temporal Lobe");
%sgtitle("Number of Electrodes Contributing to Measurement");
fontsize(gcf,30,'points');
set(gcf,'position',[0,0,2500,1000]);

%% Cued Direction 
% % Average data
% [~, m, data_mat] = braindata2.xLeadData2ROI(bds, angs_cued, chanNames, 'minSubj',1); % this is the wrong way to average circular coordinates
% angs_cued_avg = circ_mean_nan(data_mat,2);
% % want to flip across y axis so 0 is 180 (backwards)
% angs_cued_avg = phasewrap(angs_cued_avg); % wrap angles from -pi to pi (prob unecessary) 
% % Combine hemispheres
% nROI = length(m) / 2;
% ml   = m(1:nROI);
% mr   = m(1+nROI:end);
% 
% % Get index of rh start
% rhi  = sum(ml) + 1;
% 
% % Get ROI vertices on fsaverage
% VPerROI_lh = bda.roic2roi(bda.roi.ROIC_lh(ml), 'lh');
% VPerROI_rh = bda.roic2roi(bda.roi.ROIC_rh(mr), 'rh');
% VPerROI    = [VPerROI_lh; VPerROI_rh];
% 
% % Plot
% bp = bda.ezplot(); 
% cmap = phasemap();
% colormap(cmap);
% bp.plotRegionsData(VPerROI, angs_cued_avg(m,:), 'rh_begin',rhi);
% 
% phasebar('rad');
% title("Wave Direction from Frontal to Anterior Temporal Lobe");
% fontsize(gcf,30,'points');