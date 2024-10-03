%% add paths
addpath(genpath('/Volumes/FRNU-1/people/Rahil/ZaghloulCodebase'));
addpath(genpath('/Volumes/FRNU-1/people/Rahil/Scripts/macro_tws_pal'));
addpath(genpath('/Volumes/FRNU-1/people/Rahil/Scripts/circstat-matlab-master'));
addpath(genpath('/Volumes/FRNU-1/people/Rahil/Scripts/(7) Analysis/functions'));

clear; close all;

patID = 'NIH032';
sessID = 'session_1';

load('/Volumes/Rahil_FRNU/Scripts/all_pres.mat');
evs = all_pres;

folder_name = fullfile('/Volumes/Rahil_FRNU/Data', patID,sessID);

try
    ind = find(strcmp(patID,{evs.sub}) & strcmp(sessID,{evs.sess}));  % may need to change this line

    clust = evs(ind);
    root = '/Volumes/FRNU-1/data/eeg';

    freq = clust.mfreq;

    % Find direction (angle) of wave 1000 to 2000 ms
    wave_angle = clust.ang; %electrode x time x event

    
    % Seperately average cued and uncued events -- CHANGE FOR SPECIFIC TIME AND EVENT
    cued_angle = squeeze(mean(wave_angle(:,2001:3000,clust.cued),3, 'omitnan'));
    uncued_angle = squeeze(mean(wave_angle(:,2001:3000,clust.n_cued),3, 'omitnan'));
    cued_angle = squeeze(mean(cued_angle(:,clust.wave_logic(2001:3000)),2, 'omitnan'));
    uncued_angle = squeeze(mean(uncued_angle(:,clust.wave_logic(2001:3000)),2, 'omitnan'));

    cued_ang = convert_angle_to_rad_in_yzplane(clust.cltal, cued_angle)';
    uncued_ang = convert_angle_to_rad_in_yzplane(clust.cltal, uncued_angle)';

    figure()
    %% Brain Plotter
    handle = subplot(1,2,1);
    bd = braindata2(patID, root);
    bp_cued = bd.ezplotEcog([],handle);

    bp_cued.setOpacity(0.1);

    bp_cued.clearPoints();
    t = bd.tal.xyz;
    bp_cued.plotPoint(t, 'color', [0.2 0.2 0.2], 'radius', 0.5);

    if (sum(bd.tal.xyz.x > 0)/length(bd.tal.xyz.x)) > .5 % if over half of electrodes are on right side
        bp_cued.view('right');
    else
        bp_cued.view('left');
    end

    clust_tal = clust.cltal;
    freq = clust.mfreq;
    hold on;
    for el = 1:length(clust_tal)

        x = clust_tal(el,1);                          % X coordinate of arrow start
        y = clust_tal(el,2);                          % Y coordinate of arrow start
        z = clust_tal(el,3);                          % Z coordinate of arrow start
        theta = cued_ang(el);                   % Angle of arrow, from x-axis
        L = 1;                          % Length of arrow
        yEnd = y+L*cos(theta);          % X coordinate of arrow end
        zEnd = z+L*sin(theta);          % Y coordinate of arrow end
        arrow3([x y z], [x yEnd zEnd], 'b0',1,2.5);
    end
    title("Cued (Averaged Across Events and Encoding Epoch)")
    % uncued plotting
    handle = subplot(1,2,2);
    bp_uncued = bd.ezplotEcog([],handle);

    bp_uncued.setOpacity(0.1);

    bp_uncued.clearPoints();
    t = bd.tal.xyz;
    bp_uncued.plotPoint(t, 'color', [0.2 0.2 0.2], 'radius', 0.5);

    if (sum(bd.tal.xyz.x > 0)/length(bd.tal.xyz.x)) > .5 % if over half of electrodes are on right side
        bp_uncued.view('right');
    else
        bp_uncued.view('left');
    end

    hold on;
    for el = 1:length(clust_tal)

        x = clust_tal(el,1);                          % X coordinate of arrow start
        y = clust_tal(el,2);                          % Y coordinate of arrow start
        z = clust_tal(el,3);                          % Z coordinate of arrow start
        theta = uncued_ang(el);                   % Angle of arrow, from x-axis
        L = 1;                          % Length of arrow
        yEnd = y+L*cos(theta);          % X coordinate of arrow end
        zEnd = z+L*sin(theta);          % Y coordinate of arrow end
        arrow3([x y z], [x yEnd zEnd], 'b0',1,2.5);
    end
    title("Not Cued (Averaged Across Events and Encoding Epoch)")

    fontsize(gcf,28,'points')
    out_dir = fullfile('/Volumes/Rahil_FRNU/Figures/Patient/Brain_Plots');

    if ~exist(out_dir, 'dir') % make folder if doesn't exist already
        mkdir(out_dir);
    end 
    
    saveas(gcf, fullfile(out_dir, append(patID, '_', sessID, '_bp_angle.png')));
catch e
    fprintf("Error on patient %s, %s\n", patID, sessID);
    fprintf(1,'The identifier was: %s\n',e.identifier);
    fprintf(1,'The error was: %s\n',e.message);
    close;
end