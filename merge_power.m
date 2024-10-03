%% 3c. Zip together output of biowulf_swarm_pow.m into arrays for each event type

% grab sessions of interest 
clear;
close all;
load('/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct1.mat')

waveletFreqs = logspace(log10(2),log10(32),200);

for sess_i = 1:length(s)
    folder_name = fullfile('/Volumes/frnu-1/Rahil/Data', s(sess_i).pat_ID, s(sess_i).session);

    fprintf("Processing session %d \n", sess_i);
    patID = s(sess_i).pat_ID;

    try
        load(fullfile(folder_name, "sess_struct_pc.mat"));
    
        tal_struct= load(fullfile('/Volumes/FRNU/data/eeg',patID,'/tal/roi/lead_ROIC_LUT_monopolar.mat'));
    catch
        fprintf("File not found for %s, %s \n",patID, s(sess_i).session);
        continue;
    end
    %electrodes x trials x freq for each event type

    %% pres
    try
        data_path = fullfile(folder_name,'power','pres');
        load(data_path);

        % Populate corresponding position
        pres_pow = power;
        clear power;
    catch
        disp("No presentation events");
    end
    %% recog
    try
        data_path = fullfile(folder_name,'power','recog');
        load(data_path);

        % Populate corresponding position
        recog_pow = power;
        clear power;
    catch
        disp("No recognition events");
    end
    %% fr
    try
        data_path = fullfile(folder_name,'power','fr');
        load(data_path);

        % Populate corresponding position
        fr_pow = power;
        clear power;
    catch
        disp("No free recall events");
    end
    
    try
        power_struct.pres_pow=squeeze(nanmean(pres_pow,2));
    catch
        disp("Not populating presentation struct");
    end
    try
        power_struct.recog_pow=squeeze(nanmean(recog_pow,2)); 
    catch 
        disp("Not populating recognition struct");
    end
    try
        power_struct.fr_pow=squeeze(nanmean(fr_pow,2)); 
    catch
        disp("Not populating free recall struct");
    end
    try
        for el = 1:length(sess_struct.chan_names)
            tal{el,1} = tal_struct.lead_roi_lut(strcmp(tal_struct.lead_roi_lut.chanName, sess_struct.chan_names(el)),:).chanName;
            tal{el,2} = tal_struct.lead_roi_lut(strcmp(tal_struct.lead_roi_lut.chanName, sess_struct.chan_names(el)),:).x;
            tal{el,3} = tal_struct.lead_roi_lut(strcmp(tal_struct.lead_roi_lut.chanName, sess_struct.chan_names(el)),:).y;
            tal{el,4} = tal_struct.lead_roi_lut(strcmp(tal_struct.lead_roi_lut.chanName, sess_struct.chan_names(el)),:).z;
        end
        power_struct.tal=tal;
    catch
        fprintf('No electrodes in %s, %s \n', patID, s(sess_i).session);
        power_struct.tal=[];
    end

    power_struct.subject=patID;

    save(fullfile(folder_name, "pow.mat"), 'power_struct', '-v7.3');
    clear power_struct; clear pres_pow; clear recog_pow; clear fr_pow; clear tal; 
end
