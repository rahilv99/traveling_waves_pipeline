%% 5c. Merge processed data from biowolf into wave_info struct for for analysis in future scripts
clear;
% grab sessions of interest 
load("/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct2.mat");
s=s2;

for sess_i = 1:length(s)
    %get patID
    patID = s(sess_i).pat_ID;
    sessID = s(sess_i).session;

    folder_name = fullfile('/Volumes/frnu-1/Rahil/Data', patID, sessID);

    check = sprintf("Processing patient %s, session %s", patID, sessID);
    disp(check) 

    try
        load(fullfile(folder_name, 'Processing', "allwave.mat"));
    catch
        out = sprintf("allwave.mat file not found for patient %s session %s", patID, sessID);
        disp(out)
        continue;
    end
    % Wave_info files were not populated correctly by preprocessing
    % Some other patients exist in the wave_info files, removing here
    for i = 1:length(wave_info)
        % Check if sub field equals current patient
        if ~strcmp(wave_info(i).sub, patID)
            % If no, delete all following rows and continue
            while i <= length(wave_info)
                wave_info(i) = [];
            end
            break;
        end
    end
    for ct = 1: length(wave_info)
        
        dims = wave_info(ct).dims; % chan x event x time
        % Initialize variables to be returned
        ang=zeros(dims(1),dims(3),dims(2));
        spatial_freq=zeros(dims(1),dims(3),dims(2));
        cl_corr=zeros(dims(1),dims(3),dims(2));

        for el = 1:wave_info(ct).elecs

            data_path = fullfile(folder_name,'Processing', num2str(ct), num2str(el));
            % Attempt to load file, if DNE continue to next electrode
            try
                load(data_path);
            catch
                out = sprintf('Data for cluster %d electrode %d not found. That is okay. It was removed because it does not meet requirements.', ct, el);
                disp(out)
                continue;
            end
            try
                ang(el,:,:) = el_data.ang_el; 
                spatial_freq(el,:,:) = el_data.spatial_freq; 
                cl_corr(el,:,:) = el_data.cl_corr; 
            catch
                out = sprintf("Array size mismatch for electrode %d in cluster %d", el, ct);
                disp(out);
            end
        end
        try
            wave_info(ct).CLcorrelation=cl_corr;
            wave_info(ct).spatial_frequency=spatial_freq;
            wave_info(ct).ang=ang;
        catch
            out = sprintf("Could not populate wave_info file for patient %d, %d", patID, sessID);
            disp(out);
        end

        clear ang; clear spatial_freq; clear cl_corr;
    end

    save(fullfile(folder_name, "allwave.mat"), 'wave_info', '-v7.3');

    clear wave_info; 
end