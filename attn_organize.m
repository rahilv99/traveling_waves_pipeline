%% 2. attn_organize.m adapted from uma's pal_organize.m and audrey's pc1_organizing_pal.m

% organizes raw EEG data for power analysis, taking in sessions structures
% from folder 1) and using Julio’s load_eeg function;
% outputs sess_struct_pc.mat structures (dimensions channels x events x time)
% that contain data during PAL encoding/retrieval time-locked to study onset and response vocalization.
%% ERRORS
% starting data organizing for NIH090 session 8 index 401
% Found 48 evs_pres, 72 evs_recog, and 0 evs_fr
% NIH090 no bad chans
% Taking location parameter /Volumes/FRNU/data/eeg to be the rootEEGdir which contains the subject folder
% Warning: FileList is empty!
% > In GetFilesandChannels (line 440)
% In load_eegfile (line 295)
% In attn_organize (line 81)
% Finished loading eeg file
% Completed presentation events 3
% Completed recognition events 3
% Completed free recall events 3
% filled out sess_struct for NIH090_session_8\n
% rawr

%% setting-up and initializing

clear
close all
restoredefaultpath;
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/ZaghloulCodebase'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/macro_tws_attention'));
projfolder='/Volumes/Rahil_FRNU/Scripts';
eeg_path = '/Volumes/FRNU/data/eeg';

% loading sessions structures (created separately)
load("/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct2.mat")
s=s2;

%% extract + organize raw eeg data for each session: events, channels, time

% initialize time window information, in ms
pres_before = 2000; pres_after = 3000; % cue comes on for half a second and half second gap between * and word. word on for 1s
recog_before = 1000; recog_after = 5000; % goes off screen when patient responds, so liberally looking 5s after recog cue onset
fr_before = 1000; fr_after = 41000; % get 41 seconds after free recall onset
buffer = 1000; % buffer before and after time window

for sess = 212:length(s)
    patID = s(sess).pat_ID;
    session = s(sess).session;
    folder_name = fullfile('/Volumes/frnu-1/Rahil/Data', patID, session);
    fprintf("starting data organizing for " + patID + " session " + string(s(sess).sess_num) + "\n")
    try
        % load + clean events
        load(s(sess).events_file);
        if length(events)>200

            % get events for presentation, recognition memory cue, and freerecall. add to sessions structure that will be saved in folder for every session
            % Define presentation and retrieval events structures

            evs_pres = events(ismember({events.type},{'WORD_PRESENT'}) & [events.isSample] == 1);
            evs_recog = events(ismember({events.type},{'FORCED_CHOICE_0ALT'}) & [events.isTest]==1);
            evs_recog = sortStruct(evs_recog,{'blockCount','testListIndex'},{'ascend','ascend'});

            % make copy with RTs
            evs_fr = events(ismember({events.type},{'FREE_RECALL_START'}));

            % Get field values as cell array
            fieldValues = {evs_fr.FRannotateCounts};

            % Find indices of empty elements
            emptyIndexes = find(cellfun(@(x) isempty(x) || ~isnumeric(x) && numel(x) == 0 || isnumeric(x) && numel(x) > 0 && x(1) == 0, fieldValues));

            logicalIndex = logical(ones(1, numel(evs_fr)));
            logicalIndex(emptyIndexes) = 0;

            % Create new struct with excluded element
            evs_fr = evs_fr(logicalIndex);

            fprintf('Found %d evs_pres, %d evs_recog, and %d evs_fr \n',size(evs_pres,2),size(evs_recog,2),size(evs_fr,2))
            sess_struct.evs_pres = evs_pres;
            sess_struct.evs_recog = evs_recog;
            sess_struct.evs_fr = evs_fr;

            % grab timestamp of raw eeg
            timestamp = split(evs_pres(1).eegfile,'/'); timestamp = timestamp{end};

            % subject level bad chans
            try
                [allbadchans] = sub_bad_chans(s(sess).pat_ID, s, s(sess).task);

                % load eeg data (global average and monopolar and subdural only)
                [eeg_ts, Fs, channelMap, chan_names, sess_info, eegdir] = load_eegfile(s(sess).pat_ID, timestamp, 0, Inf...
                    , 'loc_detrend', 1, 'rmline', 1, 'rem_sat', 10, 'ref', 'global_avg', 'location', eeg_path, 'bad_chans', allbadchans,'filter','gen_Lowpass_FIR');
            catch
                disp([s(sess).pat_ID, ' no bad chans'])

                % load eeg data (global average and monopolar and subdural only)
                [eeg_ts, Fs, channelMap, chan_names, sess_info, eegdir] = load_eegfile(s(sess).pat_ID, timestamp, 0, Inf...
                    , 'loc_detrend', 1, 'rmline', 1, 'rem_sat', 10, 'ref', 'global_avg', 'location', eeg_path,'filter','gen_Lowpass_FIR');

            end

            disp('Finished loading eeg file')
            % initialize for raw timeseries
            raw.pres = NaN(length(chan_names), size(evs_pres,2), pres_before + pres_after + (2*buffer)); % presentation raw timeseries, initialize channels x events x time
            raw.recog = NaN(length(chan_names), size(evs_recog,2), recog_before + recog_after + (2*buffer)); % recognition raw timeseries, initialize channels x events x time
            raw.fr = NaN(length(chan_names), size(evs_fr,2), fr_before + fr_after + (2*buffer)); % fr raw timeseries, initialize channels x events x time

            % assign eegoffsets for presentation + response
            pres_eegoffsets = [evs_pres.eegoffset]; % eegoffsets for presentation, 1 x # of events
            recog_eegoffsets = [evs_recog.eegoffset]; % eegoffsets for response, 1 x # of events
            fr_eegoffsets = [evs_fr.eegoffset]; % eegoffsets for cue, 1 x # of events

            % loop thru all channels and events and assign sliced time series to matrix
            for chan_i = 1:length(chan_names)
                for event_i = 1:size(evs_pres, 2)
                    % PRESENTATION (study pair stimulus)
                    p_start = pres_eegoffsets(event_i) - pres_before - buffer; % start time for event window + buffer
                    p_end = pres_eegoffsets(event_i) + pres_after + buffer - 1; % end time for event window + buffer

                    % fill matrices
                    raw.pres(chan_i, event_i, :) = eeg_ts(p_start:p_end, chan_i); % for every channel, fill out sliced time series for that event

                end
            end
            fprintf('Completed presentation events %d \n', event_i)
            for chan_i = 1:length(chan_names)
                for event_i = 1:size(evs_recog, 2)

                    % RESPONSE (vocalization)
                    r_start = recog_eegoffsets(event_i) - recog_before - buffer; % start time for event window + buffer
                    r_end = recog_eegoffsets(event_i) + recog_after + buffer - 1; % end time for event window + buffer
                    raw.recog(chan_i, event_i, :) = eeg_ts(r_start:r_end, chan_i);

                end
            end
            fprintf('Completed recognition events %d \n', event_i)
            for chan_i = 1:length(chan_names)
                for event_i = 1:size(evs_fr, 2)

                    f_start = fr_eegoffsets(event_i) - fr_before - buffer; % start time for event window + buffer
                    f_end = fr_eegoffsets(event_i) + fr_after + buffer - 1; % end time for event window + buffer

                    raw.fr(chan_i, event_i, :) = eeg_ts(f_start:f_end, chan_i);
                end
            end
            fprintf('Completed free recall events %d \n', event_i)


            % assign raw eeg matrices to sess_struct
            sess_struct.raw_pres = raw.pres;
            sess_struct.raw_recog = raw.recog;
            sess_struct.raw_fr = raw.fr;
            sess_struct.chan_names = chan_names;
            sess_struct.timestamp = timestamp;

            disp("filled out sess_struct for " + s(sess).pat_ID + "_" + s(sess).session + "\n");


            % save to folder: sess_struct + mean_spectral.png
            if ~exist(folder_name, 'dir') % make folder if doesn't exist already
                mkdir(folder_name);
                disp("folder created for " + s(sess).pat_ID + "_" + s(sess).session);
            end
            save(fullfile(folder_name, "sess_struct_pc"), 'sess_struct', '-v7.3'); % save sess_struct to folder
            disp("folder already exists for " + s(sess).pat_ID + "_" + s(sess).session);
        end
    catch e
        fprintf("Error on patient %s, %s\n", patID, session);
        fprintf(1,'The identifier was: %s\n',e.identifier);
        fprintf(1,'There was an error! The message was: %s\n',e.message);
    end
end
%
%