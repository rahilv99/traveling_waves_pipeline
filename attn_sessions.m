%% 1. attn_sessions.m: organizing behavioral sessions structure adapted from pal_sessions.m and Audrey's dc0_sessions.m
 
% - input: grabs from eeg_path and makes structure for behavioral PAL events
% - makes the sessions structure that will be used to cut up trials for PAL
% - uses key function getBehEv (KS's function to get patient IDs and folders of interest from EEG)
% - can customize what sessions you want here (some of these sorted by ZX and SB)
% - output: a structure (pal_sessions_structure.mat) with behavioral information for the sessions you want to look at (this structure will be used in future scripts)
clear;
close all;
%% set-up
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/ZaghloulCodebase'));
addpath(genpath('/Volumes/Rahil_FRNU/Scripts/macro_tws_attention'));

% local eeg root dir
eeg_path = '/Volumes/FRNU/data/eeg'; 

% make structure for sessions of PAL (but NOT pa3)
[s] = getBehEv({'attentionTask';'attentionFeedback'},eeg_path); % this function gets patient IDs and folders of interest from EEG

% sort structure in ascending order
s = sortStruct(s,{'pat_num'},{'ascend'});
% grab sessions of interest 

% initialize session + patient number info
allsubs = unique([s.pat_num]); % grab unique patient numbers

%% loop through all subjects and get sessions of interest info
for pat_num = allsubs
    
    % all original sessions
    log_sub_index = ([s.pat_num] == pat_num); % getting index for current subject sessions
    sess_num_array = [s.sess_num];
    orig_sessions = sess_num_array(log_sub_index); % grabbing all sessions for current subject
    
    % manual sorting by ZX & SB
%     if pat_num==38
%         sessions=[1]; % 38 session 0 is practice
%     elseif pat_num==37
%         sessions=[0 1 5];
%     elseif pat_num==42
%         sessions=[0 1 2];
%     elseif pat_num==43
%         sessions=[2 3 4 5]; % 43 has no session 0, and session 1 is practice
%     elseif pat_num==36
%         sessions=[0 1 2];
%     elseif pat_num==34
%         sessions=[0];
%     elseif pat_num==76
%         sessions = [0 1];
%     elseif pat_num==29
%         sessions=[0 1 2 3];
%     elseif pat_num==39
%         sessions=[0 1 2 4];
%     elseif pat_num==42
%         sessions=[0 1 2];
%     elseif pat_num==46
%         sessions=[0 1 2];
%     elseif pat_num==47
%         sessions=[0 1 2];
%     elseif pat_num==50 % 0 1 2 are repeated sessions
%         sessions=[0 3];
%     % exclude these patients (30 + 45 spanish; 49 too few correct trials; 65/71/83 all depths)
%     elseif pat_num==30 || pat_num==45 || pat_num==49 || pat_num==65 || pat_num==71 || pat_num==83 
%         s = s(~([s.pat_num] == pat_num));
%         continue
%     else
        sessions = orig_sessions; % for all other sessions, grab all corresponding sessions
%     end
    
    % only grab first two sessions!  
%     if length(sessions) >= 2 % grab only first 2 sessions if more than 2 
%         sessions = sessions(1:2);
%     end
    
    % exclude sessions not wanted from s
    sessions_exclude = setdiff(orig_sessions, sessions)
    s = s(~([s.pat_num] == pat_num & ismember([s.sess_num], sessions_exclude)))

end


 save('/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct_full.mat','s');


%% get subjects and sessions of choice structure
% basically saves a struct to make it easier to load in subject/session info for future scripts

% subjects_of_choice = [26,28,29,32,34,35,36,37,39,40,41,42,46,47,50,52,55,75,76,62];
subjects_of_choice = unique([s.pat_num]);
sessions_of_choice = [];

soc_i = 1;
for i = 1:length(s)
    if ismember(s(i).pat_num, subjects_of_choice)
        sessions_of_choice(soc_i) = i;
        soc_i = soc_i + 1;
    end
end

soc.subjects_of_choice = subjects_of_choice;
soc.sessions_of_choice = sessions_of_choice;


%save('/Volumes/Rahil_FRNU/Scripts/soc_1.mat', '-struct', 'soc');
 
