%% 3a. Biowolf swarm script generator for pal_power (stage 3 of data pipeline)
% Originally from Noah
% Date: 6/14/2024 - Adapted by Rahil Verma (rahil.verma@duke.edu)
% Creates the text file for biowulf to read and execute
clear;
close all;

% set biowulf directory (make sure you have CMND + K’ed already) to:
% smb://hpcdrive.nih.gov/data
% THIS IS WHERE THE TEXT FILE WILL BE WRITTEN TO -> CHECK THAT THIS IS YOUR
% COMPUTER EXTENSION FOR FRNU/people/Rahil/Scripts/(3) power
%biowulfSaveDir = 'Y:\Rahil\Scripts\Code\Biowolf';
biowulfSaveDir = '/Volumes/frnu-1/Rahil/Power';
cd(biowulfSaveDir)

%% Fill out important information here:
txtFile.name = 'vermar5_attn_pipeline_power.swarm.txt';
txtFile.jobName = 'vermar5_attn_pipeline_power';
txtFile.numTasks = 2; % should always be 2 * number of parallel workers (plus a few more to not cause delays or crashes)
txtFile.GB = 6; % some of the 'norm' partitions have a maximum memory of 247 GB. Keep this number under 247 for smaller wait times. % usually 6 - if reruns, increase to 9
% if you need to go above 350 GB, use the 'largemem' partition
txtFile.scratch = '--gres=lscratch:1'; % it's a good idea to ask for the same amount of /lscratch/ space as you need for the function itself
txtFile.function = 'attn_pow_biowolf_v1'; % the name of the function
txtFile.time = '01:30:00'; % I'm pretty sure the format is: DAYS:HOURS:MINUTES:SECONDS

txtFile.introText = 'matlab -nodisplay -r';
txtFile.outroText = 'exit;" ';
txtFile.dir = '/data/FRNU/Rahil/Power'; % Biowulf location where you store your script

fileID = fopen(txtFile.name, 'w');

% write first few required lines
fprintf(fileID, '#SWARM --job-name=%s\n', txtFile.jobName);
fprintf(fileID, '#SWARM --module=matlab\n');
fprintf(fileID, '#SWARM -t %d -g %d --time %s %s\n', ...
    txtFile.numTasks, txtFile.GB, txtFile.time, txtFile.scratch);

fprintf(fileID, '\n \n'); % line break

% grab sessions of interest 
load("/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct2.mat")
s = s2;
%% Here is where you loop over conditions to create the appropriate inputs to your function.

for sess_i = 1:length(s)

    folder_name = append('/data/FRNU/Rahil/Data/', s(sess_i).pat_ID, '/', s(sess_i).session);

    fprintf('Session %d \n', sess_i);

    ev_type = 'pres';
    functText = sprintf("%s('%s', '%s'); exit;", ...
        txtFile.function, folder_name, ev_type);
    fprintf(fileID, '%s %scd(%s%s%s); %s%s\n', ...
        txtFile.introText, '"', "'", txtFile.dir,  "'", functText, '"');

    ev_type = 'recog';
    functText = sprintf("%s('%s', '%s'); exit;", ...
        txtFile.function, folder_name, ev_type);
    fprintf(fileID, '%s %scd(%s%s%s); %s%s\n', ...
        txtFile.introText, '"', "'", txtFile.dir,  "'", functText, '"');

    ev_type = 'fr';
    functText = sprintf("%s('%s', '%s'); exit;", ...
        txtFile.function, folder_name, ev_type);
    fprintf(fileID, '%s %scd(%s%s%s); %s%s\n', ...
        txtFile.introText, '"', "'", txtFile.dir,  "'", functText, '"');
end
fprintf(fileID, '\n');

% save out and close the txtfile
fclose(fileID);
disp('done :) Keep on doing great things in life')

%% Information to include in your biowulf function
% Instead of starting a parallel processing operation with parfor, you NEED
% to change where matlab updates its information about the parallel job.
% Reading and writing files to /lscratch/ has been found (by Noah, the wise and
% beautiful) to be MUCH faster than keeping it on your home directory or
% the data directory on biowulf. Use the following few lines (uncomment)
% Similar to how you should throw your dataWhen you need to start a parallel processing script with
% pc = parcluster('local');
% pc.JobStorageLocation = strcat('/lscratch/', getenv('SLURM_JOB_ID'));
% n_workers = 16; parpool(pc, n_workers, 'IdleTimeout', 240)
% n_workers = 16; parpool(pc, n_workers, 'IdleTimeout', 240)