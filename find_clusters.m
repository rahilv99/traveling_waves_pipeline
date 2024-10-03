%% 4. find_clusters_pal.m adapted from honghui's code to find cluseters in sternberg wm task
% load power structure for each session (channel x frequency) separately
% for encoding, cue, and retrieval. use tal to cluster adjacent electrodes
% with peak frequencies within 2Hz

%% setting up + initializing

% set path
clear; 
close all;

addpath(genpath('/Volumes/Rahil_FRNU/Scripts/ZaghloulCodebase'));
%addpath(genpath("Z:\people\Rahil\Scripts\ZaghloulCodebase"));

addpath(genpath('/Volumes/Rahil_FRNU/Scripts'));
%addpath(genpath("Z:\people\Rahil\Scripts"));

projfolder='/Volumes/Rahil_FRNU/Scripts/(4) clusters';
%projfolder="Z:\people\Rahil\Scripts";

% grab sessions of interest 
load("/Volumes/Rahil_FRNU/Scripts/attn_sessions_struct1.mat");

%load("Z:\people\Rahil\Scripts\allFileNames.mat");

% set up timing
% time_conds = {[-2500:5499], [-9000:2499], [-2500:5499]}; % for presentation, response, and cue

%% loop through sessions + load  power structs


data_path = '/Volumes/frnu-1/Rahil/Data';
%data_path = 'Y:\Rahil\Data';
waveletFreqs = logspace(log10(2),log10(32),200); % 200 wavelets from 2-32 Hz, log space
resamp=1000;
w  = 6;


for sess_i = 1:length(s)
    patID = s(sess_i).pat_ID;
    sess = s(sess_i).session;
    folder_name = fullfile(data_path, patID, sess);

    fprintf('Processing session %d \n', sess_i);

    try
        load(fullfile(folder_name, "pow.mat"));
        mono_elecs = readtable(fullfile('/Volumes/FRNU/data/eeg',patID,'/tal/atlas/atlas_monopolar_simple.csv'));
    catch
        fprintf('Files not found for %s, %s \n', patID, sess);
        continue;
    end

    ClusterN=0;

    %     clear ConnectM electrode_peak_frequency

    % load power and electrode infomation
    if ~isempty(power_struct.tal)
        %     result=AllGrid(s);
        %     pow=result.Pow_bad'; % get the power
        %     tal=result.tal; % get the electrode infos
        clear electrode_peak_frequency;
        n=length(power_struct.tal);  % count number of electrodes
        electrode_peak_frequency = NaN(n,3);
        Tal=[cell2mat(power_struct.tal(:,2)) cell2mat(power_struct.tal(:,3)) cell2mat(power_struct.tal(:,4))];

        % Populate power information by event type
        for i =1:n
            try
                pres_power = squeeze(power_struct.pres_pow(i,:));
                electrode_peak_frequency(i,1)=getpeak2(pres_power,waveletFreqs,2,32,1);% find peak frequencies of each electrode
            catch
                disp("No presentation events");
                break;
            end
        end

        for i =1:n
            try
                recog_power = squeeze(power_struct.recog_pow(i,:));
                electrode_peak_frequency(i,2)=getpeak2(recog_power,waveletFreqs,2,32,1);% find peak frequencies of each electrode
            catch
                disp("No recognition events");
                break;
            end
        end

        for i =1:n
            try
                fr_power = squeeze(power_struct.recog_pow(i,:));
                electrode_peak_frequency(i,3)=getpeak2(fr_power,waveletFreqs,2,32,1);% find peak frequencies of each electrode
            catch
                disp("No free recall events");
                break;
            end
        end

        steplength=1;
        windowLength=2;
        v=zeros(1,length(1:steplength:30));
        nPeaks=1;
        window=(2:steplength:30);
        distance=25;%distance threshold between electrodes to get neighboring electrodes
        %% Presentation
        if exist('pres_power', 'var')
            for k=window % loop through each window
                for i=1:n % loop through each electrode
                    ff=[electrode_peak_frequency(i,1)];
                    ff(ff<k)=[];
                    ff(ff>k+windowLength)=[];
                    if ~isempty(ff)     % Find if there's peak in the window for that electrode.
                        v(nPeaks)=v(nPeaks)+1;  % count the number of electrode which has a peak frequency in the window
                    end
                end
                nPeaks=nPeaks+1; % move to next window

            end
            %     x=localmax(v);  % find which frequencies, there are more electrode oscillating together.
            [~,locs]=findpeaks(v);  % find which frequencies, there are more electrode oscillating together.
            for k=window(locs) % loop through the frequencies there are local maximum electrode oscillating together
                ConnectM=zeros(n,n);  % create a Connectom matrix
                for i=1:n-1
                    for j=i+1:n
                        fi=[electrode_peak_frequency(i,1)]; % peaks of electrode i
                        fi(fi<k)=[];
                        fi(fi>k+windowLength)=[];  % check for eletrode i, if there's a peak in this window

                        fj=[electrode_peak_frequency(j,1)]; % peaks of electrode j
                        fj(fj<k)=[];
                        fj(fj>k+windowLength)=[];  % check for eletrode j, if there's a peak in this window
                        if (~isempty(fi))&&(~isempty(fj))&&norm(Tal(i,:)-Tal(j,:))<distance    % if both electrode have peak in this window and they are close.
                            ConnectM(i,j)=1; % set the values to 1 to matrix(i,j)
                            ConnectM(j,i)=1;
                        end
                    end
                end
                DG=sparse(ConnectM);  % calculate a sparse matrix using the Connectom matrix
                [S,C] = graphconncomp(DG);   % get groups of spatial clustered electrodes
                for i= 1:S     % for each group
                    if length(C(C==i))>4  % if there are 4 or more electrode in this group, we count they are oscillation cluster.
                        disp("Found presentation cluster!");
                        ClusterN=ClusterN+1;
                        sess_clust(ClusterN).patient=patID; % find the subject
                        sess_clust(ClusterN).pres_tal=Tal(C==i,:); % find the electrode info
                        sess_clust(ClusterN).pres_channame=[power_struct.tal{C==i,1}]; % find the electrode number
                        sess_clust(ClusterN).pres_elecnum = find(ismember(mono_elecs.chanName, [power_struct.tal{C==i,1}]) == 1);
                        ff=[electrode_peak_frequency(C==i,1)];   % get the frequency
                        ff(ff<k)=[];
                        ff(ff>k+windowLength)=[];

                        sess_clust(ClusterN).pres_freq=mean(ff);    % the center frequency of the group is the mean of each electrode.
                        dist_elec_mat = zeros(n,n);
                        for el1=1:n-1
                            for el2=el1+1:n
                                % if both electrode have peak in this window and they are close.
                                dist_elec_mat(el1,el2)=norm(Tal(el1,:)-Tal(el2,:)); % set the values to 1 to matrix(i,j)
                                dist_elec_mat(el2,el1)=norm(Tal(el1,:)-Tal(el2,:));
                            end
                        end
                        sess_clust(ClusterN).pres_dist_elec = dist_elec_mat(C == i,C == i);
                        % count next cluster
                    end
                end
            end
            disp("Presentation events completed");
        end
        %% recog clusters
        if exist('recog_power', 'var')
            nPeaks=1;

            for k=window % loop through each window
                for i=1:n % loop through each electrode
                    ff=[electrode_peak_frequency(i,2)];
                    ff(ff<k)=[];
                    ff(ff>k+windowLength)=[];
                    if ~isempty(ff)     % Find if there's peak in the window for that electrode.
                        v(nPeaks)=v(nPeaks)+1;  % count the number of electrode which has a peak frequency in the window
                    end
                end
                nPeaks=nPeaks+1; % move to next window

            end
            %     x=localmax(v);  % find which frequencies, there are more electrode oscillating together.
            [~,locs]=findpeaks(v);  % find which frequencies, there are more electrode oscillating together.
            for k=window(locs) % loop through the frequencies there are local maximum electrode oscillating together
                ConnectM=zeros(n,n);  % create a Connectom matrix
                for i=1:n-1
                    for j=i+1:n
                        fi=[electrode_peak_frequency(i,2)]; % peaks of electrode i
                        fi(fi<k)=[];
                        fi(fi>k+windowLength)=[];  % check for eletrode i, if there's a peak in this window

                        fj=[electrode_peak_frequency(j,2)]; % peaks of electrode j
                        fj(fj<k)=[];
                        fj(fj>k+windowLength)=[];  % check for eletrode j, if there's a peak in this window
                        if (~isempty(fi))&&(~isempty(fj))&&norm(Tal(i,:)-Tal(j,:))<distance    % if both electrode have peak in this window and they are close.
                            ConnectM(i,j)=1; % set the values to 1 to matrix(i,j)
                            ConnectM(j,i)=1;
                        end
                    end
                end
                DG=sparse(ConnectM);  % calculate a sparse matrix using the Connectom matrix
                [S,C] = graphconncomp(DG);   % get groups of spatial clustered electrodes
                for i= 1:S     % for each group
                    if length(C(C==i))>4  % if there are 4 or more electrode in this group, we count they are oscillation cluster.
                        disp("Found recognition cluster!");
                        ClusterN=ClusterN+1;
                        sess_clust(ClusterN).patient=patID ; % find the subject
                        sess_clust(ClusterN).recog_tal=Tal(C==i,:); % find the electrode info
                        sess_clust(ClusterN).recog_channame=[power_struct.tal{C==i,1}]; % find the electrode number
                        sess_clust(ClusterN).recog_elecnum = find(ismember(mono_elecs.chanName, [power_struct.tal{C==i,1}]) == 1);
                        ff=[electrode_peak_frequency(C==i,2)];   % get the frequency
                        ff(ff<k)=[];
                        ff(ff>k+windowLength)=[];

                        sess_clust(ClusterN).recog_freq=mean(ff);    % the center frequency of the group is the mean of each electrode.
                        dist_elec_mat = zeros(n,n);
                        for el1=1:n-1
                            for el2=el1+1:n
                                % if both electrode have peak in this window and they are close.
                                dist_elec_mat(el1,el2)=norm(Tal(el1,:)-Tal(el2,:)); % set the values to 1 to matrix(i,j)
                                dist_elec_mat(el2,el1)=norm(Tal(el1,:)-Tal(el2,:));
                            end
                        end
                        sess_clust(ClusterN).recog_dist_elec = dist_elec_mat(C == i,C == i);
                        % count next cluster
                    end
                end
            end
            disp("Recognition events completed");
        end
        %% free recall
        if exist('fr_power', 'var')
            nPeaks=1;

            for k=window % loop through each window
                for i=1:n % loop through each electrode
                    ff=[electrode_peak_frequency(i,1)];
                    ff(ff<k)=[];
                    ff(ff>k+windowLength)=[];
                    if ~isempty(ff)     % Find if there's peak in the window for that electrode.
                        v(nPeaks)=v(nPeaks)+1;  % count the number of electrode which has a peak frequency in the window
                    end
                end
                nPeaks=nPeaks+1; % move to next window

            end
            %     x=localmax(v);  % find which frequencies, there are more electrode oscillating together.
            [~,locs]=findpeaks(v);  % find which frequencies, there are more electrode oscillating together.
            for k=window(locs) % loop through the frequencies there are local maximum electrode oscillating together
                ConnectM=zeros(n,n);  % create a Connectom matrix
                for i=1:n-1
                    for j=i+1:n
                        fi=[electrode_peak_frequency(i,3)]; % peaks of electrode i
                        fi(fi<k)=[];
                        fi(fi>k+windowLength)=[];  % check for eletrode i, if there's a peak in this window

                        fj=[electrode_peak_frequency(j,3)]; % peaks of electrode j
                        fj(fj<k)=[];
                        fj(fj>k+windowLength)=[];  % check for eletrode j, if there's a peak in this window
                        if (~isempty(fi))&&(~isempty(fj))&&norm(Tal(i,:)-Tal(j,:))<distance    % if both electrode have peak in this window and they are close.
                            ConnectM(i,j)=1; % set the values to 1 to matrix(i,j)
                            ConnectM(j,i)=1;
                        end
                    end
                end
                DG=sparse(ConnectM);  % calculate a sparse matrix using the Connectom matrix
                [S,C] = graphconncomp(DG);   % get groups of spatial clustered electrodes
                for i= 1:S     % for each group
                    if length(C(C==i))>4  % if there are 4 or more electrode in this group, we count they are oscillation cluster.
                        disp("Found free recall cluster!");
                        ClusterN=ClusterN+1;
                        sess_clust(ClusterN).patient=patID ; % find the subject
                        sess_clust(ClusterN).fr_tal=Tal(C==i,:); % find the electrode info
                        sess_clust(ClusterN).fr_channame=[power_struct.tal{C==i,1}]; % find the electrode number
                        cl_el_num  = find(ismember(mono_elecs.chanName, [power_struct.tal{C==i,1}]) == 1);
                        sess_clust(ClusterN).fr_elecnum =cl_el_num;
                        ff=[electrode_peak_frequency(C==i,3)];   % get the frequency
                        ff(ff<k)=[];
                        ff(ff>k+windowLength)=[];

                        sess_clust(ClusterN).fr_freq=mean(ff);
                        dist_elec_mat = zeros(n,n);
                        for el1=1:n-1
                            for el2=el1+1:n
                                % if both electrode have peak in this window and they are close.
                                dist_elec_mat(el1,el2)=norm(Tal(el1,:)-Tal(el2,:)); % set the values to 1 to matrix(i,j)
                                dist_elec_mat(el2,el1)=norm(Tal(el1,:)-Tal(el2,:));
                            end
                        end
                        sess_clust(ClusterN).fr_dist_elec = dist_elec_mat(C == i,C == i);

                        % the center frequency of the group is the mean of each electrode.
                        % count next cluster
                    end
                end
            end
            disp("Free Recall events completed");
        end
    end
    if exist("sess_clust",'var')
        file_name = fullfile(folder_name, "clusters.mat");
        save(file_name, 'sess_clust', '-v7.3');
        disp('Cluster saved');
    end 

    clear sess_clust; clear pres_power; clear recog_power; clear fr_power;
end


