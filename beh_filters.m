%% 6. Retrieve patients that perform reasonably well on task
clear;
% eeg_path = 'C:\Users\rahil\Box\AttentionWaves\Data';
eeg_path = '/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Data';

load("/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Scripts/attn_sessions_struct_full.mat");
% load("C:\Users\rahil\Box\AttentionWaves\Scripts\attn_sessions_struct_full.mat");

%if only runnining missing sub, load relevant struct
% ldpath='/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Scripts/';
% load([ldpath 'soc.mat']);
% load([ldpath 'cued_accuracy.mat']);%for c_accuracy
% load([ldpath 'n_cued_accuracy.mat']);%for n_accuracy
% load([ldpath 'foil_accuracy.mat']);%for f_accuracy

soc = [];
c_accuracy = [];
n_accuracy = [];
f_accuracy = [];

soc_i = 1;
cnt=0;%set missing dat counter
for sess_i = 1:length(s)
    keep = true;

    patID = s(sess_i).pat_ID;
    sessID = s(sess_i).session;

    folder_name = fullfile(eeg_path,patID, sessID);

    if ~isfile(fullfile(folder_name, "allwave.mat"))
        fprintf('Patient %s,%s does not have the required files \n', patID,sessID);
        continue;
    end

    fprintf("Processing patient %s, %s \n", patID, sessID);
    
    %ks added- note which subjects have allwave but no sess_struct_pc
    if ~isfile(fullfile(folder_name, "sess_struct_pc.mat"))
        cnt=cnt+1;
        missingDat.sub{cnt}=patID;
        missingDat.sess{cnt}=sessID;
        missingDat.sIx(cnt)=sess_i;
        continue;
    else
    load(fullfile(folder_name, 'sess_struct_pc'));
    end

    %% Check if cued RT is faster than uncued (correct)

    correct = [sess_struct.evs_pres.responseCorrect]; %get correct responses
    correct(isnan(correct)) = 0; %remove NaN values

    asterisk = [strcmp({sess_struct.evs_pres.asteriskStr},'*before')]; %get cued responses
    no_asterisk = [strcmp({sess_struct.evs_pres.asteriskStr},'*none')]; %get uncued responses

    cued_rt = find(asterisk & correct); %get overlap
    
    RT_cued = NaN(length(cued_rt),1);
    for i = 1:length(cued_rt)
        index = cued_rt(i);
        RT_cued(i) = sess_struct.evs_pres(index).RT; 
    end

    % Find non asterisk events that were correct
    uncued_rt = find(no_asterisk & correct);
    
    RT_uncued = NaN(length(uncued_rt),1);
    for i = 1:length(uncued_rt)
        index = uncued_rt(i);
        RT_uncued(i) = sess_struct.evs_pres(index).RT; 
    end

    avg_cuedRT = mean(RT_cued);
    avg_uncuedRT = mean(RT_uncued);

    if avg_cuedRT > avg_uncuedRT
        keep = false;
        out = sprintf("Uncued faster than cued, %d uncued, %d cued", avg_uncuedRT, avg_cuedRT);
        disp(out);
    end



    %% Check cued vs uncued accuracy
    cued_total = find(asterisk);
    uncued_total = find(no_asterisk);

    cued_accuracy = 100*(length(cued_rt)/length(cued_total));
    uncued_accuracy = 100*(length(uncued_rt)/length(uncued_total));

    if uncued_accuracy > cued_accuracy
        keep = false;
        out = sprintf("Uncued more accurate than cued, %f uncued and %f cued", uncued_accuracy, cued_accuracy);
        disp(out);
    end

    

    %% Check FOIL vs uncued accuracy 
    foil = [strcmp({sess_struct.evs_recog.asteriskStr},'foil')];
    c = [sess_struct.evs_recog.responseCorrect]; %get correct responses for foil
    c(isnan(c)) = 0; %remove NaN values

    foil_total = find(foil);
    foil_yes = find(foil & c);
    foil_accuracy = (1 - length(foil_yes)/length(foil_total))*100;
    foil_yes = (length(foil_yes)/length(foil_total))*100;

    if foil_accuracy > uncued_accuracy
        keep = false;
        disp("More false positives than true positives");
        out = sprintf("Foil accuracy is %f while uncued is %f", foil_accuracy, uncued_accuracy);
        disp(out);
    end

    if keep == true
        soc(soc_i) = sess_i;
        c_accuracy(soc_i) = cued_accuracy;
        n_accuracy(soc_i) = uncued_accuracy;
        f_accuracy(soc_i) = foil_yes;
        soc_i=soc_i+1;
        out = sprintf("Adding patient %s, %s to soc", patID, sessID);
        disp(out)
    end
    if keep == false
        out = sprintf("Discarding patient %s, %s", patID, sessID);
        disp(out)
    end
end

% save('C:\Users\rahil\Box\AttentionWaves\Scripts\soc.mat','soc');
save('/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Scripts/soc.mat','soc');

c_mean = mean(c_accuracy);
n_mean = mean(n_accuracy);
f_mean = mean(f_accuracy);
fprintf(" Avg cued accuracy: %d \n Avg not cued accuracy: %d \n Avg FOIL '% said': %d \n",c_mean, n_mean, f_mean);

save('/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Scripts/cued_accuracy.mat','c_accuracy');
save('/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Scripts/n_cued_accuracy.mat','n_accuracy');
save('/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Scripts/foil_accuracy.mat','f_accuracy');

% save('C:\Users\rahil\Box\AttentionWaves\Scripts\cued_accuracy.mat','c_accuracy');
% save('C:\Users\rahil\Box\AttentionWaves\Scripts\n_cued_accuracy.mat','n_accuracy');
% save('C:\Users\rahil\Box\AttentionWaves\Scripts\foil_accuracy.mat','f_accuracy');

%save log of subj with "allwaves" but not sess_struct_pc.mat
% save('/Users/sundbykk/Library/CloudStorage/Box-Box/AttentionWaves/Scripts/missingDat_Log.mat','missingDat')