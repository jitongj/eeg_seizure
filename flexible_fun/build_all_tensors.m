list=[1,2,3,4,5,6,7,8,9,10,11,14,20,21,22,23];
% win_num*win_length <=20, win_length*256 must be an integral
win_num = 20;
win_length = 1; % in sec

% pre_length <= 298 (pre_win < pre_length)
% pre_length <=94 ||  pre_length - pre_win >= 101, pre_length <=130 || pre_length - pre_win >= 138
pre_length = 3.5*60; % in sec
pre_win = 0.5*60; % in sec

% wname is the analytic wavelet used to compute the CWT. 
% Valid options for wname are "morse", "amor", and "bump".
wname = 'morse';


for patient=list
    summary01=sprintf('chb%02d-summary.txt',patient);
    seizure_label01=label_seizure(summary01);
    eeg_alllabel01 = label_alldata(summary01);
    eegseizure01=read_seizure(patient);
    eeg_alldata=read_alleeg(patient);
    s_det(patient,seizure_label01,eegseizure01,eeg_alldata,win_num,win_length,wname);
    ns_det(patient,seizure_label01,eeg_alllabel01,eegseizure01,eeg_alldata,win_num,win_length,wname);
    ns_pred(patient,seizure_label01,eeg_alllabel01,eegseizure01,eeg_alldata,pre_win,wname);
    s_pre(patient,seizure_label01,eeg_alllabel01,eegseizure01,eeg_alldata,pre_length,pre_win,wname);
end
