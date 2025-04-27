list=[1,2,3,4,5,6,7,8,9,10,11,14,20,21,22,23];
% win_num*win_length <=20, win_length*256 must be an integer
win_num = 20;
win_length = 1; % in sec

% pre_length <= 298 (pre_win < pre_length)
% pre_length <=94 ||  pre_length - pre_win >= 101, pre_length <=130 || pre_length - pre_win >= 138
pre_length = 3.5*60; % in sec
pre_win = 0.5*60; % in sec

% wname is the analytic wavelet used to compute the CWT. 
% Valid options for wname are "morse", "amor", and "bump".
wname = 'morse';

% rate must be a divisor of 256.
rate = 4;

for patient=list
    summary=sprintf('chb%02d-summary.txt',patient);
    seizure_label=label_seizure(summary);
    all_data_label = label_alldata(summary);
    seizure=read_seizure(patient);
    all_data=read_alldata(patient);
    seizure_detection_tensor(patient,seizure_label,seizure,all_data,win_num,win_length,wname);
    nonseizure_detection_tensor(patient,seizure_label,all_data_label,seizure,all_data,win_num,win_length,wname,rate);
    seizure_prediction_tensor(patient,seizure_label,all_data_label,seizure,all_data,pre_length,pre_win,wname,rate);
    nonseizure_prediction_tensor(patient,seizure_label,all_data_label,seizure,all_data,pre_win,wname,rate);
end
