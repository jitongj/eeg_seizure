%% These are the basic settings that can be changed within the range in description
% win_num*win_length <=20, win_length*256 must be an integer
win_num = 20;
win_length = 1; % in sec

% pre_length <= 298 (pre_win < pre_length)
% pre_length <=94 ||  pre_length - pre_win >= 101, pre_length <=130 || pre_length - pre_win >= 138
pre_length = 3.5*60; % in sec
pre_win = 0.5*60; % in sec

%% Create datasets for detection and prediction windows of seizure and nonseizure
list=[1,2,3,4,5,6,7,8,9,10,11,14,20,21,22,23];
for patient=list
    summary=sprintf('chb%02d-summary.txt',patient);
    seizure_label=label_seizure(summary); % find the file number of seizure file
    all_data_label = label_alldata(summary); % find the file number of signal file
    seizure=read_seizure(patient); % read the seizure files
    all_data=read_alldata(patient); % read all signal files

    seizure_detection_window(patient,seizure_label,seizure,all_data,win_num,win_length);
    nonseizure_detection_window(patient,seizure_label,all_data_label,seizure,all_data,win_num,win_length);
    nonseizure_prediction_window(patient,seizure_label,all_data_label,seizure,all_data,pre_win);
    seizure_prediction_window(patient,seizure_label,all_data_label,seizure,all_data,pre_length,pre_win);
end

%% Create datasets for feature detection and prediction of seizure and nonseizure
list=[1,2,3,4,5,6,7,8,9,10,11,14,20,21,22,23];
for patient=list
    % Load prediction and detection data of seizure and nonseizure
    sz_det=sprintf('Sz_det_win%02d',patient);
    nonsz_det=sprintf('NonSz_det_win%02d',patient);
    sz_pre=sprintf('Sz_pre_win%02d',patient);
    nonsz_pre=sprintf('NonSz_pre_win%02d',patient);
    load(sz_det);
    load(nonsz_det);
    load(sz_pre);
    load(nonsz_pre);

    det_feature(patient,win_num,sz_det_win,nonsz_det_win); % contains seizure and nonseizure
    pre_feature(patient,sz_pred_win,nonsz_pred_win); % contains seizure and nonseizure

end






  
