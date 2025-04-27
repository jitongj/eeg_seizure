function seizure_prediction_window(patient,seizure_label,all_data_label,seizure,all_data,pre_length,pre_win)

%% Delete unused and empty channels
[~, all_data]=delete_chan(seizure, all_data, patient);

%% Read all data's channels
all_file = [all_data_label(:).chan];

%% Get the cutting of windows. If the seizure starts less than pre_length, looking the prediction point in the previous file
for k = 1:size(seizure_label,1)
    seizure_file = seizure_label(k,1);
    seizure_start = seizure_label(k,2);
    sei_position = find(all_file == seizure_file);
    if seizure_start < pre_length
        pre_file = all_data(patient).patient.data(sei_position-1).eeg_data;
        
        % Calculate the time gap between the two fiels.
        min1 = all_data_label(sei_position).start(2);
        min2 = all_data_label(sei_position-1).end(2);
        sec1 = all_data_label(sei_position).start(3);
        sec2 = all_data_label(sei_position-1).end(3);
        time1 = min1*60+sec1; %sec
        time2 = min2*60+sec2; %sec
        gap = time1 - time2; %sec
        
        remain = pre_length - gap - seizure_start; %sec: positive
        
        Length = size(pre_file,2)/256; %time length in sec
        prediction_start = Length-remain; 
        prediction_end = prediction_start + pre_win;
        sz_pred_win(k).window = pre_file(:,prediction_start*256: prediction_end*256 - 1);
    else
        pre_file = all_data(patient).patient.data(sei_position).eeg_data;
        
        prediction_start = seizure_label(k,2)-pre_length; 
        prediction_end = prediction_start + pre_win; 
        sz_pred_win(k).window = pre_file(:,prediction_start*256: prediction_end*256 - 1);
    end
        
end

filename=sprintf('Sz_pre_win%02d',patient);
save(filename,'sz_pred_win');


end