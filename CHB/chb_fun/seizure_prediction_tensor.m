function seizure_prediction_tensor(patient,seizure_label,all_data_label,seizure,all_data,pre_length,pre_win,wname,rate)
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
        prediction_start = Length-remain; %3.5mis
        prediction_end = prediction_start + pre_win; %30sec window
        sz_pred_win(k).window = pre_file(:,prediction_start*256: prediction_end*256 - 1);
    else
        pre_file = all_data(patient).patient.data(sei_position).eeg_data;
        
        prediction_start = seizure_label(k,2)-pre_length; %3.5mis
        prediction_end = prediction_start + pre_win; %30sec window  
        sz_pred_win(k).window = pre_file(:,prediction_start*256: prediction_end*256 - 1);
    end
        
end

%% Build the cell of 4D tensors

Sz_pre_tensor = {};
for m =1: size(sz_pred_win,2)
    
    % Cut off windows to 30 1sec windows
    counter=1;
    for w=1:size(sz_pred_win(m).window,2)/256
        sz_1sec(counter).window=sz_pred_win(m).window(:,(w*256-255):(w*256));
        counter=counter+1;
    end 
    
    % Build 3D tensor for each 1sec winow 
    Tensor_sz = {};
    for i =1: size(sz_1sec,2)
        
        % Try with one sample to get the size of 3d tensor
        mean1_seizure = mean(sz_1sec(i).window(1,:));
        std1_seizure = std(sz_1sec(i).window(1,:));  
        z1 = abs(cwt((sz_1sec(i).window(1,:)-mean1_seizure)/std1_seizure,wname)); 
        time = size(sz_1sec(i).window,2)/rate;
        freq = size(z1,1);
        chan = size(sz_1sec(i).window,1);
        
        % Build the empty 3D tensor
        for n= 1:chan
            mean_sz = mean(sz_1sec(i).window(n,:));
            std_sz = std(sz_1sec(i).window(n,:));
            if std_sz == 0
                z= abs(cwt(sz_1sec(i).window(n,:)-mean_sz,wname));
            else
                z= abs(cwt((sz_1sec(i).window(n,:)-mean_sz)/std_sz,wname));
            end

            % Downsample the time-freq domain
            z = downsample(z.',rate);
            z = z.';
            
            % Build the empty 3D tensor
            time = size(z,2);
            if n = 1
            tensor_seizure = NaN(chan, time, freq);
            end

            % Fill the tensor
            for j = 1:time
                for k = 1:51
                    tensor_seizure(n,j,k) = z(k, j);
                    
                end
            end
        end
        Tensor_sz(i,1) = {tensor_seizure};
    end  
    
    % Build the 4D tensor with the 30 1sec windows in the prediction window
    chan = size(tensor_seizure,1);
    time = size(tensor_seizure,2);
    freq = size(tensor_seizure,3);
    Tensor4D_sz = NaN(size(sz_1sec,2), chan, time, freq);
    for l = 1:size(sz_1sec,2)
        for i = 1:chan
            for j = 1:time
                for k = 1:freq
                    Tensor4D_sz (l,i,j,k) = Tensor_sz{l,1}(i,j,k);
                
                end
            end
        end
    end
    
    % Build the cell of all 4D tensor
    Sz_pre_tensor(m,1) = {Tensor4D_sz};

        
end

filename=sprintf('Sz_pre_tensor%02d',patient);
save(filename,'Sz_pre_tensor');


end
