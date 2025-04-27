function s_pre(patient,seizure_label01,eeg_alllabel01,eegseizure01,eeg_alldata,pre_length,pre_win,wname)

%% read all data from the patient

% pre_length <= 298 (pre_win < pre_length)
% pre_length <=94 ||  pre_length - pre_win >= 101, pre_length <=130 || pre_length - pre_win >= 138



%patient=k; % patient                                                                                
% summary01=sprintf('chb%02d-summary.txt',patient);
% seizure_label01=label_seizure(summary01);
% eeg_alllabel01 = label_alldata(summary01);
% 
% eegseizure01=read_seizure(patient);
% eeg_alldata=read_alleeg(patient);
 % Deleting unused and empty channels
[~, eeg_alldata]=delete_chan(eegseizure01, eeg_alldata, patient);
all_file = [eeg_alllabel01(:).chan];

%% get the cutting of windows. if the seizure starts less than 3.5 mins, looking the prediction point in the previous file
for k = 1:size(seizure_label01,1)
    seizure_file = seizure_label01(k,1);
    seizure_start = seizure_label01(k,2);
    sei_position = find(all_file == seizure_file);
    if seizure_start < pre_length
        pre_file = eeg_alldata(patient).patient.data(sei_position-1).eeg_data;
        
        min1 = eeg_alllabel01(sei_position).start(2);
        min2 = eeg_alllabel01(sei_position-1).end(2);
        sec1 = eeg_alllabel01(sei_position).start(3);
        sec2 = eeg_alllabel01(sei_position-1).end(3);
        time1 = min1*60+sec1; %sec
        time2 = min2*60+sec2; %sec
        gap = time1 - time2; %sec
        remain = pre_length - gap - seizure_start; %sec: positive
        
        Length = size(pre_file,2)/256; %time length in sec
        prediction_start = Length-remain; %3.5mis
        prediction_end = prediction_start + pre_win; %30sec window
        prediction(k).window = pre_file(:,prediction_start*256: prediction_end*256 - 1);
    else
        pre_file = eeg_alldata(patient).patient.data(sei_position).eeg_data;
        
        prediction_start = seizure_label01(k,2)-pre_length; %3.5mis
        prediction_end = prediction_start + pre_win; %30sec window  
        prediction(k).window = pre_file(:,prediction_start*256: prediction_end*256 - 1);
    end
        
end

%% build the cell of 4D tensors

Tensor4D_pre_seizure = {};
for m =1: size(prediction,2)
    
    %cut off windows to 30 1sec windows
    counter01_2=1;
    for w=1:size(prediction(m).window,2)/256
        seizure_1sec(counter01_2).window=prediction(m).window(:,(w*256-255):(w*256));
        counter01_2=counter01_2+1;
    end 
    
    %build 3D tensor for each 1sec winow 
    Tensor_seizure01 = {};
    for i =1: size(seizure_1sec,2)
        
        % the following is to use one channel as a sample to get the size
        % in freq of cwt matrix, so that we can build a empty 3D tensor
        mean1_seizure01 = mean(seizure_1sec(i).window(1,:));
        std1_seizure01 = std(seizure_1sec(i).window(1,:));  
        z01_1 = abs(cwt((seizure_1sec(i).window(1,:)-mean1_seizure01)/std1_seizure01,wname)); 
        time01 = size(seizure_1sec(i).window,2);
        freq01 = size(z01_1,1);
        chan01 = size(seizure_1sec(i).window,1);
        
        % build the empty 3D tensor, and then fill it
        tensor_seizure01 = NaN(chan01, time01, freq01);
        for n= 1:chan01
            mean_seizure01 = mean(seizure_1sec(i).window(n,:));
            std_seizure01 = std(seizure_1sec(i).window(n,:));
            if std_seizure01 == 0
                z01= abs(cwt(seizure_1sec(i).window(n,:)-mean_seizure01,wname));
            else
                z01= abs(cwt((seizure_1sec(i).window(n,:)-mean_seizure01)/std_seizure01,wname));
            end
            for j = 1:time01
                for k = 1:51
                    tensor_seizure01(n,j,k) = z01(k, j);
                    
                end
            end
        end
        Tensor_seizure01(i,1) = {tensor_seizure01};
    end  
    
    %build the 4D tensor with the 30 1sec windows in the prediction window
    chan = size(tensor_seizure01,1);
    time = size(tensor_seizure01,2);
    freq = size(tensor_seizure01,3);
    Tensor4D_seizure = NaN(size(seizure_1sec,2), chan, time, freq);
    for l = 1:size(seizure_1sec,2)
        for i = 1:chan
            for j = 1:time
                for k = 1:freq
                    Tensor4D_seizure (l,i,j,k) = Tensor_seizure01{l,1}(i,j,k);
                
                end
            end
        end
    end
    
    %build the cell of all 4D tensor
    Tensor4D_pre_seizure(m,1) = {Tensor4D_seizure};

        
end

filename=sprintf('Spre%02d',patient);
save(filename,'Tensor4D_pre_seizure');


end