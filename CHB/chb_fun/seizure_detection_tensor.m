function seizure_detection_tensor(patient,seizure_label,seizure,all_data,win_num,win_length,wname,rate)

%% Delete unsued and empty channels.
[seizure, ~]=delete_chan(seizure, all_data, patient);
%% Randomly pick a seizure.
% Calculate seizure length
for i = 1:size(seizure_label,1)
     seizure_length(i,1)= seizure_label(i,3)-seizure_label(i,2);
end

% Randomly pick a seizure
k = randperm(size(seizure_label,1),1); %row in siezure label matrix

% If the seizure length is smaller than the total windows' length, randomly, pick another one again
while seizure_length(k,1) < win_num*win_length 
    k = randperm(size(seizure_label,1),1);
end
    
%% Cut windows from the seizure data.
a = size(seizure(patient).patient.data(k).seizure_data);
total_length = a(1,2); %point
counter=1;

for i=patient:patient
    for m = 1:win_num
        n = abs(m - win_num - 1);
        if m == 1
           w(counter) = randperm(total_length - n*256*win_length,1); %point
           start_point = w(counter); %point
           end_point = start_point+256*win_length-1; %point
           sz_det_win(counter).window=seizure(i).patient.data(k).seizure_data(:,start_point:end_point);
        else
            Length = total_length - end_point;
            w(counter) = randperm(Length - n*256*win_length,1); 
            start_point = w(counter) + end_point; %point
            end_point = start_point+256*win_length-1;
            sz_det_win(counter).window=seizure(i).patient.data(k).seizure_data(:,start_point:end_point);
        end
        counter=counter+1;
    end
end


%% Build seizure detection tensor
Sz_det_tensor = {};
for i =1: size(sz_det_win,2)

    % Try with one sample to get the size of tensor
    mean1_sz = mean(sz_det_win(i).window(1,:));
    std1_sz = std(sz_det_win(i).window(1,:));  
    z1 = abs(cwt((sz_det_win(i).window(1,:)-mean1_sz)/std1_sz,wname));

    time = size(sz_det_win(i).window(1,:),2)/rate;
    freq = size(z1,1);
    chan = size(sz_det_win(1).window,1);

    % Build the 3d tensor
    for n= 1:chan

        % Normalize and do the cwt
        mean_sz = mean(sz_det_win(i).window(n,:));
        std_sz = std(sz_det_win(i).window(n,:));
        if std_sz == 0
            z= abs(cwt(sz_det_win(i).window(n,:)-mean_sz,wname));
        else
            z= abs(cwt((sz_det_win(i).window(n,:)-mean_sz)/std_sz,wname));
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
            for k = 1:freq
                tensor_seizure(n,j,k) = z(k, j);
            end
        end
    end
    Sz_det_tensor(i,1) = {tensor_seizure};
end               

filename=sprintf('Sz_det_tensor%02d',patient);
save(filename,'Sz_det_tensor');


end
