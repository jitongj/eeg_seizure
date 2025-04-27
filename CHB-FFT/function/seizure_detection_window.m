function seizure_detection_window(patient,seizure_label,seizure,all_data,win_num,win_length)

%% Delete unsued and empty channels.
[seizure, ~]=delete_chan(seizure, all_data, patient);
%% Randomly pick a seizure.

% Calculate all seizures' length.
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

filename=sprintf('Sz_det_win%02d',patient);
save(filename,'sz_det_win');


end