function s_det(patient,seizure_label01,eegseizure01,eeg_alldata,win_num,win_length,wname)

% Win_num*win_length <=20, win_length*256 must be an integral
%win_num = 20;
%win_length = 1; % in sec
% Select patient =1
%patient=k; % patient                                                                        % Patienten Nummer                                                                      % edf file           
% summary01=sprintf('chb%02d-summary.txt',patient);
% seizure_label01=label_seizure(summary01);
% %% Load Data
% eegseizure01=read_seizure(patient);
% eeg_data01=read_eeg(patient);
 %Deleting unused and empty channels
[eegseizure01, ~]=delete_chan(eegseizure01, eeg_alldata, patient);
%% choose the seizure period
%calculate seizure length
for i = 1:size(seizure_label01,1)
     seizure_length(i,1)= seizure_label01(i,3)-seizure_label01(i,2);
end

k = randperm(size(seizure_label01,1),1); %row in siezure label matrix
while seizure_length(k,1) < win_num*win_length 
    k = randperm(size(seizure_label01,1),1);
end
    
%%
a = size(eegseizure01(patient).patient.data(k).seizure_data);
total_length = a(1,2); %point
counter=1;

for i=patient:patient
    for m = 1:win_num
        n = abs(m - win_num - 1);
        if m == 1
           w(counter) = randperm(total_length - n*256*win_length,1); %point
           start_point = w(counter); %point
           end_point = start_point+256*win_length-1; %point
           seizure01(counter).window=eegseizure01(i).patient.data(k).seizure_data(:,start_point:end_point);
        else
            Length = total_length - end_point;
            w(counter) = randperm(Length - n*256*win_length,1); 
            start_point = w(counter) + end_point; %point
            end_point = start_point+256*win_length-1;
            seizure01(counter).window=eegseizure01(i).patient.data(k).seizure_data(:,start_point:end_point);
        end
        counter=counter+1;
    end
end


%%
Tensor_Seizure01 = {};
for i =1: size(seizure01,2)
    mean1_seiz01 = mean(seizure01(i).window(1,:));
    std1_seiz01 = std(seizure01(i).window(1,:));  
    z01_1 = abs(cwt((seizure01(i).window(1,:)-mean1_seiz01)/std1_seiz01,wname)); 
    time01 = size(seizure01(i).window(1,:),2);
    freq01 = size(z01_1,1);
    chan01 = size(seizure01(1).window,1);
    tensor_seizure01 = NaN(chan01, time01, freq01);
    for n= 1:chan01
        mean_seiz01 = mean(seizure01(i).window(n,:));
        std_seiz01 = std(seizure01(i).window(n,:));
        if std_seiz01 == 0
            z01= abs(cwt(seizure01(i).window(n,:)-mean_seiz01,wname));
        else
            z01= abs(cwt((seizure01(i).window(n,:)-mean_seiz01)/std_seiz01,wname));
        end 
        for j = 1:time01
            for k = 1:freq01
                tensor_seizure01(n,j,k) = z01(k, j);
            end
        end
    end
    Tensor_Seizure01(i,1) = {tensor_seizure01};
end               

filename=sprintf('S%02d',patient);
save(filename,'Tensor_Seizure01');


end
