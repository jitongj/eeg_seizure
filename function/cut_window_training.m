%% Fensterung der Seizure Abschnitte in 1sek 
function [eeg seizure]=cut_window_training(patient, eegseizure,eeg_data,eeg_size, seizure_size)
counter=1;
if patient==15  %% Overlapping for seizure files 
    ol=10;
else 
    ol=1;
end 
%% Ictal sequences with a new window every sample
for i=patient:patient
%      for k=length(eegseizure(i).patient.data)-seizure_size+1:length(eegseizure(i).patient.data)
         for k=1:seizure_size
        for w=1:ol:length(eegseizure(i).patient.data(k).seizure_data)-256     
            seizure(counter).window=eegseizure(i).patient.data(k).seizure_data(1:23,(w):(w+255));
            counter=counter+1;
        end 

     end 
end 
%% Interictal with a new window every 256samples (1sek)
counter=1;
for i=patient:patient
%     for k=length(eeg_data(i).patient.data)-eeg_size+1:length(eeg_data(i).patient.data) 
       for k=1:eeg_size      
        for w=1:floor(length(eeg_data(i).patient.data(k).eeg_data)/256 -1)
            eeg(counter).window=eeg_data(i).patient.data(k).eeg_data(1:23,(w*256-255):(w*256)); 
            counter=counter+1;
        end 
    end 
end 

end 