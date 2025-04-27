function [seizure all_data]=delete_chan(seizure, all_data, patient)
  
        if patient==4 || patient==9
            % delete the extra channel 24th
            loeschen=[24];
            if patient==4
                starting_file_eeg=7;
                starting_file_seizure=2;
            else
                starting_file_eeg=2;
                starting_file_seizure=1;
            end
            
            for i=starting_file_seizure:length(seizure(patient).patient.data) 
                for k=length(loeschen):-1:1
                    seizure(patient).patient.data(i).seizure_data(loeschen(k),:)=[];
                end 

            end 

            for i=starting_file_eeg:length(all_data(patient).patient.data)
                for k=length(loeschen):-1:1
                    all_data(patient).patient.data(i).eeg_data(loeschen(k),:)=[];
                end 
                % eeg_data(patient).patient.data(i).eeg_data(24:end,:)=[];
            end  
    
            
            
        elseif (patient==11 || patient==14 || patient==20 || patient==21 || patient==22)

            loeschen=[5 10 13 18 23];
            if patient==11
                starting_file_eeg=2;
                starting_file_seizure=1;
            else
                starting_file_eeg=1;
                starting_file_seizure=1;
            end
            
            for i=starting_file_seizure:length(seizure(patient).patient.data) 
                for k=length(loeschen):-1:1
                    seizure(patient).patient.data(i).seizure_data(loeschen(k),:)=[];
                end 
                
                %change the order of channel
                A = seizure(patient).patient.data(i).seizure_data;
                seizure(patient).patient.data(i).seizure_data(9:16,:) = A(11:18,:);
                seizure(patient).patient.data(i).seizure_data(17:18,:) = A(9:10,:);
                  
            end 

            for i=starting_file_eeg:length(all_data(patient).patient.data)
                for k=length(loeschen):-1:1
                    all_data(patient).patient.data(i).eeg_data(loeschen(k),:)=[];
                end 
               
                % change the order of channel
                B = all_data(patient).patient.data(i).eeg_data;
                all_data(patient).patient.data(i).eeg_data(9:16,:) = B(11:18,:);
                all_data(patient).patient.data(i).eeg_data(17:18,:) = B(9:10,:);
  
            end  
 
        end     
end