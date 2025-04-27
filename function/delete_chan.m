function [eegseizure eeg_data]=delete_chan(eegseizure, eeg_data, patient)
    if (patient==1 || patient==2 || patient==3 || patient==4 || patient==5 || patient==6 || patient==7 || patient==8 || patient==9 || patient==10 || patient==11 || patient==14 || patient==20  || patient==21  || patient==22  || patient==23)
        
        if patient==4 || patient==9
            loeschen=[24];
            if patient==4
                starting_file_eeg=7;
                starting_file_seizure=2;
            else
                starting_file_eeg=2;
                starting_file_seizure=1;
            end
            
            for i=starting_file_seizure:length(eegseizure(patient).patient.data) 
                for k=length(loeschen):-1:1
                    eegseizure(patient).patient.data(i).seizure_data(loeschen(k),:)=[];
                end 
                % eeg_data(patient).patient.data(i).eeg_data(24:end,:)=[];
            end 

            for i=starting_file_eeg:length(eeg_data(patient).patient.data)
                for k=length(loeschen):-1:1
                    eeg_data(patient).patient.data(i).eeg_data(loeschen(k),:)=[];
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
            
            for i=starting_file_seizure:length(eegseizure(patient).patient.data) 
                for k=length(loeschen):-1:1
                    eegseizure(patient).patient.data(i).seizure_data(loeschen(k),:)=[];
                end 
                
                %change the order
                A = eegseizure(patient).patient.data(i).seizure_data;
                eegseizure(patient).patient.data(i).seizure_data(9:16,:) = A(11:18,:);
                eegseizure(patient).patient.data(i).seizure_data(17:18,:) = A(9:10,:);
                  
            end 

            for i=starting_file_eeg:length(eeg_data(patient).patient.data)
                for k=length(loeschen):-1:1
                    eeg_data(patient).patient.data(i).eeg_data(loeschen(k),:)=[];
                end 
               
                % change the order
                B = eeg_data(patient).patient.data(i).eeg_data;
                eeg_data(patient).patient.data(i).eeg_data(9:16,:) = B(11:18,:);
                eeg_data(patient).patient.data(i).eeg_data(17:18,:) = B(9:10,:);
  
            end  
 
        end     
    end
end
