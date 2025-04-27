function eeg_alldata=read_alldata(patient)
    for k=patient:patient
        b=sprintf('chb%02d-summary.txt',k)
        label= label_alldata(b);

        for i=1:size(label,2)
            eeg.info(i,1)=label(i).chan;            
            a=sprintf('chb%02d_%02d.edf',k,label(i).chan);
            [z q]=edfread(a);
            eeg.data(i).eeg_data=q(:,:);
        end
        eeg_alldata(k).patient=eeg;
        clear eeg;
    end 
end