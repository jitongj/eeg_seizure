function ns_pred(patient,seizure_label01,eeg_alllabel01,eegseizure01,eeg_alldata,pre_win,wname)

% Select patient =1
%patient=k; % patient     
% summary01=sprintf('chb%02d-summary.txt',patient);
% seizure_label01=label_seizure(summary01);
% eeg_alllabel01 = label_alldata(summary01);
% % Load Data
% eegseizure01=read_seizure(patient);
% eeg_alldata=read_alleeg(patient);
 %% Deleting unused and empty channels
[~, eeg_alldata]=delete_chan(eegseizure01, eeg_alldata, patient);
all_file = [eeg_alllabel01(:).chan];

%% find all seizure file's position, as well as its continous front and back files' position.
b = [];%position i.e. row
for k = 1:size(seizure_label01,1)
    seizure_file = seizure_label01(k,1);
    sei_position = find(all_file == seizure_file); 
    if sei_position - 1 > 0
        index = all_file(sei_position - 1);
        if index == seizure_file - 1 % if the files are continous, should include the fornt file
            b(end+1) = sei_position - 1;
        end
    end
    b(end+1) = sei_position;
    if sei_position + 1 <= length(all_file)
        index = all_file(sei_position + 1);
        if index == seizure_file + 1 % if the files are continous, should include the fornt file
            b(end+1) = sei_position + 1;
        end
    end
end

%Remove duplicate positions
c = unique(b);%position i.e. row

%% find the position of the file which is less than 1h. If these short files happen to be the previous or next file in the seizure file, we need to include one more file
short = [];%position i.e. row
for i = 1:length(all_file)
    file_start = eeg_alllabel01(i).start;
    file_end = eeg_alllabel01(i).end;
    each_file_length = file_end(1)*60*60+file_end(2)*60+file_end(3)-(file_start(1)*60*60+file_start(2)*60+file_start(3));
    if each_file_length < 3600
        short(end+1) = i;
    end
end


if ~isempty(short)
    for i = length(short)
        if ~(ismember(all_file(short(i)),seizure_label01(:,1)))
            if ismember(all_file(short(i)-1),seizure_label01(:,1))
                c(end+1) = short(i)+1;
            elseif ismember(all_file(short(i)+1),seizure_label01(:,1))
                c(end+1) = short(i)-1;
            end
        end
    end
end
c = unique(c);%position i.e. row

%%Delete the seizure files and the files that are at least an hour before and after it, then the remaining eeg files are all at least 1 hour away from any seizure.
d = [1:length(all_file)];%position i.e. row
for i = 1:length(c)
    d(find(d==c(i))) = [];
end


%% Randomly pick one eeg file from the remaining
k = randperm(length(d),1);
eeg_file = eeg_alldata(patient).patient.data(k).eeg_data;

%%Randomly pick the same amount as the seizures' amount for prediction windows in this file
a = size(eeg_file);
total_length = a(1,2); %point
counter=1;

number = size(seizure_label01,1);


for i=patient:patient
    for m = 1:number
        n = abs(m - number - 1);
        if m == 1
           w(counter) = randperm(total_length - n*pre_win*256,1); %point 30sec each 
           start_point = w(counter); %point
           end_point = start_point+pre_win*256-1; %point
           eeg01(counter).window=eeg_alldata(i).patient.data(k).eeg_data(:,start_point:end_point);
        else
            Length = total_length - end_point;
            w(counter) = randperm(Length - n*pre_win*256,1); 
            start_point = w(counter) + end_point; %point
            end_point = start_point+pre_win*256-1;
            eeg01(counter).window=eeg_alldata(i).patient.data(k).eeg_data(:,start_point:end_point);
        end
        counter=counter+1;
    end
end

%%
Tensor4D_pre_eeg = {};
for m =1: size(eeg01,2)
    
    %cut off windows to 30 1sec windows
    counter01_2=1;
    for w=1:size(eeg01(m).window,2)/256
        eeg01_1sec(counter01_2).window=eeg01(m).window(:,(w*256-255):(w*256));
        counter01_2=counter01_2+1;
    end 
    
    %build 3D tensor
    Tensor_eeg01 = {};
    for i =1: size(eeg01_1sec,2)
        mean1_eeg01 = mean(eeg01_1sec(i).window(1,:));
        std1_eeg01 = std(eeg01_1sec(i).window(1,:));  
        z01_1 = abs(cwt((eeg01_1sec(i).window(1,:)-mean1_eeg01)/std1_eeg01,wname)); 
        time01 = size(eeg01_1sec(i).window,2);
        freq01 = size(z01_1,1);
        chan01 = size(eeg01_1sec(1).window,1);
        tensor_eeg01 = NaN(chan01, time01, freq01);
        for n= 1:chan01
            mean_eeg01 = mean(eeg01_1sec(i).window(n,:));
            std_eeg01 = std(eeg01_1sec(i).window(n,:));
            if std_eeg01 == 0
                z01= abs(cwt(eeg01_1sec(i).window(n,:)-mean_eeg01,wname));
            else
                z01= abs(cwt((eeg01_1sec(i).window(n,:)-mean_eeg01)/std_eeg01,wname));
            end
            for j = 1:time01
                for k = 1:freq01
                    tensor_eeg01(n,j,k) = z01(k, j);
                end
            end
        end
        Tensor_eeg01(i,1) = {tensor_eeg01};
    end  
    
    Tensor4D_eeg = NaN(size(eeg01_1sec,2), chan01, time01, freq01);
    for l = 1:size(eeg01_1sec,2)
        for i = 1:chan01
            for j = 1:time01
                for k = 1:freq01
                    Tensor4D_eeg (l,i,j,k) = Tensor_eeg01{l,1}(i,j,k);
                
                end
            end
        end
    end
    
    Tensor4D_pre_eeg(m,1) = {Tensor4D_eeg};

        
end

filename=sprintf('NSpre%02d',patient);
save(filename,'Tensor4D_pre_eeg');

end