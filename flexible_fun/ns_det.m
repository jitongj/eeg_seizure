function ns_det(patient,seizure_label01,eeg_alllabel01,eegseizure01,eeg_alldata,win_num,win_length,wname)


% Select patient =1
%patient=k; % patient     
% summary01=sprintf('chb%02d-summary.txt',patient);
% seizure_label01=label_seizure(summary01);
% eeg_alllabel01 = label_alldata(summary01);
% % Load Data
% eegseizure01=read_seizure(patient);
% eeg_alldata=read_alleeg(patient);
 % Deleting unused and empty channels
[~, eeg_alldata]=delete_chan(eegseizure01, eeg_alldata, patient);
all_file = [eeg_alllabel01(:).chan];

%% find seizure and its front and back position
b = [];%position i.e. row
for k = 1:size(seizure_label01,1)
    seizure_file = seizure_label01(k,1);
    sei_position = find(all_file == seizure_file); 
    if sei_position - 1 > 0
        index = all_file(sei_position - 1);
        if index == seizure_file - 1 % if the files are cont, should include the fornt file
            b(end+1) = sei_position - 1;
        end
    end
    b(end+1) = sei_position;
    if sei_position + 1 <= length(all_file)
        index = all_file(sei_position + 1);
        if index == seizure_file + 1 % if the files are cont, should include the fornt file
            b(end+1) = sei_position + 1;
        end
    end
end

c = unique(b);%position i.e. row
%% find the <1h file's position
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
%%
d = [1:length(all_file)];%position i.e. row
for i = 1:length(c)
    d(find(d==c(i))) = [];
end
%%
k = randperm(length(d),1);
eeg_file = eeg_alldata(patient).patient.data(k).eeg_data;

%%
a = size(eeg_file);
total_length = a(1,2); %point
counter=1;

for i=patient:patient
    for m = 1:win_num
        n = abs(m - win_num-1);
        if m == 1
           w(counter) = randperm(total_length - n*256*win_length,1); %point
           start_point = w(counter); %point
           end_point = start_point+256*win_length-1; %point
           eeg01(counter).window=eeg_alldata(i).patient.data(k).eeg_data(:,start_point:end_point);
        else
            Length = total_length - end_point;
            w(counter) = randperm(Length - n*256*win_length,1); 
            start_point = w(counter) + end_point; %point
            end_point = start_point+256*win_length-1;
            eeg01(counter).window=eeg_alldata(i).patient.data(k).eeg_data(:,start_point:end_point);
        end
        counter=counter+1;
    end
end

%%
Tensor_eeg01 = {};
for i =1: size(eeg01,2)
    mean1_eeg01 = mean(eeg01(i).window(1,:));
    std1_eeg01 = std(eeg01(i).window(1,:));  
    z01_1 = abs(cwt((eeg01(i).window(1,:)-mean1_eeg01)/std1_eeg01,wname)); 
    time01 = size(eeg01(i).window(1,:),2);
    freq01 = size(z01_1,1);
    chan01 = size(eeg01(1).window,1);
    tensor_eeg01 = NaN(chan01, time01, freq01);
    for n= 1:chan01
        mean_eeg01 = mean(eeg01(i).window(n,:));
        std_eeg01 = std(eeg01(i).window(n,:));
        if std_eeg01 == 0
            z01= abs(cwt(eeg01(i).window(n,:)-mean_eeg01,wname));
        else
            z01= abs(cwt((eeg01(i).window(n,:)-mean_eeg01)/std_eeg01,wname));
        end
        for j = 1:time01
            for k = 1:freq01
                tensor_eeg01(n,j,k) = z01(k, j);
            end
        end
    end
    Tensor_eeg01(i,1) = {tensor_eeg01};
end

filename=sprintf('NS%02d',patient);
save(filename,'Tensor_eeg01');


end