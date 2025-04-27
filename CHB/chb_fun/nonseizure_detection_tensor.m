function nonseizure_detection_tensor(patient,seizure_label,all_data_label,seizure,all_data,win_num,win_length,wname,rate)

%% Delete unsued and empty channels
[~, all_data]=delete_chan(seizure, all_data, patient);

%% Read all data's channels
all_file = [all_data_label(:).chan];

%% Find the position of seizure files and all files within one hour before and after them
% Find seizure's position and its previous and next positions and include all of them if they are continuous.
b = [];%position i.e. row
for k = 1:size(seizure_label,1)
    seizure_file = seizure_label(k,1);
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
% Remove duplicated positions
c = unique(b);%position i.e. row
% Find the <1h file's position. 
% If they are the previous or next consecutive files of the seizure files,
% include an additonal file to make sure all the files winthin 1h of seizure files are included.
short = [];%position i.e. row
for i = 1:length(all_file)
    file_start = all_data_label(i).start;
    file_end = all_data_label(i).end;
    each_file_length = file_end(1)*60*60+file_end(2)*60+file_end(3)-(file_start(1)*60*60+file_start(2)*60+file_start(3));
    if each_file_length < 3600
        short(end+1) = i;
    end
end


if ~isempty(short)
    for i = length(short)
        if ~(ismember(all_file(short(i)),seizure_label(:,1)))
            if ismember(all_file(short(i)-1),seizure_label(:,1))
                c(end+1) = short(i)+1;
            elseif ismember(all_file(short(i)+1),seizure_label(:,1))
                c(end+1) = short(i)-1;
            end
        end
    end
end
c = unique(c);%position i.e. row
%% Delete all seizure files and all files within one hour of them
% The remain are pure non seizure data files that are at least 1h away from any seizure.
d = [1:length(all_file)];%position i.e. row
for i = 1:length(c)
    d(find(d==c(i))) = [];
end
%% Randomly pick a nonseizure file.
k = randperm(length(d),1);
eeg_file = all_data(patient).patient.data(k).eeg_data;

%% Cut windows from the nonseizure data.
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
           nonsz_det_win(counter).window=all_data(i).patient.data(k).eeg_data(:,start_point:end_point);
        else
            Length = total_length - end_point;
            w(counter) = randperm(Length - n*256*win_length,1); 
            start_point = w(counter) + end_point; %point
            end_point = start_point+256*win_length-1;
            nonsz_det_win(counter).window=all_data(i).patient.data(k).eeg_data(:,start_point:end_point);
        end
        counter=counter+1;
    end
end

%% Build Nonseizure detection tensor
Nonsz_det_tensor = {};
for i =1: size(nonsz_det_win,2)

    % Try with one sample to get the size of tensor
    mean1_nonsz = mean(nonsz_det_win(i).window(1,:));
    std1_nonsz = std(nonsz_det_win(i).window(1,:));  
    z1 = abs(cwt((nonsz_det_win(i).window(1,:)-mean1_nonsz)/std1_nonsz,wname)); 
    time = size(nonsz_det_win(i).window(1,:),2)/rate;
    freq = size(z1,1);
    chan = size(nonsz_det_win(1).window,1);

    % Build the 3d tensor
    for n= 1:chan

        % Normalize and do the cwt
        mean_nonsz = mean(nonsz_det_win(i).window(n,:));
        std_nonsz = std(nonsz_det_win(i).window(n,:));
        if std_nonsz == 0
            z= abs(cwt(nonsz_det_win(i).window(n,:)-mean_nonsz,wname));
        else
            z= abs(cwt((nonsz_det_win(i).window(n,:)-mean_nonsz)/std_nonsz,wname));
        end

        % Downsample the time-freq domain
        z = downsample(z.',rate);
        z = z.';
        
        % Build the empty 3D tensor
        time = size(z,2);
        if n = 1
        tensor_nonseizure = NaN(chan, time, freq);
        end


        % Fill the tensor
        for j = 1:time
            for k = 1:freq
                tensor_nonseizure(n,j,k) = z(k, j);
            end
        end
    end
    Nonsz_det_tensor(i,1) = {tensor_nonseizure};
end

filename=sprintf('NonSz_det_tensor%02d',patient);
save(filename,'Nonsz_det_tensor');


end
