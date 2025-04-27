function nonseizure_prediction_window(patient,seizure_label,all_data_label,seizure,all_data,pre_win)
 %% Delete unused and empty channels
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
%Remove duplicated positions
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

% Randomly pick the same amount as the seizures' amount for prediction windows in the nonseizure data.
a = size(eeg_file);
total_length = a(1,2); %point
counter=1;

number = size(seizure_label,1);


for i=patient:patient
    for m = 1:number
        n = abs(m - number - 1);
        if m == 1
           w(counter) = randperm(total_length - n*pre_win*256,1); %point 30sec each 
           start_point = w(counter); %point
           end_point = start_point+pre_win*256-1; %point
           nonsz_pred_win(counter).window=all_data(i).patient.data(k).eeg_data(:,start_point:end_point);
        else
            Length = total_length - end_point;
            w(counter) = randperm(Length - n*pre_win*256,1); 
            start_point = w(counter) + end_point; %point
            end_point = start_point+pre_win*256-1;
            nonsz_pred_win(counter).window=all_data(i).patient.data(k).eeg_data(:,start_point:end_point);
        end
        counter=counter+1;
    end
end


filename=sprintf('NonSz_pre_win%02d',patient);
save(filename,'nonsz_pred_win');

end