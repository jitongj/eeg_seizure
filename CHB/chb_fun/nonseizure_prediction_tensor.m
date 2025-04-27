function nonseizure_prediction_tensor(patient,seizure_label,all_data_label,seizure,all_data,pre_win,wname,rate)
 %% Deleting unused and empty channels
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

% Remove duplicate positions
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
counter1=1;

number = size(seizure_label,1);


for i=patient:patient
    for m = 1:number
        n = abs(m - number - 1);
        if m == 1
           w(counter1) = randperm(total_length - n*pre_win*256,1); %point 30sec each 
           start_point = w(counter1); %point
           end_point = start_point+pre_win*256-1; %point
           nonsz_pre_win(counter1).window=all_data(i).patient.data(k).eeg_data(:,start_point:end_point);
        else
            Length = total_length - end_point;
            w(counter1) = randperm(Length - n*pre_win*256,1); 
            start_point = w(counter1) + end_point; %point
            end_point = start_point+pre_win*256-1;
            nonsz_pre_win(counter1).window=all_data(i).patient.data(k).eeg_data(:,start_point:end_point);
        end
        counter1=counter1+1;
    end
end

%% Build seizure detection tensor

NonSz_pre_tensor = {};
for m =1: size(nonsz_pre_win,2)
    
    % Cut off windows to 30 1sec windows
    counter2=1;
    for w=1:size(nonsz_pre_win(m).window,2)/256
        nonsz_1sec(counter2).window=nonsz_pre_win(m).window(:,(w*256-255):(w*256));
        counter2=counter2+1;
    end 
    
    % Build 3D tensor for each 1sec winow 
    Tensor_nonsz = {};
    for i =1: size(nonsz_1sec,2)

        % Try with one sample to get the size of 3d tensor
        mean1_nonsz = mean(nonsz_1sec(i).window(1,:));
        std1_nonsz = std(nonsz_1sec(i).window(1,:));  
        z1 = abs(cwt((nonsz_1sec(i).window(1,:)-mean1_nonsz)/std1_nonsz,wname)); 
        time = size(nonsz_1sec(i).window,2)/rate;
        freq = size(z1,1);
        chan = size(nonsz_1sec(1).window,1);

        % Build the empty 3D tensor
        for n= 1:chan
            mean_nonsz = mean(nonsz_1sec(i).window(n,:));
            std_nonsz = std(nonsz_1sec(i).window(n,:));
            if std_nonsz == 0
                z= abs(cwt(nonsz_1sec(i).window(n,:)-mean_nonsz,wname));
            else
                z= abs(cwt((nonsz_1sec(i).window(n,:)-mean_nonsz)/std_nonsz,wname));
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
        Tensor_nonsz(i,1) = {tensor_nonseizure};
    end  
    
    % Build the 4D tensor with the 30 1sec windows in the prediction window
    chan = size(tensor_nonseizure,1);
    time = size(tensor_nonseizure,2);
    freq = size(tensor_nonseizure,3);
    Tensor4D_nonsz = NaN(size(nonsz_1sec,2), chan, time, freq);
    for l = 1:size(nonsz_1sec,2)
        for i = 1:chan
            for j = 1:time
                for k = 1:freq
                    Tensor4D_nonsz (l,i,j,k) = Tensor_nonsz{l,1}(i,j,k);
                
                end
            end
        end
    end
    % Build the cell of all 4D tensor 
    NonSz_pre_tensor(m,1) = {Tensor4D_nonsz};

        
end

filename=sprintf('NonSz_pre_tensor%02d',patient);
save(filename,'NonSz_pre_tensor');

end
