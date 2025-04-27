function label= label_alldata(filename)
text_data=fopen(filename);
chb_data=textscan(text_data, '%s','delimiter','\t');
fclose(text_data);
channel_number=strfind(chb_data{1,1},'File Name:');
index_file_all=find(not(cellfun('isempty',channel_number)));
label = struct('chan',{},'start',{},'end',{});
for i=1:length(index_file_all)
    num = str2num(regexprep(chb_data{1,1}{index_file_all(i)}, '\D+', ' '));
    num_on=str2num(regexprep(chb_data{1,1}{index_file_all(i)+1}, '\D+', ' '));
    num_off=str2num(regexprep(chb_data{1,1}{index_file_all(i)+2}, '\D+', ' '));
    label(i).start=num_on;
    label(i).end=num_off;
    label(i).chan=num(2);
end 
clear num; clear num_on; clear num_off;
end