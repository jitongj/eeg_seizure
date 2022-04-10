% Seizure Start Time: 2972 seconds
% Seizure End Time: 3053 seconds

% basic infor
data_chb02 = pop_loadset('chb02_16+.set');
chan=23;

%% create seperated 5s tensor without seizure, choose a time period from 1000s to 1005s

% first try with one example to obtain the size of time and frequency after CWT

% get the mean of first channel
mean1_02_swo = mean(double(data_chb02.data(1,256*1000+1:256*1005))); 
% get the std of first channel
std1_02_swo = std(double(data_chb02.data(1,256*1000+1:256*1005)));
% first standardize, then do wavelet transformation 
z1_02 = abs(cwt((double(data_chb02.data(1,256*1000+1:256*1005))-mean1_02_swo)/std1_02_swo));
% get the size of time
time_02 = size(data_chb02.data(1,256*1000+1:256*1005),2);
% get the size of frequency
freq_02 = size(z1_02,1);
 
% build the seperated 5s tensor without seizure
N_chb02 = NaN(chan, time_02, freq_02);
 for i= 1:chan
     % get the mean of each channel
     mean02_swo =  mean(double(data_chb02.data(i,256*1000+1:256*1005)));
     % get the mean of each channel
     std02_swo = std(double(data_chb02.data(i,256*1000+1:256*1005)));
     % standardize, then do wavelet transformation
     z_02 = abs(cwt((double(data_chb02.data(i,256*1000+1:256*1005))-mean02_swo)/std02_swo));
     % build tensor
     for j = 1:time_02
         for k = 1:freq_02
             N_chb02(i, j, k) = z_02(k, j);
         end
     end
 end
N_chb02 = fmt(N_chb02);

% check the low rank

relerr02=[];
con02 =[];
% randomly initialized CPD models with different ranks ranging from 1 to 20 for 5 times
for i = 1:20
    for n = 1:5
        Uhat02 = cpd(N_chb02,i);
        % calculate the relative error
        relerr02(n,i) = frobcpdres(N_chb02, Uhat02)/frob(N_chb02);
        % calculate the core consistency
        con02(n,i) = corcond(N_chb02,Uhat02,0);
    end
end

% plot the core consistency to find the low rank
subplot(2,1,1);
plot(1:20,con02(1,:),'-');
hold on
plot(1:20,con02(2,:),'-');
hold on
plot(1:20,con02(3,:),'-');
hold on
plot(1:20,con02(4,:),'-');
hold on
plot(1:20,con02(5,:),'-')
axis([1,20,-100,100]);

% plot the relative error to find the low rank
subplot(2,1,2); 
plot(1:20,relerr02(1,:),'*')
hold on
plot(1:20,relerr02(2,:),'*')
hold on
plot(1:20,relerr02(3,:),'*')
hold on
plot(1:20,relerr02(4,:),'*')
hold on
plot(1:20,relerr02(5,:),'*')

% result: low rank = 2

% ran CPD 4 times to check the uniqueness of the low rank = 2

% ran randomly the 1st CPD
Uhat02_r1 = cpd(N_chb02,2);
%after decomposition, multiply with the standard deviation to preserve topographic information
for j = 1:chan
    Uhat02_r1{1}(j,1) = Uhat02_r1{1}(j,1)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
    Uhat02_r1{1}(j,2) = Uhat02_r1{1}(j,2)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
end

% ran randomly the 2nd CPD
Uhat02_r2 = cpd(N_chb02,2);
for j = 1:chan
    Uhat02_r2{1}(j,1) = Uhat02_r2{1}(j,1)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
    Uhat02_r2{1}(j,2) = Uhat02_r2{1}(j,2)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
end

% ran randomly the 3rd CPD
Uhat02_r3 = cpd(N_chb02,2);
for j = 1:chan
    Uhat02_r3{1}(j,1) = Uhat02_r3{1}(j,1)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
    Uhat02_r3{1}(j,2) = Uhat02_r3{1}(j,2)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
end

% ran randomly the 4th CPD
Uhat02_r4 = cpd(N_chb02,2);
for j = 1:chan
    Uhat02_r4{1}(j,1) = Uhat02_r4{1}(j,1)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
    Uhat02_r4{1}(j,2) = Uhat02_r4{1}(j,2)*std(double(data_chb02.data(j,256*1000+1:256*1005)));
end

% compare pairwise of the component using Pearson's correlation coefficient
for k = 1:3 %3 mode
    s = size(Uhat02_r1{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r1{k}(i,j),Uhat02_r2{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat02_r1{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat02_r1{k}(i,j),Uhat02_r3{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat02_r1{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat02_r1{k}(i,j),Uhat02_r4{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat02_r2{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r2{k}(i,j),Uhat02_r3{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat02_r2{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r2{k}(i,j),Uhat02_r4{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat02_r3{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r3{k}(i,j),Uhat02_r4{k}(i,j))~= 1
                disp('3&4: solution is not unique');
                break
            end
        end
    end
end   
disp('3&4: solution is unique');

% result: unique

%% create 5s tensor with seizure, choose a time period from 2972s to 2977s

% first try with one example to obtain the size of time and frequency after CWT like above
mean1_02_sw = mean(double(data_chb02.data(1,256*2972+1:256*2977)));
std1_02_sw = std(double(data_chb02.data(1,256*2972+1:256*2977)));
z1_02_se = abs(cwt((double(data_chb02.data(1,256*2972+1:256*2977))-mean1_02_sw)/std1_02_sw));
time_02_se = size(data_chb02.data(1,256*2972+1:256*2977),2);
freq_02_se = size(z1_02_se,1);

% build the seperated 5s tensor with seizure like above
N_chb02_se = NaN(chan, time_02_se, freq_02_se);
 for i= 1:chan
     mean02_sw = mean(double(data_chb02.data(i,256*2972+1:256*2977)));
     std02_sw = std(double(data_chb02.data(i,256*2972+1:256*2977)));
     z_02_sw = abs(cwt((double(data_chb02.data(i,256*2972+1:256*2977))-mean02_sw)/std02_sw));
     for j = 1:time_02_se
         for k = 1:freq_02_se
             N_chb02_se(i, j, k) = z_02_sw(k, j);
         end
     end
 end
N_chb02_se = fmt(N_chb02_se);

% check the low rank with core consistency and relative error like above
% the result is 2
relerr02_se=[];
con02_se =[];
for i = 1:20
    for n = 1:5
        Uhat02_se = cpd(N_chb02_se,i);
        relerr02_se(n,i) = frobcpdres(N_chb02_se, Uhat02_se)/frob(N_chb02_se);
        con02_se(n,i) = corcond(N_chb02_se,Uhat02_se,0);
    end
end

subplot(2,1,1);
plot(1:20,con02_se(1,:),'-');
hold on
plot(1:20,con02_se(2,:),'-');
hold on
plot(1:20,con02_se(3,:),'-');
hold on
plot(1:20,con02_se(4,:),'-');
hold on
plot(1:20,con02_se(5,:),'-')
axis([1,20,-100,100]);

subplot(2,1,2); 
plot(1:20,relerr02_se(1,:),'*')
hold on
plot(1:20,relerr02_se(2,:),'*')
hold on
plot(1:20,relerr02_se(3,:),'*')
hold on
plot(1:20,relerr02_se(4,:),'*')
hold on
plot(1:20,relerr02_se(5,:),'*')


% check the uniqueness of the low rank = 2 with correlation coefficient like above
% the result is unique
Uhat02_r1_se = cpd(N_chb02_se,2);
for j = 1:chan
    Uhat02_r1_se{1}(j,1) = Uhat02_r1_se{1}(j,1)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
    Uhat02_r1_se{1}(j,2) = Uhat02_r1_se{1}(j,2)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
end

Uhat02_r2_se = cpd(N_chb02_se,2);
for j = 1:chan
    Uhat02_r2_se{1}(j,1) = Uhat02_r2_se{1}(j,1)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
    Uhat02_r2_se{1}(j,2) = Uhat02_r2_se{1}(j,2)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
end

Uhat02_r3_se = cpd(N_chb02_se,2);
for j = 1:chan
    Uhat02_r3_se{1}(j,1) = Uhat02_r3_se{1}(j,1)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
    Uhat02_r3_se{1}(j,2) = Uhat02_r3_se{1}(j,2)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
end

Uhat02_r4_se = cpd(N_chb02_se,2);
for j = 1:chan
    Uhat02_r4_se{1}(j,1) = Uhat02_r4_se{1}(j,1)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
    Uhat02_r4_se{1}(j,2) = Uhat02_r4_se{1}(j,2)*std(double(data_chb02.data(j,256*2972+1:256*2977)));
end

for k = 1:3 
    s = size(Uhat02_r1_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r1_se{k}(i,j),Uhat02_r2_se{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 
    s = size(Uhat02_r1_se{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat02_r1_se{k}(i,j),Uhat02_r3_se{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 
    s = size(Uhat02_r1_se{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat02_r1_se{k}(i,j),Uhat02_r4_se{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 
    s = size(Uhat02_r2_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r2_se{k}(i,j),Uhat02_r3_se{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3
    s = size(Uhat02_r2_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r2_se{k}(i,j),Uhat02_r4_se{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 
    s = size(Uhat02_r3_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat02_r3_se{k}(i,j),Uhat02_r4_se{k}(i,j))~= 1
                disp('3&4: solution is not unique');
                break
            end
        end
    end
end   
disp('3&4: solution is unique');


%% plot channel, time, freq with rank 1 & 2 of the two 5s tensor

% plot the channel of tesnor with and withour seizure with rank 1
subplot(2,1,1);
plot(1:23,-Uhat02_r1{1}(:,1),'*g');
hold on
plot(1:23,Uhat02_r1_se{1}(:,1),'*r');
xticks(1:1:23);
legend('without seizure', 'with seizure')
title('Channel with r = 1')
grid on
% plot the channel of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(1:23,Uhat02_r1{1}(:,2),'*g');
hold on
plot(1:23,Uhat02_r1_se{1}(:,2),'*r');
xticks(1:1:23);
legend('without seizure', 'with seizure')
title('Channel with r = 2')
grid on

% plot the time of tesnor with and withour seizure with rank 1
t1 = 1:5*256;
subplot(2,1,1);
plot(t1/256,Uhat02_r1{2}(:,1),'g');
hold on
plot(t1/256,Uhat02_r1_se{2}(:,1),'r');
legend('without seizure', 'with seizure')
title('Time with r = 1')
grid on
% plot the time of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(t1/256,Uhat02_r1{2}(:,2),'g');
hold on
plot(t1/256,Uhat02_r1_se{2}(:,2),'r');
legend('without seizure', 'with seizure')
title('Time with r = 2')
grid on

% plot the frequency of tesnor with and withour seizure with rank 1
subplot(2,1,1);
plot(1:75,-Uhat02_r1{3}(:,1),'g');
hold on
plot(1:75,Uhat02_r1_se{3}(:,1),'r');
legend('without seizure', 'with seizure')
title('Frequency with r = 1')
grid on
% plot the requency of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(1:75,Uhat02_r1{3}(:,2),'g');
hold on
plot(1:75,Uhat02_r1_se{3}(:,2),'r');
legend('without seizure', 'with seizure')
title('Frequency with r = 2')
grid on

%% create 10s tensor from 2967s to 2977s, the seizure starts from 2972s

% first try with one example to obtain the size of time and frequency after CWT like above
mean1_02_com = mean(double(data_chb02.data(1,256*2967+1:256*2977)));
std1_02_com = std(double(data_chb02.data(1,256*2967+1:256*2977)));
z1_02_com = abs(cwt((double(data_chb02.data(1,256*2967+1:256*2977))-mean1_02_com)/std1_02_com));
time_02_com = size(data_chb02.data(1,256*2967+1:256*2977),2);
freq_02_com = size(z1_02_com,1);

% build the 10s tenosr with seizure in last 5s
N_chb02_com = NaN(chan, time_02_com, freq_02_com);
 for i= 1:chan
     mean02_com = mean(double(data_chb02.data(1,256*2967+1:256*2977)));
     std02_com = std(double(data_chb02.data(1,256*2967+1:256*2977)));
     z_02_com = abs(cwt((double(data_chb02.data(i,256*2967+1:256*2977))-mean02_com)/std02_com));
     for j = 1:time_02_com
         for k = 1:freq_02_com
             N_chb02_com(i, j, k) = z_02_com(k, j);
         end
     end
 end
N_chb02_com = fmt(N_chb02_com);

% check the low rank like above, the result is 2
relerr02_com=[];
con02_com =[];
for i = 1:20
    for n = 1:5
        Uhat02_com = cpd(N_chb02_com,i);
        relerr02_com(n,i) = frobcpdres(N_chb02_com, Uhat02_com)/frob(N_chb02_com);
        con02_com(n,i) = corcond(N_chb02_com,Uhat02_com,0);
    end
end

subplot(2,1,1);
plot(1:20,con02_com(1,:),'-');
hold on
plot(1:20,con02_com(2,:),'-');
hold on
plot(1:20,con02_com(3,:),'-');
hold on
plot(1:20,con02_com(4,:),'-');
hold on
plot(1:20,con02_com(5,:),'-')
axis([0 20 -100 100]);

subplot(2,1,2); 
plot(1:20,relerr02_com(1,:),'*')
hold on
plot(1:20,relerr02_com(2,:),'*')
hold on
plot(1:20,relerr02_com(3,:),'*')
hold on
plot(1:20,relerr02_com(4,:),'*')
hold on
plot(1:20,relerr02_com(5,:),'*')

% check the uniqueness, the result is unique
Uhat02_r1_com = cpd(N_chb02_com,2);
for j = 1:chan
    Uhat02_r1_com{1}(j,1) = Uhat02_r1_com{1}(j,1)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
    Uhat02_r1_com{1}(j,2) = Uhat02_r1_com{1}(j,2)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
end

Uhat02_r2_com = cpd(N_chb02_com,2);
for j = 1:chan
    Uhat02_r2_com{1}(j,1) = Uhat02_r2_com{1}(j,1)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
    Uhat02_r2_com{1}(j,2) = Uhat02_r2_com{1}(j,2)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
end

Uhat02_r3_com = cpd(N_chb02_com,2);
for j = 1:chan
    Uhat02_r3_com{1}(j,1) = Uhat02_r3_com{1}(j,1)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
    Uhat02_r3_com{1}(j,2) = Uhat02_r3_com{1}(j,2)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
end

Uhat02_r4_com = cpd(N_chb02_com,2);
for j = 1:chan
    Uhat02_r4_com{1}(j,1) = Uhat02_r4_com{1}(j,1)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
    Uhat02_r4_com{1}(j,2) = Uhat02_r4_com{1}(j,2)*std(double(data_chb02.data(j,256*2967+1:256*2977)));
end

for k = 1:3 
    s_com = size(Uhat02_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat02_r1_com{k}(i,j),Uhat02_r2_com{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 
    s_com = size(Uhat02_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1:s_com(1,1)
            if corrcoef(Uhat02_r1_com{k}(i,j),Uhat02_r3_com{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 
    s_com = size(Uhat02_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1:s_com(1,1)
            if corrcoef(Uhat02_r1_com{k}(i,j),Uhat02_r4_com{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 
    s_com = size(Uhat02_r2_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat02_r2_com{k}(i,j),Uhat02_r3_com{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 
    s_com = size(Uhat02_r2_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat02_r2_com{k}(i,j),Uhat02_r4_com{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 
    s_com = size(Uhat02_r3_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat02_r3_com{k}(i,j),Uhat02_r4_com{k}(i,j))~= 1
                disp('3&4: solution is not unique');
                break
            end
        end
    end
end   
disp('3&4: solution is unique');

%% plot channel, time, freq with rank 1 & 2 of the 10s tensor

% plot the channel of the tensor with rank 1
plot(1:23,Uhat02_r1_com{1}(:,1),'*c');
hold on
% plot the channel of the tensor with rank 2
plot(1:23,Uhat02_r1_com{1}(:,2),'*m');
xticks(1:1:23);
legend('rank = 1', 'rank = 2')
title('Channel of com-tensor')
grid on

% plot the time of the tensor with rank 1
t =1:2560;
plot(t/256,Uhat02_r1_com{2}(:,1),'c');
hold on
% plot the time of the tensor with rank 2
plot(t/256,Uhat02_r1_com{2}(:,2),'m');
legend('rank = 1', 'rank = 2')
title('Time of com-tensor')
grid on

% plot the frequnecy of the tensor with rank 1
plot(1:85,Uhat02_r1_com{3}(:,1),'c');
hold on
% plot the frequency of the tensor with rank 2
plot(1:85,Uhat02_r1_com{3}(:,2),'m');
legend('rank = 1', 'rank = 2')
title('Frequency of com-tensor')
grid on


%% try with the BTD instead of CWT in 10s tensor ( seizure in last 5s)

% L1 = 1, L2 = 2 in time and frequency mode
L  = [1 2 2];
% do the BTD
Uhat02_btd_com = ll1(N_chb02_com, L);

% plot the channel
plot(1:23,Uhat02_btd_com{1}{1}(:,1),'*c');
xticks(1:1:23);
legend('rank = 1')
title('BTD: Channel of com-tensor')
grid on

% plot the time
t =1:2560;
plot(t/256,Uhat02_btd_com{2}{2}(:,1),'c');
hold on
plot(t/256,Uhat02_btd_com{2}{2}(:,2),'m');
legend('rank = 1', 'rank = 2')
title('BTD: Time of com-tensor')
grid on

% plot the frequency
plot(1:85,Uhat02_btd_com{3}{3}(:,1),'c');
legend('rank = 1')
title('BTD: Frequency of com-tensor')
grid on







