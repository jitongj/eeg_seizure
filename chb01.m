% Seizure Start Time: 2996 seconds
% Seizure End Time: 3036 seconds

% basic info
data_chb01 = pop_loadset('chb01_03.set');
chan = 23;

%% create 5s tensor without seizure, choose a time period from 500s to 505s

% first try with one example to obtain the size of time and frequency after CWT

% get the mean of first channel
mean1_01_swo = mean(double(data_chb01.data(1,256*500+1:256*505)));
% get the std of first channel
std1_01_swo = std(double(data_chb01.data(1,256*500+1:256*505))); 
% first standardize, then do wavelet transformation 
z1_01 = abs(cwt((double(data_chb01.data(1,256*500+1:256*505))-mean1_01_swo)/std1_01_swo));
% get the size of time
time_01 = size(data_chb01.data(1,256*500+1:256*505),2);
% get the size of frequency
freq_01 = size(z1_01,1);

% build the seperated 5s tensor without seizure
N_chb01 = NaN(chan, time_01, freq_01);
 for i= 1:chan
     % get the mean of each channel
     mean_swo = mean(double(data_chb01.data(i,256*500+1:256*505)));
     % get the std of each channel
     std_swo = std(double(data_chb01.data(i,256*500+1:256*505)));
     % standardize, then do wavelet transformation
     z_01= abs(cwt((double(data_chb01.data(i,256*500+1:256*505))-mean_swo)/std_swo)); 
     % build tensor
     for j = 1:time_01
         for k = 1:freq_01
             N_chb01(i, j, k) = z_01(k, j);
         end
     end
 end
N_chb01 = fmt(N_chb01);




%check the low rank
relerr01=[];
con01 =[];

% randomly initialized CPD models with different ranks ranging from 1 to 20 for 5 times
for i = 1:20
    for n = 1:5
        Uhat01 = cpd(N_chb01,i);
        % calculate the relative error
        relerr01(n,i) = frobcpdres(N_chb01, Uhat01)/frob(N_chb01);
        % calculate the core consistency
        con01(n,i) = corcond(N_chb01,Uhat01,0);
    end
end

% plot the core consistency to find the low rank
subplot(2,1,1);
plot(1:20,con01(1,:),'-');
hold on
plot(1:20,con01(2,:),'-');
hold on
plot(1:20,con01(3,:),'-');
hold on
plot(1:20,con01(4,:),'-');
hold on
plot(1:20,con01(5,:),'-')

% plot the relative error to find the low rank
subplot(2,1,2); 
plot(1:20,relerr01(1,:),'*')
hold on
plot(1:20,relerr01(2,:),'*')
hold on
plot(1:20,relerr01(3,:),'*')
hold on
plot(1:20,relerr01(4,:),'*')
hold on
plot(1:20,relerr01(5,:),'*')

% result: low rank = 2


% ran CPD 4 times to check the uniqueness of the low rank = 2

% ran randomly the 1st CPD
Uhat01_r1 = cpd(N_chb01,2); 
%after decomposition, multiply with the standard deviation to preserve topographic information
for j = 1:chan
    Uhat01_r1{1}(j,1) = Uhat01_r1{1}(j,1)*std(double(data_chb01.data(j,256*500+1:256*505)));
    Uhat01_r1{1}(j,2) = Uhat01_r1{1}(j,2)*std(double(data_chb01.data(j,256*500+1:256*505)));
end

% ran randomly the 2nd CPD
Uhat01_r2 = cpd(N_chb01,2);
for j = 1:chan
    Uhat01_r2{1}(j,1) = Uhat01_r2{1}(j,1)*std(double(data_chb01.data(j,256*500+1:256*505)));
    Uhat01_r2{1}(j,2) = Uhat01_r2{1}(j,2)*std(double(data_chb01.data(j,256*500+1:256*505)));
end

% ran randomly the 3rd CPD
Uhat01_r3 = cpd(N_chb01,2);
for j = 1:chan
    Uhat01_r3{1}(j,1) = Uhat01_r3{1}(j,1)*std(double(data_chb01.data(j,256*500+1:256*505)));
    Uhat01_r3{1}(j,2) = Uhat01_r3{1}(j,2)*std(double(data_chb01.data(j,256*500+1:256*505)));
end

% ran randomly the 4th CPD
Uhat01_r4 = cpd(N_chb01,2);
for j = 1:chan
    Uhat01_r4{1}(j,1) = Uhat01_r4{1}(j,1)*std(double(data_chb01.data(j,256*500+1:256*505)));
    Uhat01_r4{1}(j,2) = Uhat01_r4{1}(j,2)*std(double(data_chb01.data(j,256*500+1:256*505)));
end

% compare pairwise of the component using Pearson's correlation coefficient
for k = 1:3 %3 mode
    s = size(Uhat01_r1{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r1{k}(i,j),Uhat01_r2{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat01_r1{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat01_r1{k}(i,j),Uhat01_r3{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat01_r1{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat01_r1{k}(i,j),Uhat01_r4{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat01_r2{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r2{k}(i,j),Uhat01_r3{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat01_r2{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r2{k}(i,j),Uhat01_r4{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat01_r3{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r3{k}(i,j),Uhat01_r4{k}(i,j))~= 1
                disp('3&4: solution is not unique');
                break
            end
        end
    end
end   
disp('3&4: solution is unique');

% result: unique


%% create 5s tensor with seizure, choose a time period from 2996s to 3001s

% first try with one example to obtain the size of time and frequency after CWT like above
mean1_01_sw = mean(double(data_chb01.data(1,256*2996+1:256*3001))); 
std1_01_sw = std(double(data_chb01.data(1,256*2996+1:256*3001)));
z1_01_se = abs(cwt((double(data_chb01.data(1,256*2996+1:256*3001))-mean1_01_sw)/std1_01_sw));
time_01_se = size(data_chb01.data(1,256*2996+1:256*3001),2);
freq_01_se = size(z1_01_se,1);

% build the seperated 5s tensor with seizure like above
N_chb01_se = NaN(chan, time_01_se, freq_01_se);
 for i= 1:chan
     mean01_sw = mean(double(data_chb01.data(i,256*2996+1:256*3001)));
     std01_sw = std(double(data_chb01.data(i,256*2996+1:256*3001)));
     z_01_sw = abs(cwt((double(data_chb01.data(i,256*2996+1:256*3001))-mean01_sw)/std01_sw));
     for j = 1:time_01_se
         for k = 1:freq_01_se
             N_chb01_se(i, j, k) = z_01_sw(k, j);
         end
     end
 end
N_chb01_se = fmt(N_chb01_se);


% check the low rank with core consistency and relative error like above
% the result is 2

relerr01_se=[];
con01_se =[];
for i = 1:20
    for n = 1:5
        Uhat01_se = cpd(N_chb01_se,i);
        relerr01_se(n,i) = frobcpdres(N_chb01_se, Uhat01_se)/frob(N_chb01_se);
        con01_se(n,i) = corcond(N_chb01_se,Uhat01_se,0);
    end
end

subplot(2,1,1);
plot(1:20,con01_se(1,:),'-');
hold on
plot(1:20,con01_se(2,:),'-');
hold on
plot(1:20,con01_se(3,:),'-');
hold on
plot(1:20,con01_se(4,:),'-');
hold on
plot(1:20,con01_se(5,:),'-')

subplot(2,1,2); 
plot(1:20,relerr01_se(1,:),'*')
hold on
plot(1:20,relerr01_se(2,:),'*')
hold on
plot(1:20,relerr01_se(3,:),'*')
hold on
plot(1:20,relerr01_se(4,:),'*')
hold on
plot(1:20,relerr01_se(5,:),'*')

%result: low rank = 2

% check the uniqueness of the low rank = 2 with correlation coefficient like above
% the result is unique

Uhat01_r1_se = cpd(N_chb01_se,2);
for j = 1:chan
    Uhat01_r1_se{1}(j,1) = Uhat01_r1_se{1}(j,1)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
    Uhat01_r1_se{1}(j,2) = Uhat01_r1_se{1}(j,2)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
end

Uhat01_r2_se = cpd(N_chb01_se,2);
for j = 1:chan
    Uhat01_r2_se{1}(j,1) = Uhat01_r2_se{1}(j,1)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
    Uhat01_r2_se{1}(j,2) = Uhat01_r2_se{1}(j,2)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
end

Uhat01_r3_se = cpd(N_chb01_se,2);
for j = 1:chan
    Uhat01_r3_se{1}(j,1) = Uhat01_r3_se{1}(j,1)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
    Uhat01_r3_se{1}(j,2) = Uhat01_r3_se{1}(j,2)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
end

Uhat01_r4_se = cpd(N_chb01_se,2);
for j = 1:chan
    Uhat01_r4_se{1}(j,1) = Uhat01_r4_se{1}(j,1)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
    Uhat01_r4_se{1}(j,2) = Uhat01_r4_se{1}(j,2)*std(double(data_chb01.data(j,256*2996+1:256*3001)));
end

for k = 1:3 
    s = size(Uhat01_r1_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r1_se{k}(i,j),Uhat01_r2_se{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 
    s = size(Uhat01_r1_se{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat01_r1_se{k}(i,j),Uhat01_r3_se{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 
    s = size(Uhat01_r1_se{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat01_r1_se{k}(i,j),Uhat01_r4_se{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 
    s = size(Uhat01_r2_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r2_se{k}(i,j),Uhat01_r3_se{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 
    s = size(Uhat01_r2_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r2_se{k}(i,j),Uhat01_r4_se{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 
    s = size(Uhat01_r3_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat01_r3_se{k}(i,j),Uhat01_r4_se{k}(i,j))~= 1
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
plot(1:23,Uhat01_r1{1}(:,1),'*g');
hold on
plot(1:23,-Uhat01_r1_se{1}(:,1),'*r');
xticks(1:1:23);
legend('without seizure', 'with seizure')
title('Channel with r = 1')
grid on
% plot the channel of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(1:23,-Uhat01_r1{1}(:,2),'*g');
hold on
plot(1:23,-Uhat01_r1_se{1}(:,2),'*r');
xticks(1:1:23);
legend('without seizure', 'with seizure')
title('Channel with r = 2')
grid on

% plot the time of tesnor with and withour seizure with rank 1
t1 = 1:5*256;
subplot(2,1,1);
plot(t1/256,Uhat01_r1{2}(:,1),'g');
hold on
plot(t1/256,Uhat01_r1_se{2}(:,1),'r');
legend('without seizure', 'with seizure')
title('Time with r = 1')
grid on
% plot the time of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(t1/256,Uhat01_r1{2}(:,2),'g');
hold on
plot(t1/256,Uhat01_r1_se{2}(:,2),'r');
legend('without seizure', 'with seizure')
title('Time with r = 2')
grid on

% plot the frequency of tesnor with and withour seizure with rank 1
subplot(2,1,1);
plot(1:75,Uhat01_r1{3}(:,1),'g');
hold on
plot(1:75,-Uhat01_r1_se{3}(:,1),'r');
legend('without seizure', 'with seizure')
title('Frequency with r = 1')
grid on
% plot the requency of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(1:75,-Uhat01_r1{3}(:,2),'g');
hold on
plot(1:75,-Uhat01_r1_se{3}(:,2),'r');
legend('without seizure', 'with seizure')
title('Frequency with r = 2')
grid on



%% create 10s tensor from 2991s to 3001s, the seizure starts from 2996s

% first try with one example to obtain the size of time and frequency after CWT like above
mean1_01_com = mean(double(data_chb01.data(1,256*2991+1:256*3001)));
std1_01_com = std(double(data_chb01.data(1,256*2991+1:256*3001)));
z1_01_com = abs(cwt((double(data_chb01.data(1,256*2991+1:256*3001))-mean1_01_com)/std1_01_com));
time_01_com = size(data_chb01.data(1,256*2991+1:256*3001),2);
freq_01_com = size(z1_01_com,1);

% build the 10s tenosr with seizure in last 5s
N_chb01_com = NaN(chan, time_01_com, freq_01_com);
 for i= 1:chan
     mean_com = mean(double(data_chb01.data(1,256*2991+1:256*3001)));
     std_com = std(double(data_chb01.data(1,256*2991+1:256*3001)));
     z_01_com = abs(cwt((double(data_chb01.data(i,256*2991+1:256*3001))-mean_com)/std_com));
     for j = 1:time_01_com
         for k = 1:freq_01_com
             N_chb01_com(i, j, k) = z_01_com(k, j);
         end
     end
 end
N_chb01_com = fmt(N_chb01_com);

% check the low rank like above, the result is 2
relerr_com=[];
con_com =[];
for i = 1:20
    for n = 1:5
        Uhat_com = cpd(N_chb01_com,i);
        relerr_com(n,i) = frobcpdres(N_chb01_com, Uhat_com)/frob(N_chb01_com);
        con_com(n,i) = corcond(N_chb01_com,Uhat_com,0);
    end
end

subplot(2,1,1);
plot(1:20,con_com(1,:),'-');
hold on
plot(1:20,con_com(2,:),'-');
hold on
plot(1:20,con_com(3,:),'-');
hold on
plot(1:20,con_com(4,:),'-');
hold on
plot(1:20,con_com(5,:),'-')
axis([0 20 -100 100]);

subplot(2,1,2); 
plot(1:20,relerr_com(1,:),'*')
hold on
plot(1:20,relerr_com(2,:),'*')
hold on
plot(1:20,relerr_com(3,:),'*')
hold on
plot(1:20,relerr_com(4,:),'*')
hold on
plot(1:20,relerr_com(5,:),'*')


% check the uniqueness, the result is unique
Uhat01_r1_com = cpd(N_chb01_com,2);
for j = 1:chan
    Uhat01_r1_com{1}(j,1) = Uhat01_r1_com{1}(j,1)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
    Uhat01_r1_com{1}(j,2) = Uhat01_r1_com{1}(j,2)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
end


Uhat01_r2_com = cpd(N_chb01_com,2);
for j = 1:chan
    Uhat01_r2_com{1}(j,1) = Uhat01_r2_com{1}(j,1)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
    Uhat01_r2_com{1}(j,2) = Uhat01_r2_com{1}(j,2)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
end

Uhat01_r3_com = cpd(N_chb01_com,2);
for j = 1:chan
    Uhat01_r3_com{1}(j,1) = Uhat01_r3_com{1}(j,1)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
    Uhat01_r3_com{1}(j,2) = Uhat01_r3_com{1}(j,2)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
end

Uhat01_r4_com = cpd(N_chb01_com,2);
for j = 1:chan
    Uhat01_r1_com{4}(j,1) = Uhat01_r4_com{1}(j,1)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
    Uhat01_r1_com{4}(j,2) = Uhat01_r4_com{1}(j,2)*std(double(data_chb01.data(j,256*2991+1:256*3001)));
end

for k = 1:3 
    s_com = size(Uhat01_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat01_r1_com{k}(i,j),Uhat01_r2_com{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 
    s_com = size(Uhat01_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1:s_com(1,1)
            if corrcoef(Uhat01_r1_com{k}(i,j),Uhat01_r3_com{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 
    s_com = size(Uhat01_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1:s_com(1,1)
            if corrcoef(Uhat01_r1_com{k}(i,j),Uhat01_r4_com{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 
    s_com = size(Uhat01_r2_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat01_r2_com{k}(i,j),Uhat01_r3_com{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 
    s_com = size(Uhat01_r2_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat01_r2_com{k}(i,j),Uhat01_r4_com{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 %
    s_com = size(Uhat01_r3_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat01_r3_com{k}(i,j),Uhat01_r4_com{k}(i,j))~= 1
                disp('3&4: solution is not unique');
                break
            end
        end
    end
end   
disp('3&4: solution is unique');

%% plot channel, time, freq with rank 1 & 2 of the 10s tensor

% plot the channel of the tensor with rank 1
plot(1:23,Uhat01_r1_com{1}(:,1),'*c');
hold on
% plot the channel of the tensor with rank 2
plot(1:23,Uhat01_r1_com{1}(:,2),'*m');
xticks(1:1:23);
legend('rank = 1', 'rank = 2')
title('Channel of com-tensor')
grid on

% plot the time of the tensor with rank 1
t =1:2560;
plot(t/256,Uhat01_r1_com{2}(:,1),'c');
hold on
% plot the time of the tensor with rank 2
plot(t/256,Uhat01_r1_com{2}(:,2),'m');
legend('rank = 1', 'rank = 2')
title('Time of com-tensor')
grid on

% plot the frequnecy of the tensor with rank 1
plot(1:85,Uhat01_r1_com{3}(:,1),'c');
hold on
% plot the frequency of the tensor with rank 2
plot(1:85,Uhat01_r1_com{3}(:,2),'m');
legend('rank = 1', 'rank = 2')
title('Frequency of com-tensor')
grid on


%% try with the BTD instead of CWT in 10s tensor ( seizure in last 5s)

% L1 = 1, L2 = 2 in time and frequency mode
L  = [1 2 2];
% do the BTD
Uhat01_btd_com = ll1(N_chb01_com, L);

% plot the channel
plot(1:23,Uhat01_btd_com{1}{1}(:,1),'*c');
xticks(1:1:23);
legend('rank = 1')
title('BTD: Channel of com-tensor')
grid on

% plot the time
t =1:2560;
plot(t/256,Uhat01_btd_com{2}{2}(:,1),'c');
hold on
plot(t/256,Uhat01_btd_com{2}{2}(:,2),'m');
legend('rank = 1', 'rank = 2')
title('BTD: Time of com-tensor')
grid on

% plot the frequency
plot(1:85,Uhat01_btd_com{3}{3}(:,1),'c');
legend('rank = 1')
title('BTD: Frequency of com-tensor')
grid on


