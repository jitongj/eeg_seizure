%Seizure Start Time: 362 seconds
%Seizure End Time: 414 seconds

% basic info
data_chb01 = pop_loadset('chb01_03.set');
chan = 23;

%% create 5s tensor without seizure, choose a time period from 100s to 105s

% first try with one example to obtain the size of time and frequency after CWT

% get the mean of first channel
mean1_03_swo = mean(double(data_chb03.data(1,256*100+1:256*105)));
% get the std of first channel
std1_03_swo = std(double(data_chb03.data(1,256*100+1:256*105)));
% first standardize, then do wavelet transformation 
z1_03 = abs(cwt((double(data_chb03.data(1,256*100+1:256*105))-mean1_03_swo)/std1_03_swo));
% get the size of time
time_03 = size(data_chb03.data(1,256*100+1:256*105),2);
% get the size of frequency
freq_03 = size(z1_03,1);

% build the seperated 5s tensor without seizure
N_chb03 = NaN(chan, time_03, freq_03);
 for i= 1:chan
     % get the mean of each channel
     mean03_swo =  mean(double(data_chb03.data(i,256*100+1:256*105)));
     % get the std of each channel
     std03_swo = std(double(data_chb03.data(i,256*100+1:256*105)));
     % standardize, then do wavelet transformation
     z_03 = abs(cwt((double(data_chb03.data(i,256*100+1:256*105))-mean03_swo)/std03_swo));
     % build tensor
     for j = 1:time_03
         for k = 1:freq_03
             N_chb03(i, j, k) = z_03(k, j);
         end
     end
 end
N_chb03 = fmt(N_chb03);

% check low rank
relerr03=[];
con03 =[];
% randomly initialized CPD models with different ranks ranging from 1 to 20 for 5 times
for i = 1:20
    for n = 1:5
        Uhat03 = cpd(N_chb03,i);
        % calculate the relative error
        relerr03(n,i) = frobcpdres(N_chb03, Uhat03)/frob(N_chb03);
        % calculate the core consistency
        con03(n,i) = corcond(N_chb03,Uhat03,0);
    end
end

% plot the core consistency to find the low rank
subplot(2,1,1);
plot(1:20,con03(1,:),'-');
hold on
plot(1:20,con03(2,:),'-');
hold on
plot(1:20,con03(3,:),'-');
hold on
plot(1:20,con03(4,:),'-');
hold on
plot(1:20,con03(5,:),'-')
axis([1,20,-100,100]);

% plot the relative error to find the low rank
subplot(2,1,2); 
plot(1:20,relerr03(1,:),'*')
hold on
plot(1:20,relerr03(2,:),'*')
hold on
plot(1:20,relerr03(3,:),'*')
hold on
plot(1:20,relerr03(4,:),'*')
hold on
plot(1:20,relerr03(5,:),'*')

% result: low rank = 2


% ran CPD 4 times to check the uniqueness of the low rank = 2

% ran randomly the 1st CPD
Uhat03_r1 = cpd(N_chb03,2);
%after decomposition, multiply with the standard deviation to preserve topographic information
for j = 1:chan
    Uhat03_r1{1}(j,1) = Uhat03_r1{1}(j,1)*std(double(data_chb03.data(i,256*100+1:256*105)));
    Uhat03_r1{1}(j,2) = Uhat03_r1{1}(j,2)*std(double(data_chb03.data(i,256*100+1:256*105)));
end

% ran randomly the 2nd CPD
Uhat03_r2 = cpd(N_chb03,2);
for j = 1:chan
    Uhat03_r2{1}(j,1) = Uhat03_r2{1}(j,1)*std(double(data_chb03.data(i,256*100+1:256*105)));
    Uhat03_r2{1}(j,2) = Uhat03_r2{1}(j,2)*std(double(data_chb03.data(i,256*100+1:256*105)));
end

% ran randomly the 3rd CPD
Uhat03_r3 = cpd(N_chb03,2);
for j = 1:chan
    Uhat03_r3{1}(j,1) = Uhat03_r3{1}(j,1)*std(double(data_chb03.data(i,256*100+1:256*105)));
    Uhat03_r3{1}(j,2) = Uhat03_r3{1}(j,2)*std(double(data_chb03.data(i,256*100+1:256*105)));
end

% ran randomly the 4th CPD
Uhat03_r4 = cpd(N_chb03,2);
for j = 1:chan
    Uhat03_r4{1}(j,1) = Uhat03_r4{1}(j,1)*std(double(data_chb03.data(i,256*100+1:256*105)));
    Uhat03_r4{1}(j,2) = Uhat03_r4{1}(j,2)*std(double(data_chb03.data(i,256*100+1:256*105)));
end

% compare pairwise of the component using Pearson's correlation coefficient
for k = 1:3 %3 mode
    s = size(Uhat03_r1{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r1{k}(i,j),Uhat03_r2{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r1{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat03_r1{k}(i,j),Uhat03_r3{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r1{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat03_r1{k}(i,j),Uhat03_r4{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r2{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r2{k}(i,j),Uhat03_r3{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r2{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r2{k}(i,j),Uhat03_r4{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r3{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r3{k}(i,j),Uhat03_r4{k}(i,j))~= 1
                disp('3&4: solution is not unique');
                break
            end
        end
    end
end   
disp('3&4: solution is unique');

% result: unique


%% create 5s tensor with seizure, choose a time period from 362s to 367s

% first try with one example to obtain the size of time and frequency after CWT like above
mean1_03_sw = mean(double(data_chb03.data(1,256*362+1:256*367)));
std1_03_sw = std(double(data_chb03.data(1,256*362+1:256*367)));
z1_03_se = abs(cwt((double(data_chb03.data(1,256*362+1:256*367))-mean1_03_sw)/std1_03_sw));
time_03_se = size(data_chb03.data(1,256*362+1:256*367),2);
freq_03_se = size(z1_03_se,1);

% build the seperated 5s tensor with seizure like above
N_chb03_se = NaN(chan, time_03_se, freq_03_se);
 for i= 1:chan
     mean03_sw = mean(double(data_chb03.data(i,256*362+1:256*367)));
     std03_sw = std(double(data_chb03.data(i,256*362+1:256*367)));
     z_03_sw = abs(cwt((double(data_chb03.data(i,256*362+1:256*367))-mean03_sw)/std03_sw));
     for j = 1:time_03_se
         for k = 1:freq_03_se
             N_chb03_se(i, j, k) = z_03_sw(k, j);
         end
     end
 end
N_chb03_se = fmt(N_chb03_se);



% check the low rank with core consistency and relative error like above
% the result is 2
relerr03_se=[];
con03_se =[];
for i = 1:20
    for n = 1:5
        Uhat03_se = cpd(N_chb03_se,i);
        relerr03_se(n,i) = frobcpdres(N_chb03_se, Uhat03_se)/frob(N_chb03_se);
        con03_se(n,i) = corcond(N_chb03_se,Uhat03_se,0);
    end
end

subplot(2,1,1);
plot(1:20,con03_se(1,:),'-');
hold on
plot(1:20,con03_se(2,:),'-');
hold on
plot(1:20,con03_se(3,:),'-');
hold on
plot(1:20,con03_se(4,:),'-');
hold on
plot(1:20,con03_se(5,:),'-')
axis([1,20,-100,100]);


subplot(2,1,2); 
plot(1:20,relerr03_se(1,:),'*')
hold on
plot(1:20,relerr03_se(2,:),'*')
hold on
plot(1:20,relerr03_se(3,:),'*')
hold on
plot(1:20,relerr03_se(4,:),'*')
hold on
plot(1:20,relerr03_se(5,:),'*')

%result: low rank = 2

% check the uniqueness of the low rank = 2 with correlation coefficient like above
% the result is unique
Uhat03_r1_se = cpd(N_chb03_se,2);
for j = 1:chan
    Uhat03_r1_se{1}(j,1) = Uhat03_r1_se{1}(j,1)*std(double(data_chb03.data(j,256*362+1:256*367)));
    Uhat03_r1_se{1}(j,2) = Uhat03_r1_se{1}(j,2)*std(double(data_chb03.data(j,256*362+1:256*367)));
end

Uhat03_r2_se = cpd(N_chb03_se,2);
for j = 1:chan
    Uhat03_r2_se{1}(j,1) = Uhat03_r2_se{1}(j,1)*std(double(data_chb03.data(j,256*362+1:256*367)));
    Uhat03_r2_se{1}(j,2) = Uhat03_r2_se{1}(j,2)*std(double(data_chb03.data(j,256*362+1:256*367)));
end

Uhat03_r3_se = cpd(N_chb03_se,2);
for j = 1:chan
    Uhat03_r3_se{1}(j,1) = Uhat03_r3_se{1}(j,1)*std(double(data_chb03.data(j,256*362+1:256*367)));
    Uhat03_r3_se{1}(j,2) = Uhat03_r3_se{1}(j,2)*std(double(data_chb03.data(j,256*362+1:256*367)));
end

Uhat03_r4_se = cpd(N_chb03_se,2);
for j = 1:chan
    Uhat03_r4_se{1}(j,1) = Uhat03_r4_se{1}(j,1)*std(double(data_chb03.data(j,256*362+1:256*367)));
    Uhat03_r4_se{1}(j,2) = Uhat03_r4_se{1}(j,2)*std(double(data_chb03.data(j,256*362+1:256*367)));
end


for k = 1:3 %3 mode
    s = size(Uhat03_r1_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r1_se{k}(i,j),Uhat03_r2_se{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r1_se{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat03_r1_se{k}(i,j),Uhat03_r3_se{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r1_se{k});
    for j = 1:s(1,2)
        for i = 1:s(1,1)
            if corrcoef(Uhat03_r1_se{k}(i,j),Uhat03_r4_se{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r2_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r2_se{k}(i,j),Uhat03_r3_se{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r2_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r2_se{k}(i,j),Uhat03_r4_se{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 %3 mode
    s = size(Uhat03_r3_se{k});
    for j = 1:s(1,2)
        for i = 1: s(1,1)
            if corrcoef(Uhat03_r3_se{k}(i,j),Uhat03_r4_se{k}(i,j))~= 1
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
plot(1:23,Uhat03_r1{1}(:,1),'*g');
hold on
plot(1:23,-Uhat03_r1_se{1}(:,1),'*r');
xticks(1:1:23);
legend('without seizure', 'with seizure')
title('Channel with r = 1')
grid on
% plot the channel of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(1:23,Uhat03_r1{1}(:,2),'*g');
hold on
plot(1:23,Uhat03_r1_se{1}(:,2),'*r');
xticks(1:1:23);
legend('without seizure', 'with seizure')
title('Channel with r = 2')
grid on


% plot the time of tesnor with and withour seizure with rank 1
t1 = 1:5*256;
subplot(2,1,1);
plot(t1/256,Uhat03_r1{2}(:,1),'g');
hold on
plot(t1/256,Uhat03_r1_se{2}(:,1),'r');
legend('without seizure', 'with seizure')
title('Time with r = 1')
grid on
% plot the time of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(t1/256,Uhat03_r1{2}(:,2),'g');
hold on
plot(t1/256,Uhat03_r1_se{2}(:,2),'r');
legend('without seizure', 'with seizure')
title('Time with r = 2')
grid on

% plot the frequency of tesnor with and withour seizure with rank 1
subplot(2,1,1);
plot(1:75,Uhat03_r1{3}(:,1),'g');
hold on
plot(1:75,-Uhat03_r1_se{3}(:,1),'r');
legend('without seizure', 'with seizure')
title('Frequency with r = 1')
grid on
% plot the requency of tesnor with and withour seizure with rank 2
subplot(2,1,2);
plot(1:75,Uhat03_r1{3}(:,2),'g');
hold on
plot(1:75,Uhat03_r1_se{3}(:,2),'r');
legend('without seizure', 'with seizure')
title('Frequency with r = 2')
grid on


%% create 10s tensor from 357s to 367s, the seizure starts from 362s

% first try with one example to obtain the size of time and frequency after CWT like above
mean1_03_com = mean(double(data_chb03.data(1,256*357+1:256*367)));
std1_03_com = std(double(data_chb03.data(1,256*357+1:256*367)));
z1_03_com = abs(cwt((double(data_chb03.data(1,256*357+1:256*367))-mean1_03_com)/std1_03_com));
time_03_com = size(data_chb03.data(1,256*357+1:256*367),2);
freq_03_com = size(z1_03_com,1);

% build the 10s tenosr with seizure in last 5s
N_chb03_com = NaN(chan, time_03_com, freq_03_com);
 for i= 1:chan
     mean03_com = mean(double(data_chb03.data(1,256*357+1:256*367)));
     std03_com = std(double(data_chb03.data(1,256*357+1:256*367)));
     z_03_com = abs(cwt((double(data_chb03.data(i,256*357+1:256*367))-mean03_com)/std03_com));
     for j = 1:time_03_com
         for k = 1:freq_03_com
             N_chb03_com(i, j, k) = z_03_com(k, j);
         end
     end
 end
N_chb03_com = fmt(N_chb03_com);


% check the low rank like above, the result is 3
relerr03_com=[];
con03_com =[];
for i = 1:20
    for n = 1:5
        Uhat03_com = cpd(N_chb03_com,i);
        relerr03_com(n,i) = frobcpdres(N_chb03_com, Uhat03_com)/frob(N_chb03_com);
        con03_com(n,i) = corcond(N_chb03_com,Uhat03_com,0);
    end
end

subplot(2,1,1);
plot(1:20,con03_com(1,:),'-');
hold on
plot(1:20,con03_com(2,:),'-');
hold on
plot(1:20,con03_com(3,:),'-');
hold on
plot(1:20,con03_com(4,:),'-');
hold on
plot(1:20,con03_com(5,:),'-')
axis([0 20 -100 100]);

subplot(2,1,2); 
plot(1:20,relerr03_com(1,:),'*')
hold on
plot(1:20,relerr03_com(2,:),'*')
hold on
plot(1:20,relerr03_com(3,:),'*')
hold on
plot(1:20,relerr03_com(4,:),'*')
hold on
plot(1:20,relerr03_com(5,:),'*')

% check the uniqueness, the result is unique
Uhat03_r1_com = cpd(N_chb03_com,3);
for j = 1:chan
    Uhat03_r1_com{1}(j,1) = Uhat03_r1_com{1}(j,1)*std(double(data_chb03.data(j,256*357+1:256*367)));
    Uhat03_r1_com{1}(j,2) = Uhat03_r1_com{1}(j,2)*std(double(data_chb03.data(j,256*357+1:256*367)));
end

Uhat03_r2_com = cpd(N_chb03_com,3);
for j = 1:chan
    Uhat03_r2_com{1}(j,1) = Uhat03_r2_com{1}(j,1)*std(double(data_chb03.data(j,256*357+1:256*367)));
    Uhat03_r2_com{1}(j,2) = Uhat03_r2_com{1}(j,2)*std(double(data_chb03.data(j,256*357+1:256*367)));
end

Uhat03_r3_com = cpd(N_chb03_com,3);
for j = 1:chan
    Uhat03_r3_com{1}(j,1) = Uhat03_r3_com{1}(j,1)*std(double(data_chb03.data(j,256*357+1:256*367)));
    Uhat03_r3_com{1}(j,2) = Uhat03_r3_com{1}(j,2)*std(double(data_chb03.data(j,256*357+1:256*367)));
end

Uhat03_r4_com = cpd(N_chb03_com,3);
for j = 1:chan
    Uhat03_r4_com{1}(j,1) = Uhat03_r4_com{1}(j,1)*std(double(data_chb03.data(j,256*357+1:256*367)));
    Uhat03_r4_com{1}(j,2) = Uhat03_r4_com{1}(j,2)*std(double(data_chb03.data(j,256*357+1:256*367)));
end

for k = 1:3 
    s_com = size(Uhat03_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat03_r1_com{k}(i,j),Uhat03_r2_com{k}(i,j))~= 1
                disp('1&2: solution is not unique');
                break
            end
        end
    end
end
disp('1&2: solution is unique');

for k = 1:3 
    s_com = size(Uhat03_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1:s_com(1,1)
            if corrcoef(Uhat03_r1_com{k}(i,j),Uhat03_r3_com{k}(i,j))~= 1
                disp('1&3: solution is not unique');
                break
            end
        end
    end
end
disp('1&3: solution is unique');

for k = 1:3 
    s_com = size(Uhat03_r1_com{k});
    for j = 1:s_com(1,2)
        for i = 1:s_com(1,1)
            if corrcoef(Uhat03_r1_com{k}(i,j),Uhat03_r4_com{k}(i,j))~= 1
                disp('1&4: solution is not unique');
                break
            end
        end
    end
end
disp('1&4: solution is unique');

for k = 1:3 
    s_com = size(Uhat03_r2_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat03_r2_com{k}(i,j),Uhat03_r3_com{k}(i,j))~= 1
                disp('2&3: solution is not unique');
                break
            end
        end
    end
end   
disp('2&3: solution is unique');

for k = 1:3 
    s_com = size(Uhat03_r2_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat03_r2_com{k}(i,j),Uhat03_r4_com{k}(i,j))~= 1
                disp('2&4: solution is not unique');
                break
            end
        end
    end
end             
disp('2&4: solution is unique');

for k = 1:3 
    s_com = size(Uhat03_r3_com{k});
    for j = 1:s_com(1,2)
        for i = 1: s_com(1,1)
            if corrcoef(Uhat03_r3_com{k}(i,j),Uhat03_r4_com{k}(i,j))~= 1
                disp('3&4: solution is not unique');
                break
            end
        end
    end
end   
disp('3&4: solution is unique');

%% plot channel, time, freq with rank 1 & 2 of the 10s tensor

% plot the channel of the tensor with rank 1
plot(1:23,-Uhat03_r1_com{1}(:,1),'*c');
hold on
% plot the channel of the tensor with rank 2
plot(1:23,Uhat03_r1_com{1}(:,2),'*m');
hold on
% plot the channel of the tensor with rank 3
plot(1:23,Uhat03_r1_com{1}(:,3),'*b');
xticks(1:1:23);
legend('rank = 1', 'rank = 2','rank = 3')
title('Channel of com-tensor')
grid on

% plot the time of the tensor with rank 1
t =1:2560;
plot(t/256,Uhat03_r1_com{2}(:,1),'c');
hold on
% plot the time of the tensor with rank 2
plot(t/256,Uhat03_r1_com{2}(:,2),'m');
hold on
% plot the time of the tensor with rank 3
plot(t/256,Uhat03_r1_com{2}(:,3),'b');
legend('rank = 1', 'rank = 2','rank = 3')
title('Time of com-tensor')
grid on

% plot the frequnecy of the tensor with rank 1
plot(1:85,-Uhat03_r1_com{3}(:,1),'c');
hold on
% plot the frequency of the tensor with rank 2
plot(1:85,Uhat03_r1_com{3}(:,2),'m');
hold on
% plot the frequency of the tensor with rank 3
plot(1:85,Uhat03_r1_com{3}(:,3),'b');
legend('rank = 1', 'rank = 2','rank = 3')
title('Frequency of com-tensor')
grid on


%% try with the BTD instead of CWT in 10s tensor ( seizure in last 5s)

% L1 = 1, L2 = 2 in time and frequency mode
L  = [1 2 2];
% do the BTD
Uhat03_btd_com = ll1(N_chb03_com, L);

% plot the channel
plot(1:23,Uhat03_btd_com{1}{1}(:,1),'*c');
xticks(1:1:23);
legend('rank = 1')
title('BTD: Channel of com-tensor')
grid on

% plot the time
t =1:2560;
plot(t/256,Uhat03_btd_com{2}{2}(:,1),'c');
hold on
plot(t/256,Uhat03_btd_com{2}{2}(:,2),'m');
legend('rank = 1', 'rank = 2')
title('BTD: Time of com-tensor')
grid on

% plot the frequency
plot(1:85,Uhat03_btd_com{3}{3}(:,1),'c');
legend('rank = 1')
title('BTD: Frequency of com-tensor')
grid on





