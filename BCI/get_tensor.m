function [left_train right_train left_test right_test] = get_tensor(subject_name,wname)
    load (subject_name)
    chan = 62;
    time = 1000;
    
    trail_s1d1 = size(session1_data1,3);
    a_s1d1 = session1_data1(:,:,1);  
    mean1_s1d1 = mean(a_s1d1(1,:));
    std1_s1d1 = std(a_s1d1(1,:));
    z1_s1d1 = abs(cwt((a_s1d1(1,:)-mean1_s1d1)/std1_s1d1,wname)); 
    freq_s1d1 = size(z1_s1d1,1);
    tensor_s1d1 = NaN(chan, time, 23);
    Tensor_s1d1 = {};
    for i = 1:trail_s1d1
        b = session1_data1(:,:,i);
        for n= 1:chan
            mean_s1d1 = mean(b(n,:));
            std_s1d1 = std(b(n,:));
            if std_s1d1 == 0
                z_s1d1= abs(cwt(b(n,:)-mean_s1d1,wname));
            else
                z_s1d1= abs(cwt((b(n,:)-mean_s1d1)/std_s1d1,wname));
            end 
            
            end_row = freq_s1d1 - 7;
            start_row = end_row - 22;
            z = z_s1d1(start_row:end_row,:);
            
            for j = 1:time
                for k = 1:23
                    tensor_s1d1(n,j,k) = z(k, j);
                end
            end
        end
        Tensor_s1d1(i,1) = {tensor_s1d1};
    end
    left_train = Tensor_s1d1;
    
    
    
    trail_s1d2 = size(session1_data2,3);
    a_s1d2 = session1_data2(:,:,1);
    mean1_s1d2 = mean(a_s1d2(1,:));
    std1_s1d2 = std(a_s1d2(1,:));
    z1_s1d2 = abs(cwt((a_s1d2(1,:)-mean1_s1d2)/std1_s1d2,wname)); 
    freq_s1d2 = size(z1_s1d2,1);
    tensor_s1d2 = NaN(chan, time, 23);
    Tensor_s1d2 = {};
    for i = 1:trail_s1d2
        b = session1_data2(:,:,i);
        for n= 1:chan
            mean_s1d2 = mean(b(n,:));
            std_s1d2 = std(b(n,:));
            if std_s1d2 == 0
                z_s1d2= abs(cwt(b(n,:)-mean_s1d2,wname));
            else
                z_s1d2= abs(cwt((b(n,:)-mean_s1d2)/std_s1d2,wname));
            end 
            
            end_row = freq_s1d1 - 7;
            start_row = end_row - 22;
            z = z_s1d2(start_row:end_row,:);
            
            
            for j = 1:time
                for k = 1:23
                    tensor_s1d2(n,j,k) = z(k, j);
                end
            end
        end
        Tensor_s1d2(i,1) = {tensor_s1d2};
    end
    right_train = Tensor_s1d2;
    
    
    
    trail_s2d1 = size(session2_data1,3);
    a_s2d1 = session2_data1(:,:,1);
    mean1_s2d1 = mean(a_s2d1(1,:));
    std1_s2d1 = std(a_s2d1(1,:));
    z1_s2d1 = abs(cwt((a_s2d1(1,:)-mean1_s2d1)/std1_s2d1,wname)); 
    freq_s2d1 = size(z1_s2d1,1);
    tensor_s2d1 = NaN(chan, time, 23);
    Tensor_s2d1 = {};
    for i = 1:trail_s2d1
        b = session2_data1(:,:,i);
        for n= 1:chan
            mean_s2d1 = mean(b(n,:));
            std_s2d1 = std(b(n,:));
            if std_s2d1 == 0
                z_s2d1= abs(cwt(b(n,:)-mean_s2d1,wname));
            else
                z_s2d1= abs(cwt((b(n,:)-mean_s2d1)/std_s2d1,wname));
            end 
            
            end_row = freq_s1d1 - 7;
            start_row = end_row - 22;
            z = z_s2d1(start_row:end_row,:);
            
            for j = 1:time
                for k = 1:23
                    tensor_s2d1(n,j,k) = z(k, j);
                end
            end
        end
        Tensor_s2d1(i,1) = {tensor_s2d1};
    end
    left_test = Tensor_s2d1;
    
    
    
    trail_s2d2 = size(session2_data2,3);
    a_s2d2 = session2_data2(:,:,i);
    mean1_s2d2 = mean(a_s2d2(1,:));
    std1_s2d2 = std(a_s2d2(1,:));
    z1_s2d2 = abs(cwt((a_s2d2(1,:)-mean1_s2d2)/std1_s2d2,wname)); 
    freq_s2d2 = size(z1_s2d2,1);
    tensor_s2d2 = NaN(chan, time, 23);
    Tensor_s2d2 = {};
    for i = 1:trail_s2d2
        b = session2_data2(:,:,1);
        for n= 1:chan
            mean_s2d2 = mean(b(n,:));
            std_s2d2 = std(b(n,:));
            if std_s2d2 == 0
                z_s2d2= abs(cwt(b(n,:)-mean_s2d2,wname));
            else
                z_s2d2= abs(cwt((b(n,:)-mean_s2d2)/std_s2d2,wname));
            end 
            
            end_row = freq_s1d1 - 7;
            start_row = end_row - 22;
            z = z_s2d2(start_row:end_row,:);
            
            for j = 1:time
                for k = 1:23
                    tensor_s2d2(n,j,k) = z(k, j);
                end
            end
        end
        Tensor_s2d2(i,1) = {tensor_s2d2};
    end
    right_test = Tensor_s2d2;
end
