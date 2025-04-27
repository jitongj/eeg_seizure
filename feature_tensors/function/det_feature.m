function det_feature(patient,win_num,sz_det_win,nonsz_det_win)
feature_det = {};
for k = 1:(win_num*2)
    if k <(win_num+1)
        a = sz_det_win(k).window;
    else
        a = nonsz_det_win(k-win_num).window;
    end

    tensor = NaN(size(a,1),4,5,4);%chan*wave(db1-db4)*feature*level(1-4)

    for n = 1:size(a,1)%chan
        b= a(n,:);
        for i=1:4 %wave
            wavename=sprintf('db%02d',i);
            for j = 1:4%level
                [c,l] = wavedec(b,j,wavename);
                approx = appcoef(c,l,wavename);
            
                %mean
                mean_ca = mean(approx);
        
                %shannon entropy
                shen_ca = wentropy(approx,'shannon');
        
                %kurtosis
                kurt_ca = kurtosis(approx);

                %standard deviation
                std_ca = std(approx);
                
                %zero-crossing
                zc_ca = zerocrossrate(approx);
          
                tensor(n,i,:,j) = [mean_ca;shen_ca;kurt_ca;std_ca;zc_ca];
    
            end
        
        end
    end

    feature_det(k,1) = {tensor};
    if k <(win_num+1)
        feature_det(k,2) = {'seizure'};
    else
        feature_det(k,2) = {'non_seizure'};
    end
end

filename=sprintf('Feat_det%02d',patient);
save(filename,'feature_det');
end

