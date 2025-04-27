function fft_SeiPred_tensor(patient,sz_pred_win,SamplingRate,rate)
fft_SeiPred_tensor = {};
for i =1: size(sz_pred_win,2)
    FilteredChannelData1 = sz_pred_win(i).window(1,:);
    Demeaned_FilteredlData1 = FilteredChannelData1 - mean(FilteredChannelData1);

    BinsSec1     = 2 ; %in s
    DesiredRoll1 = 100; %in ms

    %Convert DesiredRoll to s
    Rolling1 = DesiredRoll1*0.001;
    %Calculate the Total Time in the recording
    TotalTime1 = floor(length(Demeaned_FilteredlData1)/SamplingRate);
    %Calculate the middle time value for each bin
    tmid1 = (0+(BinsSec1/2)):Rolling1:(TotalTime1-(BinsSec1/2));
    %Define the frequency                      
    Fs1 = SamplingRate;
    spec1 = [];
    %Define the indexes used to generate the bins
    Index1 = 1:(Rolling1*SamplingRate):(length(Demeaned_FilteredlData1));

    for ii = 1:(length(tmid1)-1)
        %Select the data segment that is going to be used to calculate the FFT
        LeftEnd1  =  Index1(ii);
        RightEnd1 =  Index1(ii+(BinsSec1/Rolling1));
        x1 = Demeaned_FilteredlData1(LeftEnd1:RightEnd1);
        %Calculate the FFT for the desired time segment
    
        %Calculate the hamming window
        N1 = length(x1);
        ham1 = hamming(N1);
        %Multiply the data by the hamming window
        hx1 = ham1.*x1;
        %Perform the FFT on the data
        xdft1 = fft(hx1);
        xdft1 = xdft1(1:round(N1/2)+1);
        psdx1 = (1/(Fs1*N1)) * abs(xdft1).^2;
        psdx1(2:end-1) = 2*psdx1(2:end-1);
        freq1 = 0:Fs1/length(hx1):Fs1/2;
        spec1(:,ii) = 10*log10(psdx1);
    end

    y1 = spec1(1:50,:);
    y1 = downsample(y1.',rate);
    y1 = y1.';


    time = size(y1,2);
    freq = size(y1,1);
    chan = size(sz_pred_win(i).window,1);



    % Build the 3d tensor
    tensor = NaN(chan, time, freq);

    for n= 1:chan

        FilteredChannelData = sz_pred_win(i).window(n,:);
        Demeaned_FilteredlData = FilteredChannelData - mean(FilteredChannelData);

        BinsSec = 2 ; %in s
        DesiredRoll = 100; %in ms

        %Convert DesiredRoll to s
        Rolling = DesiredRoll*0.001;
        %Calculate the Total Time in the recording
        TotalTime = floor(length(Demeaned_FilteredlData)/SamplingRate);
        %Calculate the middle time value for each bin
        tmid = (0+(BinsSec/2)):Rolling:(TotalTime-(BinsSec/2));
        %Define the frequency                      
        Fs = SamplingRate;
        spec = [];
        %Define the indexes used to generate the bins
        Index = 1:(Rolling*SamplingRate):(length(Demeaned_FilteredlData));

        for ii = 1:(length(tmid)-1)
            %Select the data segment that is going to be used to calculate the FFT
            LeftEnd  =  Index(ii);
            RightEnd =  Index(ii+(BinsSec/Rolling));
            x = Demeaned_FilteredlData(LeftEnd:RightEnd);
            %Calculate the FFT for the desired time segment
    
            %Calculate the hamming window
            N = length(x);
            ham = hamming(N);
            %Multiply the data by the hamming window
            hx = ham.*x;
            %Perform the FFT on the data
            xdft = fft(hx);
            xdft = xdft(1:round(N/2)+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            freq = 0:Fs/length(hx):Fs/2;
            spec(:,ii) = 10*log10(psdx);
        end

        y = spec(1:50,:);
        y = downsample(y.',rate);
        y = y.';
 
        % Fill the tensor
     
        tensor(n,:,:) = y.';
    end

    fft_SeiPred_tensor(i,1) = {tensor};
end               

filename=sprintf('fft_SeiPred_tensor02d',patient);
save(filename,'fft_SeiPred_tensor');
end


