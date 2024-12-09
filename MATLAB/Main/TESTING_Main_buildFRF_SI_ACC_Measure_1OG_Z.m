clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\pickpeaks");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
cutoff_freq = 40;
%[b,a] = Func_FilterDesign_highlow('low',cutoff_freq,4,1024);
[b,a] = Func_FilterDesign(4,cutoff_freq,4,1024);

for i_file = 1:length(mat_tile_list)

    %i_file = 15;

    if ismember(i_file,list_a )
        continue;
    end
    
    load(mat_tile_list{i_file});
    outputSignal = double(timeSeriesData.Data(9,:));
    outputSignal = filtfilt(b, a, outputSignal);
    t = timeSeriesData.Time;
    inputSignal = double(timeSeriesData.Data(19,:));
    
    
    [FRF_ROM,freq_ROM,SI_info] = Func_Proposed_SI_FRF(inputSignal,outputSignal,t,fs,cutoff_freq,2,0.5);
    Wn = SI_info(SI_info(:,2)>17,2);
    [pxx,f_psd] = pwelch(outputSignal,[],[],[],fs);
    %[pxy, ~] = cpsd(inputSignal, outputSignal,[],[],[], fs);
    %FRF_ = pxy ./ pxx;
    powerbw(pxx,f_psd)/(2*unique(Wn))  
    %for j = 1:size(SI_info,1)/2
    %    k = 2*j-1;
    %    scatter(SI_info(k,2),SI_info(k+1,1)) 
    %    hold on 
    %end
    
    %grid on
    % Compute the FFT of the input and output signals
    N = length(t);
    f = (0:N-1)*(fs/N); % Frequency vector
    
    inputFFT = fft(inputSignal);
    outputFFT = fft(outputSignal);
    
    %inputFFT  = inputFFT(f<=cutoff_freq);
    %outputFFT  = outputFFT(f<=cutoff_freq);
    
    % Compute the frequency response (H = output/input)
    H = outputFFT ./ inputFFT;
    H(1) = 0 + 1i*0;
    H_1 = H(f<=cutoff_freq);
    f_1 = f(f<=cutoff_freq);

    
    figure 
    plot(freq_ROM,abs(FRF_ROM))
    hold on 
    plot(f_1,abs(H_1))

end

   


