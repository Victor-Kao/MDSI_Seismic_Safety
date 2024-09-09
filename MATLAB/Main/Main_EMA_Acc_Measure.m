clear;
clc;
%close all;

addpath("D:\MDSI_project\MATLAB\Func");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
low_freq = 2;
high_freq = 100;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);


k = 1;
figure
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list_a)
        continue
    end

    load(mat_tile_list{i_file});
    outputSignal = double(timeSeriesData.Data(9,:));
    %outputSignal = filtfilt(b, a, outputSignal);
    inputSignal = double(timeSeriesData.Data(19,:));
    
    % Compute the Power Spectral Density (PSD) of the input signal
    [pxxInput, f] = pwelch(inputSignal, [], [], [], fs);
    % Compute the Cross-Power Spectral Density (CSD) between input and output
    [pxy, f] = cpsd(inputSignal, outputSignal, [], [], [], fs);
    
    % Compute the Frequency Response Function (FRF)
    FRF = pxy ./ pxxInput;
    
    % Compute the PSD of the output signal for comparison
    [pxxOutput, f] = pwelch(outputSignal, [], [], [], fs);

    % Compute the coherence between input and output signals
    [Cxy, f] = mscohere(inputSignal, outputSignal, [], [], [], fs);
    
    % Plot the PSD of input and output signals
    %figure;
    subplot(3,1,1);
    plot(f, 10*log10(pxxInput));
    hold on;
    plot(f, 10*log10(pxxOutput));
    title('Power Spectral Density (PSD)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend('Input Signal', 'Output Signal');
    xlim([low_freq,high_freq])
    
    % Plot the Frequency Response Function (FRF)
    subplot(3,1,2);
    plot(f, abs(FRF));
    
  
    FRF_all(k,:) = abs(FRF);
    k = k +1;

    title('Frequency Response Function (FRF)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([low_freq,high_freq])
    hold on
    
    % Plot the coherence between input and output signals
    subplot(3,1,3);
    plot(f, Cxy);
    yline(0.9)
    title('Coherence');
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    ylim([0 1]);
    xlim([low_freq,high_freq])
    hold on
end


[l,p] = Func_ConfiPlot(FRF_all,f,1);
xlim([low_freq,high_freq])