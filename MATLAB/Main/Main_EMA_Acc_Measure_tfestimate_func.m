clear;
clc;
%close all;

addpath("D:\MDSI_project\MATLAB\Func");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
low_freq = 6;
high_freq = 100;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);
i_file = 9;
figure
list_a = [1,2,4,5,6,7,8];
for i_pos = 1:15%length(mat_tile_list)
    if ismember(i_pos, list_a)
        continue
    end

    load(mat_tile_list{i_file});
    outputSignal = double(timeSeriesData.Data(i_pos,:));
    %outputSignal = filtfilt(b, a, outputSignal);
    inputSignal = double(timeSeriesData.Data(19,:));
    
    % Compute the Frequency Response Function (FRF) using tfestimate
    [FRF, f] = tfestimate(inputSignal, outputSignal, [], [], 5000, fs);
    FRF = FRF./(-((2 * pi * f).^2));
    % Compute the phase of the FRF
    phase = angle(FRF);
    
    % Plot the Magnitude of the FRF
    
    subplot(2, 1, 1);
    plot(f, abs(FRF));
    hold on
    title('Frequency Response Function (FRF) - real');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([low_freq,high_freq])
    grid on;
    
    % Plot the Phase of the FRF
    subplot(2, 1, 2);
    plot(f, imag(FRF));
    hold on
    title('Frequency Response Function (FRF) - imag');
    xlabel('Frequency (Hz)');
    ylabel('Phase (radians)');
    xlim([low_freq,high_freq])
    grid on;

    freq_ = f;
    real_ = real(FRF);
    imag_ = imag(FRF);
    %save(['D:/MDSI_project/MATLAB/Surrogate_main/FRF/FRF_test_',num2str(i_file),'_ch_',num2str(i_pos),'.mat'], 'freq_', 'real_', 'imag_');
end