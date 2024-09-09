clear;
clc;
%close all;

addpath("D:\MDSI_project\MATLAB\Func");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
% 8,10, 13 in 1 OG shows different pattern
% 8,9,10 in 2 OG show different pattern
list_a = []; 
fs = 1024;
low_freq = 2;
high_freq = 100;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);

k = 1;
figure
for i_file = 1:13
    if ismember(i_file, list_a)
        continue
    end
    
    if i_file <= 7
        i_Acc = 1;
    elseif i_file <= 10
        i_Acc = 7;
    else
        i_Acc = 13;
    end
    load(mat_tile_list{i_Acc});
    %outputSignal = double(timeSeriesData.Data(15,:));
    %outputSignal = filtfilt(b, a, outputSignal);
    inputSignal = double(timeSeriesData.Data(19,:));
    Acc_f = Func_FFT_half(inputSignal,1024);
    Acc_signal = Acc_f.s(Acc_f.f<=500);
    Acc_freq = Acc_f.f(Acc_f.f<=500);
    
    %% Menhir Data for 1 O.G. first mode
    dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer_validation_2OG";
    mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
    fs = 1000;
    load(mat_tile_list{i_file});
    az = data_T.Z_mm_s;
    az_acc = diff(az)/(1/fs);
    [az_psd,freq_psd] = pwelch(az_acc,[],[],[],1000);
    [vz_psd,v_freq_psd] = pwelch(az,[],[],[],1000);
    az_psd = az_psd(freq_psd<=500);
    freq_psd = freq_psd(freq_psd<=500);
    az = filtfilt(b, a, az);
    
    %% FRF
    Menhir_f = Func_FFT_half(transpose(az),1000);
    Men_signal = Menhir_f.s(Menhir_f.f<=500);
    Men_freq = Menhir_f.f(Menhir_f.f<=500);
    
    Men_signal_inter_r = interp1(Men_freq,real(Men_signal),Acc_freq);
    Men_signal_inter_i = interp1(Men_freq,imag(Men_signal),Acc_freq);
    %az_psd_inter = interp1(freq_psd,az_psd,Acc_freq);
    Men_signal_inter = Men_signal_inter_r + 1i*Men_signal_inter_i;
    
    res_Men = Func_PSD_FRF_COH_freq(Acc_freq,Acc_signal,Men_signal_inter);
    
    % Plot the PSD of input and output signals
    %figure;
    subplot(3,1,1);
    plot(res_Men.f, 10*log10(res_Men.Pxx_in));
    hold on;
    plot(res_Men.f, 10*log10(res_Men.Pxx_out));
    title('Power Spectral Density (PSD)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend('Input Signal', 'Output Signal');
    xlim([low_freq,high_freq])
    
    % Plot the Frequency Response Function (FRF)
    %subplot(3,1,2);
    %plot(res_Men.f, abs(res_Men.FRF));
    %title('Frequency Response Function (FRF)');
    %xlabel('Frequency (Hz)');
    %ylabel('Magnitude');
    %xlim([low_freq,high_freq])
    %hold on

    subplot(3,1,2);
    plot(v_freq_psd, vz_psd);
    title('PSD of menhir data (Velocity)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([low_freq,high_freq])
    hold on

    FRF_all(k,:) = vz_psd;
    k = k +1;
    
    %% Plot the coherence between input and output signals
    subplot(3,1,3);
    plot(freq_psd, az_psd);
    title('PSD of menhir data (Acceleration)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    %ylim([0 1]);
    xlim([low_freq,high_freq])
    hold on
end
[l,p] = Func_ConfiPlot(FRF_all,v_freq_psd,1);
xlim([low_freq,high_freq])