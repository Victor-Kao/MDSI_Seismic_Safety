clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");


low_freq = 2;
high_freq = 50;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1000);

%i_file = 5;
list_a = [10,11];

figure;
for i_file = 1:13
    if ismember(i_file, list_a)
        continue
    end
    %% Acc Data for 2 O.G. (Channel 15), first mode
    dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
    mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
    fs = 1024;
    load(mat_tile_list{i_file});
    outputSignal = double(timeSeriesData.Data(9,:));
    outputSignal = filtfilt(b, a, outputSignal);
    inputSignal = double(timeSeriesData.Data(19,:));
    Acc_f = Func_FFT_half(inputSignal,fs);
    Acc_signal = Acc_f.s(Acc_f.f<=500);
    Acc_freq = Acc_f.f(Acc_f.f<=500);
    res_Acc = Func_PSD_FRF_COH(inputSignal,outputSignal,[],[],[],1024);
    
    
    %% Menhir Data for 1 O.G. first mode
    dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer_validation_2OG";
    mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
    fs = 1000;
    load(mat_tile_list{i_file});
    az = diff(data_T.Z_mm_s)/(1/fs); % WARNING! Identify based on Velocity not ACC here 
    %az = data_T.Z_mm_s;
    az = filtfilt(b, a, az);
    [az_psd,az_freq] = pwelch(az,[],[],[],1000);
    
    %% FRF
    %Menhir_f = Func_FFT_half(transpose(az),1000);
    %Men_signal = Menhir_f.s(Menhir_f.f<=500);
    %Men_freq = Menhir_f.f(Menhir_f.f<=500);
    %Men_signal_inter_r = interp1(Men_freq,real(Men_signal),Acc_freq);
    %Men_signal_inter_i = interp1(Men_freq,imag(Men_signal),Acc_freq);
    %Men_signal_inter = Men_signal_inter_r + 1i*Men_signal_inter_i;
    %res_Men = Func_PSD_FRF_COH_freq(Acc_freq,Acc_signal,Men_signal_inter);
    
    subplot(3,1,1);
    plot(res_Acc.f, 10*log10(res_Acc.Pxx_in));
    hold on
    plot(res_Acc.f, 10*log10(res_Acc.Pxx_out));
    title('Power Spectral Density (PSD)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend('Input Signal','Acc Signal');
    xlim([low_freq,high_freq])
    
    % Plot the Frequency res_Menponse Function (FRF)
    subplot(3,1,2);
    %plot(res_Men.f, abs(res_Men.FRF)/max(abs(res_Men.FRF)));
    plot(az_freq,az_psd/max(abs(az_psd)),'r')
    hold on 
    plot(res_Acc.f, abs(res_Acc.FRF)/max(abs(res_Acc.FRF)),'b');
    legend('Menhir Signal','Acc Signal');
    title('Frequency resenponse Function (FRF)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([low_freq,high_freq])
    
    % Plot the coherence between input and output signals
    subplot(3,1,3);
    %plot(res_Men.f, res_Men.Cxy);
    %hold on
    plot(res_Acc.f, res_Acc.Cxy);
    hold on
    yline(0.9)
    legend('Menhir Signal','Acc Signal');
    title('Coherence');
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    ylim([0 1]);
    xlim([low_freq,high_freq])
end






