%% Visualization of Menhir data for OMA
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")
dir_activate_all = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Event_sort\Event_all_activated";
%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer";
mat_tile_list = Func_FindMatFiles(dir_activate_all);
fs = 1000;
Nmodes = 4;
list = [];
Table_SI = {};   
freq_SI = [];
low_freq = 2; % Lower cutoff frequency (40 Hz)
high_freq = 100; % Upper cutoff frequency (60 Hz)
filterOrder = 4; % Filter order

[b,a] = Func_FilterDesign(low_freq,high_freq,filterOrder,1000);

figure 
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list)
        continue
    end
    load(mat_tile_list{i_file});
    outputSignal = data_T.Z_mm_s_GM_1_all;
    
    outputSignal = filtfilt(b, a, outputSignal);
    inputSingal = data_T.Z_mm_s_GM_0_all;
    %inputSingal = filtfilt(b, a, inputSingal);

    out_S_f = Func_FFT_half(transpose(outputSignal),1000);
    in_S_f = Func_FFT_half(transpose(inputSingal),1000);

    %out_S_f.s = movmean(out_S_f.s,20);
    %in_S_f.s = movmean(in_S_f.s,20);
    %FRF_fft = out_S_f.s./in_S_f.s;

    res_Men = Func_PSD_FRF_COH(inputSingal,outputSignal,[],[],[],1000);
    
    subplot(3,1,1);
    plot(res_Men.f, 10*log10(res_Men.Pxx_in),"Color",'r');
    hold on;
    plot(res_Men.f, 10*log10(res_Men.Pxx_out),"Color",'b');
    
    title('Power Spectral Density (PSD)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend('Input Signal', 'Output Signal');
    xlim([low_freq,high_freq])
    
    % Plot the Frequency res_Menponse Function (FRF)
    subplot(3,1,2);
    plot(res_Men.f, abs(res_Men.FRF)/max(abs(res_Men.FRF)));
    hold on
    title('Frequency resenponse Function (FRF)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([low_freq,high_freq])
    
    % Plot the coherence between input and output signals
    subplot(3,1,3);
    %plot(out_S_f.f,abs(out_S_f.s)/max(abs(out_S_f.s)));
    plot(res_Men.f, res_Men.Cxy);
    %yline(0.9)
    hold on
    title('Coherence');
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    %ylim([0 1]);
    xlim([low_freq,high_freq])
    
end