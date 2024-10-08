clear;
clc;
close all;

%% Visualization of Menhir data for OMA
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")
addpath("D:\MDSI_project\MATLAB\Func")
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Event_sort\Event_all_activated";
%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);

figure 
k = 1;
low_freq = 0;
high_freq = 100;
for i_file = 1:length(mat_tile_list)

    load(mat_tile_list{i_file});
    FRF = Func_PSD_FRF_COH(data_T.Z_mm_s_GM_0_all,data_T.Z_mm_s_GM_1_all,[],[],5000,1000);
        %figure;


    subplot(3,1,1);
    plot(FRF.f, 10*log10(FRF.Pxx_in),'Color',[0.8500 0.3250 0.0980]);
    hold on;
    plot(FRF.f, 10*log10(FRF.Pxx_out),'Color',[0 0.4470 0.7410]);
    title('Power Spectral Density $S_{x}(f)$', 'Interpreter', 'latex');
    xlabel('Freq (Hz)', 'Interpreter', 'latex');
    ylabel('$S_{x}(f)$(DB)', 'Interpreter', 'latex');
    legend('Excitation', 'Response', 'Interpreter', 'latex');
    xlim([low_freq,high_freq])
    
    % Plot the Frequency Response Function (FRF)
    subplot(3,1,2);
    plot(FRF.f, abs(FRF.FRF));
    
  
    FRF_all(k,:) = abs(FRF.FRF);
    k = k +1;

    title('Frequency Response Function (FRF)', 'Interpreter', 'latex');
    xlabel('Freq (Hz)', 'Interpreter', 'latex');
    ylabel('Magnitude', 'Interpreter', 'latex');
    xlim([low_freq,high_freq])
    hold on
    
    % Plot the coherence between input and output signals
    subplot(3,1,3);
    plot(FRF.f, FRF.Cxy);
    yline(0.9)
    title('Coherence', 'Interpreter', 'latex');
    xlabel('Freq (Hz)', 'Interpreter', 'latex');
    ylabel('Coherence', 'Interpreter', 'latex');
    ylim([0 1]);
    xlim([low_freq,high_freq])
    hold on
end