clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain_update";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [];
fs = 1024;
low_freq = 2;
high_freq = 100;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);


k = 1;
figure
sgtitle('EMA using Impact Hammer ', 'Interpreter', 'latex','FontSize',25)
for i_file = 10:13%length(mat_tile_list)
    if ismember(i_file, list_a)
        continue
    end

    load(mat_tile_list{i_file});
    outputSignal = double(timeSeriesData.Data(14,:));
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
    
    res_Men = Func_PSD_FRF_COH(inputSignal, outputSignal,[],[],5000,1024);

    % Plot the PSD of input and output signals
    %figure;
    subplot(3,1,1);
    plot(res_Men.f, 20*log10(res_Men.Pxx_in),"Color",'r');
    hold on;
    plot(res_Men.f, 20*log10(res_Men.Pxx_out),"Color",'b');
    title('Power Spectral Density $S_{x}(f)$', 'Interpreter', 'latex','FontSize',14);
    xlabel('Freq (Hz)', 'Interpreter', 'latex','FontSize',14);
    ylabel('$S_{x}(f)$(DB)', 'Interpreter', 'latex','FontSize',14);
    legend('Excitation', 'Response', 'Interpreter', 'latex','FontSize',14);
    xlim([low_freq,high_freq])
    
    % Plot the Frequency Response Function (FRF)
    subplot(3,1,2);
    plot(res_Men.f, abs(res_Men.FRF));
    hold on
    
  
    FRF_all(k,:) = 20*log10(abs(FRF));
    k = k +1;

    title('Frequency Response Function (FRF in disp)', 'Interpreter', 'latex','FontSize',14);
    xlabel('Freq (Hz)', 'Interpreter', 'latex','FontSize',14);
    ylabel('Magnitude', 'Interpreter', 'latex','FontSize',14);
    xlim([low_freq,high_freq])
    hold on
    
    % Plot the coherence between input and output signals
    subplot(3,1,3);
    %plot(out_S_f.f,abs(out_S_f.s)/max(abs(out_S_f.s)));
    plot(res_Men.f, res_Men.Cxy);
    yline(0.9)
    title('Coherence', 'Interpreter', 'latex','FontSize',14);
    xlabel('Freq (Hz)', 'Interpreter', 'latex','FontSize',14);
    ylabel('Coherence', 'Interpreter', 'latex','FontSize',14);
    ylim([0 1]);
    xlim([low_freq,high_freq])
    hold on
end


[l,p] = Func_ConfiPlot(FRF_all,f,1);
xlim([low_freq,high_freq])