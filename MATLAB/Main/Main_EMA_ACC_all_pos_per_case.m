clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\pickpeaks");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain_update";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [];
fs = 1024;
low_freq = 5;
high_freq = 100;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);


k = 1;
figure
sgtitle('EMA using Impact Hammer ', 'Interpreter', 'latex','FontSize',25)

i_file = 9;
list_a = [1,2,4,5,6,7,8];
for i_pos = 1:15%length(mat_tile_list)
    if ismember(i_pos, list_a)
        continue
    end

    load(mat_tile_list{i_file});
    outputSignal = double(timeSeriesData.Data(i_pos,:));
    %outputSignal = filtfilt(b, a, outputSignal);
    inputSignal = double(timeSeriesData.Data(19,:));

    
    res_Men = Func_PSD_FRF_COH(inputSignal, outputSignal,[],[],5000,1024);
    
    low_freq_ = 5;
    high_freq_ = 25;
    FRF = res_Men.FRF(res_Men.f <=high_freq_);
    f = res_Men.f(res_Men.f <=high_freq_);
    FRF = FRF(f >=low_freq_);
    FRF = abs(FRF);
    f = f(f>=low_freq_);
    [peak_loc,~] = pickpeaks(abs(FRF),2,0);
    f(peak_loc)
    

    if ismember(i_pos,[9,10,11,12]) 
    subplot(2,1,1);
    % Plot the Frequency Response Function (FRF)
    plot(res_Men.f, imag(res_Men.FRF),'DisplayName', ['ID:',num2str(i_pos)]);
    hold on
    k = k +1;
    title('FRF in displacement 1.OG', 'Interpreter', 'latex','FontSize',14);
    xlabel('Freq (Hz)', 'Interpreter', 'latex','FontSize',14);
    ylabel('Magnitude', 'Interpreter', 'latex','FontSize',14);
    xlim([low_freq,high_freq])
    legend('show');  % Show legend
    yline(0)
    else
    subplot(2,1,2);
    plot(res_Men.f, imag(res_Men.FRF),'DisplayName', ['ID:',num2str(i_pos)]);
    hold on
    k = k +1;
    title('FRF in displacement 2.OG', 'Interpreter', 'latex','FontSize',14);
    xlabel('Freq (Hz)', 'Interpreter', 'latex','FontSize',14);
    ylabel('Magnitude', 'Interpreter', 'latex','FontSize',14);
    xlim([low_freq,high_freq])
    legend('show');  % Show legend
    yline(0)
    end
    
    freq_ = res_Men.f;
    real_ = real(res_Men.FRF);
    imag_ = imag(res_Men.FRF);
    %save(['D:/MDSI_project/MATLAB/Surrogate_main/FRF/FRF_test_',num2str(i_file),'_ch_',num2str(i_pos),'.mat'], 'freq_', 'real_', 'imag_');
end


xlim([low_freq,high_freq])