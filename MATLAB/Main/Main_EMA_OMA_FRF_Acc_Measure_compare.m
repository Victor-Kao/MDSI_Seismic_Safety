clear;
clc;
%close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
low_freq = 2;
high_freq = 40;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);

figure
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list_a)
        continue
    end

    load(mat_tile_list{i_file});
    outputSignal = double(timeSeriesData.Data(9,:));
    outputSignal_f = Func_FFT_half(outputSignal,1024);
    outputSignal = filtfilt(b, a, outputSignal);
    inputSignal = double(timeSeriesData.Data(19,:));

    %% EMA method
    res = Func_PSD_FRF_COH(inputSignal,outputSignal,[],[],[],fs);
    [res_1, freq] = tfestimate(inputSignal, outputSignal, [], [], [], fs);

    %% OMA methd
    [Saz,f] = pwelch(detrend(outputSignal),[],[],[],fs);         
    Nperiod = 60; % The IRF is computed for 30 period on each mode
    [fn,zeta] = modalID_singleSensor(outputSignal,Saz,f,3,fs,...
        'plotOpt',0,'Nperiod',Nperiod,'PickingMethod','auto');
    %zeta(1) = 0.04;
    %zeta(2) = 0.0197;
    [FRF_OMA,freq_OMA] = Func_buildFRF(fn,zeta,high_freq,fs);
    

    %plot(res.f,max(abs(FRF_OMA))*abs(res.FRF)/max(abs(res.FRF)),'r')
    plot(res.f,abs(res.FRF)/max(abs(res.FRF)),'r')
    hold on
    %plot(freq,abs(res_1),'--')
    plot(freq_OMA,abs(FRF_OMA)/max(abs(FRF_OMA)),'b')
    %plot(outputSignal_f.f,max(abs(FRF_OMA))*abs(outputSignal_f.s)/max(abs(outputSignal_f.s)))

    xlim([low_freq,high_freq])
end