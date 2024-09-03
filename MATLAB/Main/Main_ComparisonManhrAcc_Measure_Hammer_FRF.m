clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer_validation_1OG")
addpath("D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer_validation_2OG")

HammerSignal_path = "D:\MDSI_project\DATA_GM_RawData\hammer_menhir.mat";
load(HammerSignal_path)

list_a = [];
fs = 1024;
low_freq = 5;
high_freq = 40;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);

size_table = size(fileTable);
for i = 1:size_table(1)
    if ismember(i, list_a)
        continue
    end
    
    Acc_data_name = fileTable.FileNames(i);
    file_path = fileTable.Var2(i);
    load(string(file_path));

    index = find(strcmp(t_data(:,1), Acc_data_name));
    Acc_data = t_data(index,2);

    % Import input data 
    Acc_input = Acc_data{1}.timeSeriesData.Data(19,:);
    Acc_f = Func_FFT_half(Acc_input,1024);
    Acc_signal = Acc_f.s(Acc_f.f<=50);
    Acc_freq = Acc_f.f(Acc_f.f<=50);
    
    % Import output data
    outputSignal = data_T.Z_mm_s;
    Menhir_f = Func_FFT_half(transpose(outputSignal),1000);
    Men_signal = Menhir_f.s(Menhir_f.f<=50);
    Men_freq = Menhir_f.f(Menhir_f.f<=50);
    
    Men_signal_inter_r = interp1(Men_freq,real(Men_signal),Acc_freq);
    Men_signal_inter_i = interp1(Men_freq,imag(Men_signal),Acc_freq);
    %az_psd_inter = interp1(freq_psd,az_psd,Acc_freq);
    Men_signal_inter = Men_signal_inter_r + 1i*Men_signal_inter_i;

    Acc_output = Acc_data{1}.timeSeriesData.Data(9,:);
    Acc_output = filtfilt(b, a, double(Acc_output));
    %Acc_input = filtfilt(b, a, double(Acc_input));

    % FRF 
    FRF_Acc = Func_PSD_FRF_COH(Acc_input,Acc_output,[],[],[],1024);
    FRF_Acc_f = FRF_Acc.FRF(FRF_Acc.f<=50);
    Acc_cut_f = FRF_Acc.f(FRF_Acc.f<=50);
    FRF_Menhir = Func_PSD_FRF_COH_freq(Acc_freq,Acc_signal,Men_signal_inter);

    figure
    plot(Acc_cut_f,abs(FRF_Acc_f)/max(abs(FRF_Acc_f)),'r');
    hold on 
    plot(FRF_Menhir.f, abs(FRF_Menhir.FRF) / max(abs(FRF_Menhir.FRF)),'b' );
    xlim([0,50])
end