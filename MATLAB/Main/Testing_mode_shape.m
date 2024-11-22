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


i_file = 13;
list_a = [9,10,11,12,3,13,14,15];
mode_shape_vector = zeros(16,4);
for i = 1:8%length(mat_tile_list)
    
    i_pos = list_a(i);
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
    f = f(f>=low_freq_);
    [peak_loc,~] = pickpeaks(abs(FRF),2,0);
    peak_loc = sort(peak_loc);
    f_peak = f(peak_loc);
    mode_shape = imag(FRF);
    mode_shape_peak = mode_shape(peak_loc);
    mode_shape_vector(2*i-1,1) = f_peak(1);
    mode_shape_vector(2*i,1) = f_peak(2);
    mode_shape_vector(2*i-1,2) = mode_shape_peak(1);
    mode_shape_vector(2*i,2) = mode_shape_peak(2);
end
mode_shape_vector(:,3) = mode_shape_vector(:,2)./max(abs(mode_shape_vector(:,2)));
max_1 = max(abs(mode_shape_vector(1:8,2)));
max_2 = max(abs(mode_shape_vector(9:16,2)));
mode_shape_vector(1:8,4) = mode_shape_vector(1:8,2)./max_1;
mode_shape_vector(9:16,4) = mode_shape_vector(9:16,2)./max_2;
phi = mode_shape_vector;
%save(['D:/MDSI_project/MATLAB/Surrogate_main/FRF/mode_shape_test_',num2str(i_file),'.mat'], 'phi');