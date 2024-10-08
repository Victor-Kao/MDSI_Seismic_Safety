clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-SSICOV-82ce27a");


%% Acc Data for 2 O.G. (Channel 15), first mode
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

[b,a] = Func_FilterDesign_highlow('high',40,2,1000);
i_file = 20;

load(mat_tile_list{i_file});
S_1 = data_T.Z_mm_s_GM_0_all;
S_2 = data_T.Z_mm_s_GM_1_all;  
S_3 = data_T.Z_mm_s_GM_2_all;  

S_1 = filtfilt(b, a, S_1);
S_2 = filtfilt(b, a, S_2);
S_3 = filtfilt(b, a, S_3);

S_1_f = Func_FFT_half(transpose(S_1),1000);
S_2_f = Func_FFT_half(transpose(S_2),1000);
S_3_f = Func_FFT_half(transpose(S_3),1000);

S_1 = Func_Resample(S_1,1000,100);
S_2 = Func_Resample(S_2,1000,100);
S_3 = Func_Resample(S_3,1000,100);

figure
plot(S_2_f.f,abs(S_2_f.s))
hold on 
plot(S_3_f.f,abs(S_3_f.s))
plot(S_1_f.f,abs(S_1_f.s))
xlim([0,high_freq])
rz = transpose([S_1,S_2,S_3]);
dt = 1/200;
Ts = 10;
        
[Nyy,N]= size(rz); %Nyy is the number of sensors and N the number of time step
       
tic
[fn0,zeta0,phi0,paraPlot] = SSICOV(rz,dt,'Ts',Ts,'Nmin',1,'Nmax',60,'eps_cluster',0.05,'eps_freq',1e-2,'eps_MAC',1e-2);
toc
disp(fn0)
disp(zeta0)
[h] = plotStabDiag(paraPlot.fn,rz(2,:),fs,paraPlot.status,paraPlot.Nmin,paraPlot.Nmax);