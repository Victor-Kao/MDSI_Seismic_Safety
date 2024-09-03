clear;
clc;
%close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-SSICOV-82ce27a");

dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Event_sort\Event_1G2G_activated";

mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
fs = 1000;
Nmodes = 4;

lowCutoffFreq = 4; % Lower cutoff frequency (40 Hz)
highCutoffFreq = 22; % Upper cutoff frequency (60 Hz)
filterOrder = 4; % Filter order


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
low_freq = 2;
high_freq = 40;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);
i_file = 2;

load(mat_tile_list{i_file});
S_1 = double(timeSeriesData.Data(10,:));
S_2 = double(timeSeriesData.Data(11,:));
S_3 = double(timeSeriesData.Data(12,:));
S_4 = double(timeSeriesData.Data(9,:));

S_1 = filtfilt(b, a, S_1);
S_2 = filtfilt(b, a, S_2);
S_3 = filtfilt(b, a, S_3);
S_4 = filtfilt(b, a, S_4);

S_3_f = Func_FFT_half(S_3,1024);
S_4_f = Func_FFT_half(S_4,1024);

S_1 = Func_Resample(S_1,1024,128);
S_2 = Func_Resample(S_2,1024,128);
S_3 = Func_Resample(S_3,1024,128);
S_4 = Func_Resample(S_4,1024,128);


plot(S_3_f.f,abs(S_3_f.s))
hold on 
plot(S_4_f.f,abs(S_4_f.s))
rz = [S_1;S_2;S_3;S_4];
dt = 1/256;
Ts = 8;
        
[Nyy,N]= size(rz); %Nyy is the number of sensors and N the number of time step
       
tic
[fn0,zeta0,phi0,paraPlot] = SSICOV(rz,dt,'Ts',Ts,'Nmin',1,'Nmax',30,'eps_cluster',0.05,'eps_freq',1e-2,'eps_MAC',1e-2);
toc
disp(fn0)
disp(zeta0)
[h] = plotStabDiag(paraPlot.fn,rz(2,:),fs,paraPlot.status,paraPlot.Nmin,paraPlot.Nmax);