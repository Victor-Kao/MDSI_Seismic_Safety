clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
%load("D:/MDSI_project/DATA_GM_RawData/CellArrayAllEvent.mat");

%%
plot_ = 1;

%% 
load("Menhir_result_1G_Z.mat");

cut_f = 150;
% Apply DBSCAN
f  = linspace(0,500,2501);
All_s = All_s(:,f<=cut_f);
All_s_smooth = movmean(All_s,5,2);
norm_All_s = normalize(All_s_smooth,2,"range");
f_s = f(f<=cut_f);

hc = Func_HC_clustering(All_s,f_s,19,1,1);

%e = trapz(All_s_smooth,2);
%s_e = All_s_smooth./e;
%s_e_e = max(s_e,[],2);
%histogram(s_e_e);