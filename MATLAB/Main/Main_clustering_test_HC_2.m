clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
load("D:/MDSI_project/DATA_GM_RawData/CellArrayAllEvent.mat");

%%
plot_ = 0;

load("Events_all_1G_X.mat");
Event_All_Floor = Events_GM;
load("Events_1G2G_1G_X.mat");
Event_1G2G = Events_GM;
Events_GM = [Event_All_Floor;Event_1G2G] ;


%load("Menhir_result_1G2G_ALL.mat");

cut_f = 150;
f  = linspace(0,500,2501);
All_s = Events_GM(:,f<=cut_f);
f_s = f(f<=cut_f);

hc = Func_HC_clustering(All_s,f_s,0,1,1);

%index = find(hc==6);
%for i = 1:length(index)
%    plot(f_s,All_s(index(i),:))
%    hold on 
%end







