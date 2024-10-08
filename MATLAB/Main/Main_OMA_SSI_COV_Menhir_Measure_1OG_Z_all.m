clear;
clc;
%close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-SSICOV-82ce27a");

%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_all = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Event_sort\Event_all_activated";
mat_tile_list = Func_FindMatFiles(dir_activate_all);
list_a = [10,11];
fs = 1000;
low_freq = 2;
high_freq = 35;
%[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);
[b,a] = Func_FilterDesign_highlow('low',high_freq,4,1000);

All_F_res = {};
All_F_damp = {};

for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list_a)
        continue
    end

    load(mat_tile_list{i_file});
    %% For 1 O.G
    S_0 = data_T.Z_mm_s_GM_0_all; 
    S_1 = data_T.Z_mm_s_GM_1_all;  
    S_2 = data_T.Z_mm_s_GM_2_all; 

    %% For 2 O.G
    %S_1 = double(timeSeriesData.Data(1,:));
    %S_2 = double(timeSeriesData.Data(2,:));
    %S_3 = double(timeSeriesData.Data(3,:));
    %S_4 = double(timeSeriesData.Data(15,:));
    
    S_0 = filtfilt(b, a, S_0);
    S_1 = filtfilt(b, a, S_1);
    S_2 = filtfilt(b, a, S_2);

    S_0 = Func_Resample(S_0,1000,100);
    S_1 = Func_Resample(S_1,1000,100);
    S_2 = Func_Resample(S_2,1000,100);

    %rz = [transpose(S_0);transpose(S_1);transpose(S_2)];
    rz = [transpose(S_0);transpose(S_1)];
    dt = 1/200;
            
    [Nyy,N]= size(rz); %Nyy is the number of sensors and N the number of time step
    %
    [fn0,zeta0,~,~] = SSICOV(rz,dt,'methodCOV',1,'Nmin',2,'Nmax',20,'eps_cluster',0.2);
    
    All_F_res{i_file}= fn0;
    All_F_damp{i_file}= zeta0;
end

x_axis = [];
y_axis = [];
z_axis = [];
%list_a = [10,11];
k = 1;
for i_file = 1:length(All_F_res)
    if ismember(i_file, list_a)
        continue
    end
    x_axis = [x_axis , All_F_res{i_file}];
    y_axis = [y_axis , All_F_damp{i_file}];
    z_axis = [z_axis , k*ones(1,length(All_F_res{i_file}))];
    k = k+1;
end

scatter3(x_axis,y_axis,z_axis,30,x_axis,'filled','o','MarkerEdgeColor','flat','MarkerFaceAlpha',0.4)
grid on
grid minor
xlim([0,high_freq]);
xlabel('$f_{r}$ (Hz)', 'Interpreter', 'latex');
ylabel('$\zeta $', 'Interpreter', 'latex');
zlabel('$i_{th}$ testing', 'Interpreter', 'latex');


%[h] = plotStabDiag(paraPlot.fn,rz(2,:),fs,paraPlot.status,paraPlot.Nmin,paraPlot.Nmax);
