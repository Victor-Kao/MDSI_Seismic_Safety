clear;
clc;
%close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-SSICOV-82ce27a");

%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain_update";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [];
fs = 1024;
low_freq = 2;
high_freq = 40;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);
%[b,a] = Func_FilterDesign_highlow('low',high_freq,4,1024);

All_F_res = {};
All_F_damp = {};

for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list_a)
        continue
    end

    load(mat_tile_list{i_file});
    %% For 1 O.G
    S_1 = double(-timeSeriesData.Data(9,:));
    S_2 = double(-timeSeriesData.Data(10,:));
    S_3 = double(-timeSeriesData.Data(11,:));
    S_4 = double(-timeSeriesData.Data(12,:));

    %% For 2 O.G
    S_5 = double(-timeSeriesData.Data(3,:));
    S_6 = double(-timeSeriesData.Data(13,:));
    S_7 = double(-timeSeriesData.Data(14,:));
    S_8 = double(-timeSeriesData.Data(15,:));
    
    S_1 = filtfilt(b, a, S_1);
    S_2 = filtfilt(b, a, S_2);
    S_3 = filtfilt(b, a, S_3);
    S_4 = filtfilt(b, a, S_4);
    S_5 = filtfilt(b, a, S_5);
    S_6 = filtfilt(b, a, S_6);
    S_7 = filtfilt(b, a, S_7);
    S_8 = filtfilt(b, a, S_8);
    
    S_1 = Func_Resample(S_1,1024,128);
    S_2 = Func_Resample(S_2,1024,128);
    S_3 = Func_Resample(S_3,1024,128);
    S_4 = Func_Resample(S_4,1024,128);
    S_5 = Func_Resample(S_5,1024,128);
    S_6 = Func_Resample(S_6,1024,128);
    S_7 = Func_Resample(S_7,1024,128);
    S_8 = Func_Resample(S_8,1024,128);
    
    rz = [S_1;S_2;S_3;S_4;S_5;S_6;S_7;S_8];
    dt = 1/256;
    %dt = 1/2048;
    Ts = 2;
            
    [Nyy,N]= size(rz); %Nyy is the number of sensors and N the number of time step
    %
    [fn0,zeta0,phi,~] = SSICOV(rz,dt,'methodCOV',1,'Ts',Ts,'Nmin',2,'Nmax',30,'eps_cluster',0.2);
    i_file
    phi
    phi_list = zeros(16,1);
    if i_file <= 6
        phi_list(1:2:end, :) = phi(1,:);
        phi_list(2:2:end, :) = phi(2,:);
        
    elseif i_file >=7 && i_file <=9
        phi_list(1:2:end, :)  = phi(1,:);
        phi_list(2:2:end, :) = phi(3,:);
    else
        phi_list(1:2:end, :) = phi(1,:);
        phi_list(2:2:end, :) = phi(2,:);
    end
    phi_list
    All_F_res{i_file}= fn0;
    All_F_damp{i_file}= zeta0;
end

x_axis = [];
y_axis = [];
z_axis = [];
%list_a = [10,11];
k = 1;
for i_file = 1:13
    if ismember(i_file, list_a)
        continue
    end
    x_axis = [x_axis , All_F_res{i_file}];
    y_axis = [y_axis , All_F_damp{i_file}];
    z_axis = [z_axis , k*ones(1,length(All_F_res{i_file}))];
    k = k+1;
end

scatter3(x_axis,y_axis,z_axis,50,x_axis,'filled','o','MarkerEdgeColor','flat','MarkerFaceAlpha',0.4)
grid on
grid minor
xlim([0,high_freq]);
xlabel('$f_{r}$ (Hz)', 'Interpreter', 'latex');
ylabel('$\zeta $', 'Interpreter', 'latex');
zlabel('$i_{th}$ testing', 'Interpreter', 'latex');




%[h] = plotStabDiag(paraPlot.fn,rz(2,:),fs,paraPlot.status,paraPlot.Nmin,paraPlot.Nmax);
