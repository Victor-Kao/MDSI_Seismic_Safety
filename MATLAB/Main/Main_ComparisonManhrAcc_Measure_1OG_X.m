clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")

%% Intialized
Nmodes = 2;
list_a = [1,2,7,10,11,13,14,15];
%list_a = [];
list_m = [];
Table_SI_Menhir = {};   
Table_SI_Acc = {}; 
freq_SI_Menhir = [];
freq_SI_Acc = [];

%% butter-filter = [Order = 6,Low = 15, High = 100, Fs = 1000]
lowCutoffFreq = 10; % Lower cutoff frequency (40 Hz)
highCutoffFreq = 50; % Upper cutoff frequency (60 Hz)
filterOrder = 6; % Filter order
% Normalize the cutoff frequencies
Wn_m = [lowCutoffFreq highCutoffFreq] / (1000 / 2);
% Design the Butterworth bandpass filter
[b_m, a_m] = butter(filterOrder, Wn_m, 'bandpass');


%% butter-filter = [Order = 6,Low = 4, High = 100, Fs = 1024]
lowCutoffFreq = 4; % Lower cutoff frequency (40 Hz)
highCutoffFreq = 45; % Upper cutoff frequency (60 Hz)
filterOrder = 6; % Filter order
% Normalize the cutoff frequencies
Wn_a = [lowCutoffFreq highCutoffFreq] / (1024 / 2);
% Design the Butterworth bandpass filter
[b_a, a_a] = butter(filterOrder, Wn_m, 'bandpass');


%% Menhir Data for 1 O.G. first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer_validation_1OG";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
fs = 1000;

for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list_m)
        continue
    end
    load(mat_tile_list{i_file});
    %az = diff(data_T.Z_mm_s)/(1/fs); % WARNING! Identify based on Velocity not ACC here 
    az = data_T.X_mm_s;
    % Apply the bandpass filter to the signal
    az = filtfilt(b_m, a_m, az);
    [Saz,f] = pwelch(detrend(az),[],[],[],fs);         
    Nperiod = 60; % The IRF is computed for 30 period on each mode
    [fn,zeta] = modalID_singleSensor(az,Saz,f,Nmodes,fs,...
        'plotOpt',0,'Nperiod',Nperiod,'PickingMethod','auto');
    Table_SI_Menhir{i_file,1} = fn;
    Table_SI_Menhir{i_file,2} = transpose(zeta);
    freq_SI_Menhir(i_file) = fn(1);
end


%% Acc Data for 1 O.G. (Channel 7), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
fs = 1024;
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list_a)
        continue
    end
    load(mat_tile_list{i_file});
    az = double(timeSeriesData.Data(7,:));
    % Normalize the cutoff frequencies
    Wn = [lowCutoffFreq highCutoffFreq] / (fs / 2);
    % Design the Butterworth bandpass filter
    [b, a] = butter(filterOrder, Wn, 'bandpass');
    % Apply the bandpass filter to the signal
    az = filtfilt(b, a, az);
    [Saz,f] = pwelch(detrend(az),[],[],[],fs);         
    Nperiod = 60; % The IRF is computed for 30 period on each mode
    [fn,zeta] = modalID_singleSensor(az,Saz,f,Nmodes,fs,...
        'plotOpt',0,'Nperiod',Nperiod,'PickingMethod','auto');
    Table_SI_Acc{i_file,1} = fn;
    Table_SI_Acc{i_file,2} = transpose(zeta);
    freq_SI_acc(i_file) = fn(1);
end


figure
for j = 1:length(Table_SI_Menhir)
    scatter(Table_SI_Menhir{j,1}(1:2),Table_SI_Menhir{j,2}(1:2),"MarkerEdgeColor","#0072BD")
    hold on 
end

for k = 1:length(Table_SI_Acc)
    if ismember(k, list_a)
        continue
    end
    scatter(Table_SI_Acc{k,1}(1:2),Table_SI_Acc{k,2}(1:2),"MarkerEdgeColor","#D95319")
    hold on 
end