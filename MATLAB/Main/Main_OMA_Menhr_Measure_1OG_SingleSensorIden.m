clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Event_sort\Event_all_activated";

mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
fs = 1000;
Nmodes = 2;

lowCutoffFreq = 4; % Lower cutoff frequency (40 Hz)
highCutoffFreq = 35; % Upper cutoff frequency (60 Hz)
filterOrder = 4; % Filter order

%% Z case 
list = [];
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list)
        continue
    end
    load(mat_tile_list{i_file});
    az = data_T.Z_mm_s_GM_1_all;
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
    Table_SI_Z{i_file,1} = fn;
    Table_SI_Z{i_file,2} = transpose(zeta);
    frequ_SI_Z(i_file) = fn(1);
end

%% X case 
list = [];
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list)
        continue
    end
    load(mat_tile_list{i_file});
    az = data_T.X_mm_s_GM_1_all;
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
    Table_SI_X{i_file,1} = fn;
    Table_SI_X{i_file,2} = transpose(zeta);
    frequ_SI_X(i_file) = fn(1);
end

%% Y case 
list = [40];
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list)
        continue
    end
    load(mat_tile_list{i_file});
    az = data_T.Y_mm_s_GM_1_all;
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
    Table_SI_Y{i_file,1} = fn;
    Table_SI_Y{i_file,2} = transpose(zeta);
    frequ_SI_Y(i_file) = fn(1);
end


figure
for j = 1:length(Table_SI_X)
    subplot(3,1,1)
    scatter(Table_SI_X{j,1},Table_SI_X{j,2})
    hold on 
    xlim([0,35])
    grid on
    grid minor
    xlabel('$f_{r}$ (Hz)', 'Interpreter', 'latex');
    ylabel('$\zeta $', 'Interpreter', 'latex');
    title('System Identification during operation, x dir', 'Interpreter', 'latex');

    subplot(3,1,2)
    scatter(Table_SI_Y{j,1},Table_SI_Y{j,2})
    hold on 
    xlim([0,35])
    grid on
    grid minor
    xlabel('$f_{r}$ (Hz)', 'Interpreter', 'latex');
    ylabel('$\zeta $', 'Interpreter', 'latex');
    title('System Identification during operation, y dir', 'Interpreter', 'latex');

    subplot(3,1,3)
    scatter(Table_SI_Z{j,1},Table_SI_Z{j,2})
    hold on 
    xlim([0,35])
    grid on
    grid minor
    xlabel('$f_{r}$ (Hz)', 'Interpreter', 'latex');
    ylabel('$\zeta $', 'Interpreter', 'latex');
    title('System Identification during operation, z dir', 'Interpreter', 'latex');
end