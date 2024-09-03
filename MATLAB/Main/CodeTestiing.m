clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");

%% Testing the import used function
%dir_1_OG = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Mp2-1.OG\events";
%csv_files_1 = Func_FindCsvFiles(dir_1_OG);
%DT_Event_1G = Func_FindDateTime(csv_files_1);
%GM_1_all = Func_ImportMenhirData2Tab(csv_files_1{1});


%% Visualization of Menhir data for OMA
addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Event_sort\Event_1G2G_activated";
%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
fs = 1000;
Nmodes = 4;
list = [];
Table_SI = {};   
freq_SI = [];
lowCutoffFreq = 2; % Lower cutoff frequency (40 Hz)
highCutoffFreq = 100; % Upper cutoff frequency (60 Hz)
filterOrder = 4; % Filter order
figure 
for i_file = 1:length(mat_tile_list)
    if ismember(i_file, list)
        continue
    end
    load(mat_tile_list{i_file});
    az = data_T.Y_mm_s_GM_1_1G2G;
    %az = pwelch(az,[],[],[],1000);
    %az = diff(data_T.Z_mm_s)/(1/fs); % WARNING! Identify based on Velocity not ACC here 
    %az = bandpass(az,[17 100],1000);
    % Normalize the cutoff frequencies
    Wn = [lowCutoffFreq highCutoffFreq] / (fs / 2);
    % Design the Butterworth bandpass filter
    [b, a] = butter(filterOrder, Wn, 'bandpass');
    % Apply the bandpass filter to the signal
    az = filtfilt(b, a, az);
    [Saz,f] = pwelch(detrend(az),[],[],[],fs);
    az_f = Func_FFT_half(transpose(az),1000);
    plot(f,Saz/max(Saz))
    %plot(az_f.f,abs(az_f.s)/max(abs(az_f.s)))
    xlim([0,100])
    hold on
end


%% Testing the system identification and plot from Menhir data
%addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")
%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Event_sort\Event_1G2G_activated";
%%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Measure12032024\Event_sort\Event_1G2G_Hammer";
%mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
%fs = 1000;
%Nmodes = 4;
%list = [];
%Table_SI = {};   
%freq_SI = [];
%lowCutoffFreq = 2; % Lower cutoff frequency (40 Hz)
%highCutoffFreq = 22; % Upper cutoff frequency (60 Hz)
%filterOrder = 4; % Filter order
%for i_file = 1:length(mat_tile_list)
%    if ismember(i_file, list)
%        continue
%    end
%    load(mat_tile_list{i_file});
%    az = data_T.X_mm_s_GM_2_1G2G;
%    %az = pwelch(az,[],[],[],1000);
%    %az = diff(data_T.Z_mm_s)/(1/fs); % WARNING! Identify based on Velocity not ACC here 
%    %az = bandpass(az,[17 100],1000);
%    % Normalize the cutoff frequencies
%    Wn = [lowCutoffFreq highCutoffFreq] / (fs / 2);
%    % Design the Butterworth bandpass filter
%    [b, a] = butter(filterOrder, Wn, 'bandpass');
%    % Apply the bandpass filter to the signal
%    az = filtfilt(b, a, az);
%    [Saz,f] = pwelch(detrend(az),[],[],[],fs);         
%    Nperiod = 60; % The IRF is computed for 30 period on each mode
%    [fn,zeta] = modalID_singleSensor(az,Saz,f,Nmodes,fs,...
%        'plotOpt',0,'Nperiod',Nperiod,'PickingMethod','auto');
%    Table_SI{i_file,1} = fn;
%    Table_SI{i_file,2} = transpose(zeta);
%    freq_SI(i_file) = fn(1);
%end
%figure
%for j = 1:length(Table_SI)
%    scatter(Table_SI{j,1},Table_SI{j,2})
%    hold on 
%end


%% Testing the system identification and plot from Accereration data
%addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")
%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
%mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
%fs = 1024;
%Nmodes = 1;
%list = [];
%Table_SI = {}; 
%freq_SI = [];
%lowCutoffFreq = 4; % Lower cutoff frequency (40 Hz)
%highCutoffFreq = 100; % Upper cutoff frequency (60 Hz)
%filterOrder = 6; % Filter order
%for i_file = 1:length(mat_tile_list)
%    if ismember(i_file, list)
%        continue
%    end
%    load(mat_tile_list{i_file});
%    az = double(timeSeriesData.Data(19,:));
%    %az = diff(data_T.Z_mm_s_GM_1_1G2G)/(1/fs); % WARNING! Identify based on Velocity not ACC here 
%    %az = bandpass(az,[17 100],1000);
%    % Normalize the cutoff frequencies
%    Wn = [lowCutoffFreq highCutoffFreq] / (fs / 2);
%    % Design the Butterworth bandpass filter
%    [b, a] = butter(filterOrder, Wn, 'bandpass');
%    % Apply the bandpass filter to the signal
%    az = filtfilt(b, a, az);
%    [Saz,f] = pwelch(detrend(az),[],[],[],fs);         
%    Nperiod = 10; % The IRF is computed for 30 period on each mode
%    [fn,zeta] = modalID_singleSensor(az,Saz,f,Nmodes,fs,...
%        'plotOpt',1,'Nperiod',Nperiod,'PickingMethod','auto');
%    Table_SI{i_file,1} = fn;
%    Table_SI{i_file,2} = transpose(zeta);
%    freq_SI(i_file) = fn(1);
%end
%figure
%for j = 1:length(Table_SI)
%    scatter(Table_SI{j,1}(1),Table_SI{j,2}(1))
%    hold on 
%end


%% Testing date time for Acc
%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
%mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
%for i_file = 1:length(mat_tile_list )
%    load(mat_tile_list{i_file});
%    timeSeriesData.Parameters.MeasurementBegin
%end


%addpath("D:\MDSI_project\MATLAB\Lib\ECheynet-modalID_singleSensor-b67ee9f")
%dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
%mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
%
%for i = 1:length(mat_tile_list)
%    load(mat_tile_list{i})
%
%    % Use fileparts to split the file path
%    [~, fileName, fileExt] = fileparts(mat_tile_list{i});
%    
%    % Concatenate the filename and extension
%    fullFileName = [fileName];
%    
%    % Display the result
%    disp(fullFileName);
%    disp(timeSeriesData.Parameters.MeasurementBegin)
%end
