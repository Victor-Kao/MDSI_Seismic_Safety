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
    az = data_T.Z_mm_s_GM_1_1G2G;
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