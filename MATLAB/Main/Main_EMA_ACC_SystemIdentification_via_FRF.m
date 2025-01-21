clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\pickpeaks");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain_update";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [];
fs = 1024;
noverlap = [];
window = [];
f = 5000;
low_freq = 4;
high_freq = 25;
[b,a] = Func_FilterDesign(low_freq,high_freq,4,1024);

% Specify index of hammer test events
i_file =1;
i_pos = 9;

load(mat_tile_list{i_file});
outputSignal = double(timeSeriesData.Data(i_pos,:));
%outputSignal = filtfilt(b, a, outputSignal);
inputSignal = double(timeSeriesData.Data(19,:));

% Compute the Power Spectral Density (PSD) of the input signal
[pxxInput, ~] = pwelch(inputSignal, window, noverlap,f, fs);
% Compute the Cross-Power Spectral Density (CSD) between input and output
[pxy, ~] = cpsd(inputSignal, outputSignal, window, noverlap, f, fs);
% Compute the Frequency Response Function (FRF)
FRF = pxy ./ pxxInput;
% Compute the PSD of the output signal for comparison
[pxxOutput, ~] = pwelch(outputSignal, window, noverlap, f, fs);
% Compute the coherence between input and output signals
[Cxy, freq] = mscohere(inputSignal, outputSignal, window, noverlap, f, fs);

FRF = FRF(freq<=high_freq);
freq = freq(freq<=high_freq);


n = 60; % Order of numerator
m = n; % Order of denominator
% Angular frequency (rad/s)
omega = 2 * pi * freq;
% Fit the FRF data
[b, a] = invfreqs(FRF, omega, n, m);



% Evaluate the fitted FRF
fitted_FRF = freqs(b, a, omega);

% Plot the original vs. fitted FRF
figure;
subplot(2,1,1);
plot(freq, abs(FRF), 'b', 'DisplayName', 'Original');
hold on;
plot(freq, abs(fitted_FRF), 'r--', 'DisplayName', 'Fitted');
xlabel('Frequency (Hz)');
ylabel('|H(j\omega)|');
legend;
grid on;

subplot(2,1,2);
plot(freq, angle(FRF), 'b', 'DisplayName', 'Original');
hold on;
plot(freq, angle(fitted_FRF), 'r--', 'DisplayName', 'Fitted');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
legend;
grid on;
