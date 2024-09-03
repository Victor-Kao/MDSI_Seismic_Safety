clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
cutoff_freq = 40;
%[b,a] = Func_FilterDesign_highlow('low',cutoff_freq,4,1024);
[b,a] = Func_FilterDesign(4,cutoff_freq,4,1024);

i_file = 3;

load(mat_tile_list{i_file});
outputSignal = double(timeSeriesData.Data(9,:));
outputSignal = filtfilt(b, a, outputSignal);
t = timeSeriesData.Time;
inputSignal = double(timeSeriesData.Data(19,:));

res = Func_PSD_FRF_COH(inputSignal,outputSignal,[],[],[],1024);


% Compute the FFT of the input and output signals
N = length(t);
f = (0:N-1)*(fs/N); % Frequency vector

inputFFT = fft(inputSignal);
outputFFT = fft(outputSignal);

%inputFFT  = inputFFT(f<=cutoff_freq);
%outputFFT  = outputFFT(f<=cutoff_freq);



% Compute the frequency response (H = output/input)
H = outputFFT ./ inputFFT;
%H(1) = 0 + 1i*0;
H_1 = H(f<=cutoff_freq);
f_1 = f(f<=cutoff_freq);


%[fit,errdb]= rationalfit(f_1,H_1,"NPoles",250);

%[wn,zeta,p] = damp(fit);
%fit_original = fit.A;
%fit.A = fit.A([137,138,145:146]);
%fit.C = fit.C([137,138,145:146]);
%[resp,freq] = freqresp(fit,f_1);

num_dof = 4;
num_eigen = 100;
error_list = zeros(num_eigen ,1);
for i = 1:num_eigen
    [~,errdb]= rationalfit(f_1,H_1,"NPoles",2*i,Tolerance=-1);
    error_list(i) = errdb;
    disp(i)
end

[minValue, minIndex] = min(error_list);
[fit,errdb]= rationalfit(f_1,H_1,"NPoles",2* minIndex);
num_eigen = minIndex;
[resp,freq] = freqresp(fit,f_1);

fit_tune = copy(fit);
fit_tune_opt = copy(fit);
fit_tune_opt.A = [];
fit_tune_opt.C = [];
min_index_list = zeros(num_dof,1);

overdamped_list = find(abs(imag(fit.A))<=0.00001);

error_fit_dof = ones(num_eigen,1);
for j = 1:num_eigen
    if ismember(j,overdamped_list(2)/2)
        disp("jump")
        continue;
    end
    fit_tune.A = fit.A(2*j-1:2*j);
    fit_tune.C = fit.C(2*j-1:2*j);
    [test_resp,~] = freqresp(fit_tune,f_1);  
    error_fit_dof(j) = sum(abs(abs(resp)-abs((test_resp))));
end

% Sort the list and get the indices
[sortedData, sortedIndices] = sort(error_fit_dof);

% Select the top five smallest values and their indices
topFiveValues = sortedData(1:num_dof);
topFiveIndices = sortedIndices(1:num_dof);


for i_idx = 1:length(topFiveIndices)
    fit_tune_opt.A = [fit_tune_opt.A;fit.A(2*topFiveIndices(i_idx)-1:2*topFiveIndices(i_idx))];
    fit_tune_opt.C = [fit_tune_opt.C;fit.C(2*topFiveIndices(i_idx)-1:2*topFiveIndices(i_idx))];
end

%fit_tune_opt.A = [fit_tune_opt.A;fit.A(2*topFiveIndices(i_idx)-1:2*topFiveIndices(i_idx))];
%fit_tune_opt.C = [fit_tune_opt.C;fit.C(2*topFiveIndices(i_idx)-1:2*topFiveIndices(i_idx))];

[resp_opt,freq_opt] = freqresp(fit_tune_opt,f_1); 
plot(error_fit_dof(1:end-1))

info = [abs(real(fit_tune_opt.A)./abs(fit_tune_opt.A)),abs(fit_tune_opt.A)/(2*pi)]




% Only use the first half of the FFT results (positive frequencies)
H = H(1:N/2+1);
f = f(1:N/2+1);


% Define the order of the numerator and denominator of the transfer function
n = 100; % Order of the numerator
m = 100; % Order of the denominator

% Use invfreqz to fit a rational transfer function model to the frequency response data
%[b, a] = invfreqz(H, 2*pi*f/fs, n, m);

% Create a transfer function model
%sys = freqz(b, a, 4097);

% Display the transfer function
%disp(sys);

% Plot the Magnitude and Phase of the Frequency Response
figure;
subplot(2, 1, 1);
plot(f, abs(H));
hold on 
%plot(f, real(H));
%plot(f, imag(H));
%plot(f,abs(sys))
plot(res.f,abs(res.FRF))
plot(freq,abs(resp))
plot(freq_opt,abs(resp_opt))


title('Estimated Frequency Response Function (FRF) - Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('FRF','PSD','RFR','OPT')
grid on;
xlim([2,cutoff_freq])

subplot(2, 1, 2);
plot(f, angle(H));
title('Estimated Frequency Response Function (FRF) - Phase');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
xlim([2,cutoff_freq])
grid on;

%figure
%plot(error_list);





