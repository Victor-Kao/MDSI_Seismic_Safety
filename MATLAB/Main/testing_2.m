clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
cutoff_freq = 100;
[b,a] = Func_FilterDesign_highlow('high',4,4,1024);

i_file = 1;

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
H(1) = 0 + 1i*0;
H_1 = H(f<=cutoff_freq);
f_1 = f(f<=cutoff_freq);

[fit,errdb]= rationalfit(f_1,H_1,"NPoles",100,Tolerance=-1);

fit_original = copy(fit);
[resp_1,~] = freqresp(fit_original,f_1);

error_list = ones(50,1);
for j = 1:50
    fit.A = fit_original.A(2*j-1:2*j);
    fit.C = fit_original.C(2*j-1:2*j);
    [resp,~] = freqresp(fit,f_1);
    error_list(j) = sum(abs(abs(resp_1)-abs((resp))));
end
%list = [1,2,19,20,94,93,82,81,86,85];
list = [50,49,94,93,64,63,56,55,2,1,6,5];
fit.A = fit_original.A(list);
fit.C = fit_original.C(list);


%[resp_1,~] = freqresp(fit_original,f_1);
[resp,~] = freqresp(fit,f_1);
error = sum(abs(abs(resp_1)-abs((resp))))

H_1 = H_1./ ((2 * pi * f_1).^2);
plot(f_1,imag(H_1))
hold on

plot(res.f,imag(res.FRF))
%plot(f_1,abs(resp))

xlim([4,100])






