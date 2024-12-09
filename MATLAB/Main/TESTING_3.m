clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
addpath("D:\MDSI_project\MATLAB\Lib\pickpeaks");


%% Acc Data for 2 O.G. (Channel 15), first mode
dir_activate_1G2G = "D:\MDSI_project\DATA_GM_RawData\DATA_ACC_Measure12032024\DATA_Hammer\Time_domain";
mat_tile_list = Func_FindMatFiles(dir_activate_1G2G);
list_a = [10,11];
fs = 1024;
cutoff_freq = 40;
%[b,a] = Func_FilterDesign_highlow('low',cutoff_freq,4,1024);
[b,a] = Func_FilterDesign(4,cutoff_freq,4,1024);

i_file = 15;

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


num_dof = 4;
num_eigen = 50;
error_list = zeros(num_eigen ,1);
for i = num_eigen:num_eigen
    [~,errdb]= rationalfit(f_1,H_1,"NPoles",2*i,Tolerance=-1);
    error_list(i) = errdb;
    disp(i)
end

[minValue, minIndex] = min(error_list);
[fit,errdb]= rationalfit(f_1,H_1,"NPoles",2* minIndex);
num_eigen = minIndex;
[resp,freq] = freqresp(fit,f_1);




%fit.A = fit.A([19,20,59,60,87,88]);
%fit.C = fit.C([19,20,59,60,87,88]);

info = [abs(real(fit.A)./abs(fit.A)),abs(fit.A)/(2*pi)];
%[resp_,freq_] = freqresp(fit,f_1);


[peaks,criterion] = pickpeaks(abs(resp),2,1);

ranges = zeros(length(peaks),2);
peak_detected = freq(peaks);

for i_range = 1:length(peaks)
    ranges(i_range,:) = [peak_detected(i_range)-0.5,peak_detected(i_range)+0.5];
    
end

% Initialize a logical array to store whether each element is in any range
in_range = false(size(info(:,2)));

% Loop through each range and update the logical array

Copy_fit = copy(fit);
Original_fit = copy(fit);
ROM_cand_fit = copy(fit);
ROM_fit = copy(fit);
ROM_fit.A = [];
ROM_fit.C = [];
for i_r = 1:size(ranges, 1)
    in_range = (info(:,2) >= ranges(i_r, 1) & info(:,2) <= ranges(i_r, 2));
    % Extract the data that falls within the specified ranges
    Copy_fit.A = fit.A(in_range);
    Copy_fit.C = fit.C(in_range);

    %Copy_fit.C(imag(Copy_fit.A) < 0.0005) = [];
    %Copy_fit.A(imag(Copy_fit.A) < 0.0005) = [];
    
    disp(['Fitting error for mode ',num2str(i_r)]);
    error_list  = ones(length(Copy_fit.A),1);
    for j = 1:length(Copy_fit.A)/2
        ROM_cand_fit.A = Copy_fit.A([2*j-1,2*j]);
        ROM_cand_fit.C = Copy_fit.C([2*j-1,2*j]);
        
        [FRF_cand,freq_cand] = freqresp(ROM_cand_fit,f_1);
        error = sum(abs(resp) - abs(FRF_cand));
        error_list(j) = error;
    end

    [~,min_err_idx] = min(error_list);
    ROM_fit.A = [ROM_fit.A; Copy_fit.A([2*min_err_idx-1,2*min_err_idx])];
    ROM_fit.C = [ROM_fit.C; Copy_fit.C([2*min_err_idx-1,2*min_err_idx])];
end

[FRF_ROM,freq_ROM] = freqresp(ROM_fit,f_1);

info = [abs(real(ROM_fit.A)./abs(ROM_fit.A)),abs(ROM_fit.A)/(2*pi)]






%[freq_(peaks)-3 ; freq_(peaks)+3]


figure
plot(f_1, abs(H_1));
hold on 
plot(freq,abs(resp));
plot(freq_ROM, abs(FRF_ROM));