function [FRF_ROM,freq_ROM,SI_info] = Func_Proposed_SI_FRF(inputSignal,outputSignal,t,fs,cutoff_freq,NumDof,SearchRange)


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


%num_dof = 4;
num_eigen = 50;
error_list = zeros(num_eigen ,1);
for i = num_eigen:num_eigen
    [~,errdb]= rationalfit(f_1,H_1,"NPoles",2*i,Tolerance=-1);
    error_list(i) = errdb;
end

[~, minIndex] = min(error_list);
[fit,~]= rationalfit(f_1,H_1,"NPoles",2* minIndex);
%num_eigen = minIndex;
[resp,freq] = freqresp(fit,f_1);


%fit.A = fit.A([19,20,59,60,87,88]);
%fit.C = fit.C([19,20,59,60,87,88]);

info = [abs(real(fit.A)./abs(fit.A)),abs(fit.A)/(2*pi)];
%[resp_,freq_] = freqresp(fit,f_1);

[peaks,~] = pickpeaks(abs(resp),NumDof,0);

ranges = zeros(length(peaks),2);
peak_detected = freq(peaks);

for i_range = 1:length(peaks)
    ranges(i_range,:) = [peak_detected(i_range)-SearchRange,peak_detected(i_range)+SearchRange];
    
end

% Initialize a logical array to store whether each element is in any range
%in_range = false(size(info(:,2)));

% Loop through each range and update the logical array

Copy_fit = copy(fit);
%Original_fit = copy(fit);
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
        
        [FRF_cand,~] = freqresp(ROM_cand_fit,f_1);
        error = sum(abs(resp) - abs(FRF_cand));
        error_list(j) = error;
    end

    [~,min_err_idx] = min(error_list);
    ROM_fit.A = [ROM_fit.A; Copy_fit.A([2*min_err_idx-1,2*min_err_idx])];
    ROM_fit.C = [ROM_fit.C; Copy_fit.C([2*min_err_idx-1,2*min_err_idx])];
end

[FRF_ROM,freq_ROM] = freqresp(ROM_fit,f_1);

SI_info = [abs(real(ROM_fit.A)./abs(ROM_fit.A)),abs(ROM_fit.A)/(2*pi)];





end