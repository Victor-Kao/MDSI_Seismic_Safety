function [FRF, freq] = Func_buildFRF(fn_list, damp_list, cutoff_freq_high, fs)
    % Define the frequency range for the FRF    
    f = 0:(1/fs):cutoff_freq_high;    
    omega = 2 * pi * f; % Convert frequency to rad/s
    FRF = zeros(1,length(f)) + 1i*zeros(1,length(f));
    for i = 1:length(fn_list)
        omega_n = 2*pi*fn_list(i);
        zeta = damp_list(i);
        H1 = 1 ./ (omega_n^2 - omega.^2 + 2j*zeta*omega_n*omega);
        FRF = FRF + H1;
    end
    freq = f;
end