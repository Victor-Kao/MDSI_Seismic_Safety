function freq_signal = Func_FFT_half(time_singal, fs)
    L = length(time_singal);
    Y = fft(time_singal);
    P2 = Y/L;
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;

    if ~isequal(size(P1), [1,length(P1)])
        P1 = transpose(P1);
    end

    if ~isequal(size(P1), size(f))
        f = transpose(f);
    end

    freq_signal = table(f, P1, 'VariableNames', {'f', 's'});
end