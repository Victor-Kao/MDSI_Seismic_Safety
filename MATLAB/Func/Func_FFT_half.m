function freq_signal = Func_FFT_half(time_singal, nfft, fs)
    if ~isequal(size(time_singal), [1,length(time_singal)])
        time_singal = transpose(time_singal);
    end

    if nfft > length(time_singal)
        time_singal = [time_singal, zeros(1, nfft - length(time_singal))];
    elseif isempty(nfft)
    else
        disp(['Warning! nfft: ', num2str(nfft), ' <= size of input, x: ',num2str(length(time_singal))]);
        disp(['Using default size, x = ',num2str(length(time_singal))]);
    end

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