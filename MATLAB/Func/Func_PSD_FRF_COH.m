function result = Func_PSD_FRF_COH(inputSignal,outputSignal,window,noverlap,f,fs)

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

    result = table(freq,pxxInput,pxxOutput,pxy,FRF,Cxy,'VariableNames', {'f','Pxx_in','Pxx_out','Pxy','FRF','Cxy'});

end