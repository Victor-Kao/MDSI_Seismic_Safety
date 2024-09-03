function result = Func_PSD_FRF_COH_freq(freq_output, inputSignal, outputSignal)
    %% Transfer menhir data (outputSignal) from velocity to acceleration
    omega = 2 * pi * freq_output; % Angular frequency
    % Convert velocity to acceleration in the frequency domain
    outputSignal = (1j * omega) .* outputSignal;

    %% Compute PDS FRF COH
    pxxInput = abs(inputSignal).^2;
    
    % Compute the Cross-Power Spectral Density (CSD) between input and output
    pxy = inputSignal.* conj(outputSignal);
    
    % Compute the Frequency Response Function (FRF)
    FRF = pxy ./ pxxInput;
    FRF(1) = 0 + 1i*0;

    % Compute the PSD of the output signal for comparison
    pxxOutput = abs(outputSignal).^2;
    
    % Compute the coherence between input and output signals
    Cxy = (abs(pxy).^2) ./ (pxxInput .* pxxOutput);
    
    result = table(freq_output,pxxInput,pxxOutput,pxy,FRF,Cxy,'VariableNames', {'f','Pxx_in','Pxx_out','Pxy','FRF','Cxy'});

end