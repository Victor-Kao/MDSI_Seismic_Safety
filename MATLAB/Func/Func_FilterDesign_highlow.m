function [b, a] = Func_FilterDesign_highlow(ftype,CutoffFreq,filterOrder,fs)
    % Normalize the cutoff frequencies
    Wn = CutoffFreq / (fs / 2);
    % Design the Butterworth bandpass filter
    [b, a] = butter(filterOrder, Wn, ftype);
end

