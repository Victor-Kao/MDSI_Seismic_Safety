function [b, a] = Func_FilterDesign(lowCutoffFreq,highCutoffFreq,filterOrder,fs)
    % Normalize the cutoff frequencies
    Wn = [lowCutoffFreq highCutoffFreq] / (fs / 2);
    % Design the Butterworth bandpass filter
    [b, a] = butter(filterOrder, Wn, 'bandpass');
end