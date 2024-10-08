function [Y] = filter_shm (X,filterCoef)
% Feature Extraction: Filters signals with FIR filter
%
% VERBOSE FUNCTION CALL:
%   [Filtered Signal Matrix] = Filter Signal (Signal Matrix, Filter
%   Coefficients)
%
% CALL METHODS:
%   [Y] = filter_shm (X, filterCoef)
%
% DESCRIPTION:
%   Filters a signal matrix using fast fft convolution with specified
%   filter impulse response coefficients. The fft of the signal and the
%   filter impulse response are zero padded to be of adequate length and
%   multiplied in the frequency domain. The inverse Fourier transform is
%   then taken. After the signal has been convolved by the filter impulse
%   response, the end is cropped where the impulse response begins to run
%   off the length of the signal such that the signal length of the input
%   signal is equal to the signal length of the output signal.
%
% OUTPUTS:
%   Y (SAMPLES, CHANNELS, INSTANCES) : filtered signal matrix
%
% INPUTS:
%   X (SAMPLES, CHANNELS, INSTANCES) : signal matrix
%
%   filterCoef (NFILTER, 1): filter tap coefficients
%
% SEE ALSO:
%   fir1_shm

% AUTHORS:
%   Luke Robinson
%
%   LA-CC-14-046
%   Copyright (c) 2010, Triad National Security, LLC 
%   All rights reserved.
%
% LAST MODIFIED:
%   7/16/2013
%
% REFFERENCES:
%   [1] Smith, S., The Scientist and Engineer's Guide to Digital Signal
%   Processing. California Technical Pub. http://www.dspguide.com/ch9/3.htm, Web.

%Get signal space dimensions
[nSignal,nChannel,nInstance] = size(X);

%Get fft length for fft filtering opperation
filterNFFT = nSignal+length(filterCoef)-1;

%Generate a filterCoeff space for fft filtering of signal space.
filterCoef = repmat(filterCoef,[1,nChannel,nInstance]);

%Filter signal space using fft filtering method
Y = ifft((fft(filterCoef,filterNFFT,1).*fft(X,filterNFFT,1)),[],1);

%Resize signal space to not include filter runoff.
Y = Y(1:nSignal,:,:);

end