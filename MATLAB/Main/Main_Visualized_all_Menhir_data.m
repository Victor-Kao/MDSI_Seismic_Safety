clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");

load("D:/MDSI_project/DATA_GM_RawData/CellArrayAllEvent.mat");

figure
j = 1;
for i = 1:length(All_event)
    
    fs = ceil(1/(All_event{i}.time(3)-All_event{i}.time(2)));
    [pxx,f] = pwelch(All_event{i}.Z_mm_s,[],[],[],fs);
    %signal_FFT = Func_FFT_half(All_event{i}.Z_mm_s,1000);
    energy = trapz(f, pxx);
    

    %plot(signal_FFT.f,abs(signal_FFT.s))
    

    if length(f) == 2049 
        energy_store(j,1) = j;
        energy_store(j,2) = energy;
        All_s(j,:) = pxx;
        Freq = f;
        j = j+1;
    end
    plot(f,pxx,'r');
    hold on 
end





%plot(Freq,mean( All_s,1),'b');

%figure
%scatter(energy_store(:,1),log(energy_store(:,2)))




%figure
%for i = 1:length(All_event)
%    signal_FFT = Func_FFT_half(All_event{i}.Z_mm_s,1000);
%    
%    plot(signal_FFT.f,imag(signal_FFT.s))
%    hold on 
%end



