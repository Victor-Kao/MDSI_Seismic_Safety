clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");

%% For 1G Menhir data
load("D:/MDSI_project/MATLAB/Main/CellArrayEvent_EG.mat");

%% For 2G Menhir data
%load("D:/MDSI_project/DATA_GM_RawData/CellArrayEvent_2G.mat");

figure
j = 1;
all_size = zeros(size(All_event));
for k = 1:length(All_event)
    all_size(k) = size(All_event{k},1);
end
%nfft = ceil(max(all_size/1000))*1000;
nfft = 5000;


for i = 1:length(All_event)
    
    fs = ceil(1/(All_event{i}.time(3)-All_event{i}.time(2)));
    [pxx,f] = pwelch(All_event{i}.Z_mm_s,[],[],nfft,fs);
    %signal_FFT = Func_FFT_half(All_event{i}.Z_mm_s,nfft,1000);
    energy = trapz(f, pxx);
    
    %length(signal_FFT.f)
    %plot(signal_FFT.f,abs(signal_FFT.s))
    %hold on

    if length(f) ~= 2000      
        energy_store(j,1) = j;
        energy_store(j,2) = energy;
        All_s(j,:) = pxx;
        Freq = f;
        j = j+1;
        plot(f,pxx/energy);
        hold on
    end
    grid on 
    xlim([0,30])

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



