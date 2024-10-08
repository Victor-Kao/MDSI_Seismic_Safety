clear;
clc;
close all;

addpath('D:\MDSI_project\MATLAB\Lib\pickpeaks');
load('cluster_1G_60hz_15_7c.mat')
load("Menhir_result_1G.mat");

histogram(hc);

cut_f = 40;
f = linspace(0,500,size(All_s,2));
All_s = All_s(:,f<=cut_f);
All_s = All_s./trapz(All_s,2);
f_s = f(f<=cut_f);

all_peaks = {};

k = 1;
for p = 1:length(unique(hc)) 
    potential_peak = zeros(10000,1);
    p_peak = zeros(10000,1);
    index = find(hc==p);
    for i = 1:length(index)
        peaks = pickpeaks(All_s(index(i),:),0.2,0);
        %plot(f_s,All_s(index(i),:))
        %hold on 
        for j = 1:length(peaks) 
            scatter(f_s(peaks(j)),All_s(index(i),peaks(j)),'ro');
            potential_peak(k,1) = f_s(peaks(j));
            k = k+1;
        end 
    end
    potential_peak(potential_peak == 0) = [];
    all_peaks{p} = potential_peak;
end


figure
for f = 1:length(unique(hc))
    x = f*ones(length(all_peaks{f}),1);
    scatter(all_peaks{f},f,'red','filled','o','MarkerFaceAlpha',0.3)
    hold on
end
grid on 
xlim([0,40])
ylim([0,length(unique(hc))+1])
