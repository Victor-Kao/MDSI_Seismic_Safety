clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
load("D:/MDSI_project/DATA_GM_RawData/CellArrayAllEvent.mat");

%[b,a] = Func_FilterDesign_highlow('low',100,4,1000);
%
%j = 1;
%for i = 1:length(All_event)
%    fs = ceil(1/(All_event{i}.time(3)-All_event{i}.time(2)));
%    outputSignal = filtfilt(b, a, All_event{i}.Z_mm_s);
%    [pxx,f] = pwelch(outputSignal,[],[],[],fs);
%    energy = trapz(f, pxx);
%    if length(f) == 2049 
%        All_s(j,:) = pxx;
%        Freq = f;
%        j = j+1;
%    end
%end

%% 
load("Menhir_result.mat");

% Set the parameters for DBSCAN
epsilon = 0.75; % Adjust this value based on your data
minPts = 2;    % Adjust this value based on your data

cut_f = 150;
% Apply DBSCAN
f  = linspace(0,500,2049);
All_s = All_s(:,f<=cut_f);
norm_All_s = normalize(All_s,2,"range");

labels = dbscan(norm_All_s, epsilon, minPts, 'Distance', 'euclidean');

f = f(f<=cut_f);

unique_labels = unique(labels);
num_clusters = length(unique_labels);

%figure
%for j = 1:length(norm_All_s)
%    plot(f,norm_All_s(j,:));
%    hold on
%end



for i = 1:num_clusters
    cluster_label = unique_labels(i);
    cluster_signals = All_s(labels == cluster_label, :);

    if cluster_label == -100
        labels_mo = dbscan(cluster_signals, 0.0001, 2, 'Distance', 'euclidean');
        uni_mo = unique(labels_mo);
        for m = 1:length(uni_mo)
            cluster_signals_mo = cluster_signals(labels_mo== uni_mo(m), :);
            figure;
            plot(f,cluster_signals_mo);
            title(['Cluster from -1: ', num2str(uni_mo(m))]);
            xlabel('Time');
            ylabel('Amplitude');
        end

    elseif cluster_label == 100
        labels_o = dbscan(cluster_signals, 0.0001, 2, 'Distance', 'euclidean');
        uni_o = unique(labels_o);
        for n = 1:length(uni_o)
            cluster_signals_o = cluster_signals(labels_o== uni_o(n), :);
            figure;
            plot(f,cluster_signals_o);
            title(['Cluster from 1: ', num2str(uni_o(n))]);
            xlabel('Time');
            ylabel('Amplitude');
        end
    else
        figure;
        plot(f,cluster_signals);
        title(['Cluster ', num2str(cluster_label)]);
        xlabel('Time');
        ylabel('Amplitude');
    end
end

figure
histogram(labels)

figure
scatter(1:length(labels),labels,'.')
