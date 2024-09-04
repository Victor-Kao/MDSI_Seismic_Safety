clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
load("D:/MDSI_project/DATA_GM_RawData/CellArrayAllEvent.mat");

%%
plot_ = 1;

%% 
load("Menhir_result.mat");


cut_f = 100;
% Apply DBSCAN
f  = linspace(0,500,2049);
All_s = All_s(:,f<=cut_f);
norm_All_nm = normalize(All_s,2,"range");
All_s_smooth = movmean(All_s,5,2);
norm_All_s = normalize(All_s_smooth,2,"range");
f_s = f(f<=cut_f);

% Step 1: Compute pairwise distances between signals (1500 signals, each of length 3000)
distances = pdist(norm_All_s, 'euclidean'); 

% Step 2: Perform hierarchical clustering using the linkage method
Z = linkage(distances, 'ward');  % 'ward' method is used here, but other methods can be tried


% Step 3: Plot the dendrogram to visualize the hierarchical clustering
figure;
dendrogram(Z,0);
title('Hierarchical Binary Cluster Tree', 'Interpreter', 'latex');
set(gca, 'XTickLabel', []);
xlabel('Leaf Nodes,~Index of Records (not shown)', 'Interpreter', 'latex');
ylabel('Distance', 'Interpreter', 'latex')

% Step 4: Create clusters from the linkage result
numClusters = 7;  % Define the number of clusters you want, adjust as needed
Distance = 14;
%clusters = cluster(Z, 'maxclust', numClusters);
clusters = cluster(Z, 'cutoff',Distance,'criterion','distance');
h = yline(Distance,'r--',['Threshold = ',num2str(Distance)]);
h.FontSize = 8;  % Set font size to 14



validIndices = true(size(clusters));
numClusters_sec = numClusters;


% Loop through each cluster and check its size
for i = 1:numClusters
    % Find the signals belonging to cluster 'i'
    clusterIndices = find(clusters == i);
    
    % Check if the number of signals in this cluster is smaller than 5
    if length(clusterIndices) < 1
        % Mark these signals for removal
        validIndices(clusterIndices) = false;
        numClusters_sec = numClusters_sec -1;
        
    end
end

% Remove the signals belonging to small clusters
norm_All_s = norm_All_s(validIndices, :);
norm_All_nm = norm_All_nm(validIndices, :);
All_s = All_s(validIndices, :);
clusters = clusters(validIndices);

numClusters = max(clusters);  % Number of clusters

if plot_ == 1
    for i = 1:numClusters
        % Find the indices of the signals that belong to cluster 'i'
        clusterIndices = find(clusters == i);
        
        if  isempty(clusterIndices)
            continue;
        end
        
        disp(i)
        
        % Create a new figure for this cluster
        figure;
        sgtitle(['Cluster ' num2str(i) ', $N =$' num2str(length(clusterIndices))], 'Interpreter', 'latex');
        
        % Plot each signal in this cluster
        subplot(2,1,1)
        
        for j = 1:length(clusterIndices)
            energy = trapz(All_s(clusterIndices(j), :));
            plot(f_s, All_s(clusterIndices(j), :)/energy);
            hold on;
        end
        title('Normalized PSD by Energy $\tilde{S}_{x}(f)$', 'Interpreter', 'latex')
        xlabel('freq (Hz)', 'Interpreter', 'latex');
        ylabel('$\tilde{S}_{x}(f)$', 'Interpreter', 'latex');
    
        subplot(2,1,2)
        %norm_All_nm = normalize(All_s(clusterIndices, :),2,"range");
        std_dev =  std(norm_All_nm(clusterIndices, :),0,1);
        plot(f_s, mean(norm_All_nm(clusterIndices, :),1),'k');
        hold on;
        plot(f_s ,mean(norm_All_nm(clusterIndices, :),1)+1*std_dev,'c');
        plot(f_s, mean(norm_All_nm(clusterIndices, :),1)+3*std_dev,'b');

        for j = 1:length(clusterIndices)
            plot(f_s, norm_All_s(clusterIndices(j), :),'r:');
        end

        plot(f_s, mean(norm_All_nm(clusterIndices, :),1),'k');
        plot(f_s ,mean(norm_All_nm(clusterIndices, :),1)+1*std_dev,'c');
        plot(f_s, mean(norm_All_nm(clusterIndices, :),1)+3*std_dev,'b');

        legend({'$\mathbf{E}[\hat{S}_{x}(f)]$','$\mathbf{E}[\hat{S}_{x}(f)] + 1\sigma$','$\mathbf{E}[\hat{S}_{x}(f)]+ 3\sigma$','$\hat{S}_{x}(f)$'}, 'Interpreter', 'latex');
        title('Normalized PSD by Maximum Amplitude $\hat{S}_{x}(f)$', 'Interpreter', 'latex')
        xlabel('freq (Hz)', 'Interpreter', 'latex');
        ylabel('$\hat{S}_{x}(f)$', 'Interpreter', 'latex');


        stddev_all(i,1) = sum(std_dev)/length(std_dev);

        hold off;
    end
end

