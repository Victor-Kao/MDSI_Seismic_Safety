function Func_HC_clustering(Input_Signal,x,Y_criteria,NumIgnoreCluster, plot_)

norm_All_nm = normalize(Input_Signal,2,"range");
Input_Signal_smooth = movmean(Input_Signal,5,2);
norm_Input_Signal = normalize(Input_Signal_smooth,2,"range");

% Step 1: Compute pairwise distances between signals (1500 signals, each of length 3000)
distances = pdist(norm_Input_Signal, 'euclidean'); 

% Step 2: Perform hierarchical clustering using the linkage method
Z = linkage(distances, 'ward');  % 'ward' method is used here, but other methods can be tried

% Step 3: Plot the dendrogram to visualize the hierarchical clustering
figure;
dendrogram(Z,0);
grid on
title('Hierarchical Binary Cluster Tree', 'Interpreter', 'latex');
set(gca, 'XTickLabel', []);
xlabel('Leaf Nodes,~Index of Records (not shown)', 'Interpreter', 'latex');
ylabel('Distance', 'Interpreter', 'latex')

% Step 4: Create clusters from the linkage result
if Y_criteria <= 0
    Distance = input('Please enter a number: ');
else
    Distance = Y_criteria;
end
clusters = cluster(Z, 'cutoff',Distance,'criterion','distance');
h = yline(Distance,'r--',['Threshold = ',num2str(Distance)]);
h.FontSize = 8;  % Set font size to 14

validIndices = true(size(clusters));
numClusters = unique(length(clusters));
numClusters_sec = numClusters;
% Loop through each cluster and check its size
for i = 1:numClusters
    % Find the signals belonging to cluster 'i'
    clusterIndices = find(clusters == i);
    
    % Check if the number of signals in this cluster is smaller than 5
    if length(clusterIndices) < NumIgnoreCluster
        % Mark these signals for removal
        validIndices(clusterIndices) = false;
        numClusters_sec = numClusters_sec -1;
    end
end

% Remove the signals belonging to small clusters
norm_Input_Signal = norm_Input_Signal(validIndices, :);
norm_All_nm = norm_All_nm(validIndices, :);
Input_Signal = Input_Signal(validIndices, :);
clusters = clusters(validIndices);
numClusters = max(clusters);  % Number of clusters


if   size(x) == size(Input_Signal(1, :))
    f_s = x;
else
    f_s = transpose(x);
end

if plot_ == 1
    for i = 1:numClusters
        % Find the indices of the signals that belong to cluster 'i'
        clusterIndices = find(clusters == i);
        
        if  isempty(clusterIndices)
            continue;
        end
        
        disp(['Number of cluster: ',num2str(i)]);
        
        % Create a new figure for this cluster
        figure;
        sgtitle(['Cluster ' num2str(i) ', $N =$' num2str(length(clusterIndices))], 'Interpreter', 'latex');
        
        % Plot each signal in this cluster
        subplot(2,1,1)

        for j = 1:length(clusterIndices)
            energy = trapz(Input_Signal(clusterIndices(j), :));
            plot(f_s, Input_Signal(clusterIndices(j), :)/energy);
            hold on;
        end
        title('Normalized PSD by Energy $\tilde{S}_{x}(f)$', 'Interpreter', 'latex')
        xlabel('freq (Hz)', 'Interpreter', 'latex');
        ylabel('$\tilde{S}_{x}(f)$', 'Interpreter', 'latex');
    
        subplot(2,1,2)
        %norm_All_nm = normalize(Input_Signal(clusterIndices, :),2,"range");
        std_dev =  std(norm_All_nm(clusterIndices, :),0,1);
        plot(f_s, mean(norm_All_nm(clusterIndices, :),1),'k');
        hold on;
        plot(f_s ,mean(norm_All_nm(clusterIndices, :),1)+1*std_dev,'c');
        plot(f_s, mean(norm_All_nm(clusterIndices, :),1)+3*std_dev,'b');

        for j = 1:length(clusterIndices)
            plot(f_s, norm_Input_Signal(clusterIndices(j), :),'r:');
        end

        plot(f_s, mean(norm_All_nm(clusterIndices, :),1),'k');
        plot(f_s ,mean(norm_All_nm(clusterIndices, :),1)+1*std_dev,'c');
        plot(f_s, mean(norm_All_nm(clusterIndices, :),1)+3*std_dev,'b');

        legend({'$\mathbf{E}[\hat{S}_{x}(f)]$','$\mathbf{E}[\hat{S}_{x}(f)] + 1\sigma$','$\mathbf{E}[\hat{S}_{x}(f)]+ 3\sigma$','$\hat{S}_{x}(f)$'}, 'Interpreter', 'latex');
        title('Normalized PSD by Maximum Amplitude $\hat{S}_{x}(f)$', 'Interpreter', 'latex')
        xlabel('freq (Hz)', 'Interpreter', 'latex');
        ylabel('$\hat{S}_{x}(f)$', 'Interpreter', 'latex');

        hold off;
    end
end

