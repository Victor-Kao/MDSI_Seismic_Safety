clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");
load("D:/MDSI_project/DATA_GM_RawData/CellArrayAllEvent.mat");

%%
plot_ = 0;

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


% Perform SVD
[U, S, V] = svd(norm_All_s, 'econ');

% Initialize a superposed component vector
superposed_component = zeros(1, size(V, 2));

V = transpose(V);
% Use the first row as the reference for alignment
reference_row = V(1,:);

% Superpose the top 5 components (aligned with the reference)
for i = 1:4
    % Check correlation with the reference row
    if corr(V(i,:)', reference_row') < 0
        % Flip the row if the correlation is negative
        superposed_component = superposed_component + (-V(i,:));  % Flip and add
    else
        % Add the row directly if it's already aligned
        superposed_component = superposed_component + V(i,:);
    end
end

% Plot the superposed result
figure;
plot(f_s,superposed_component, 'LineWidth', 2);
xlabel('Data Points');
ylabel('Amplitude');
title('Superposed Top 5 Significant Components from SVD');



% Apply SVD
%[U, S, V] = svd(norm_All_s, 'econ');  % 'econ' gives the economy-sized SVD
%
% Singular values
singular_values = diag(S);

% Plot singular values to see the significance of components
figure;
plot(singular_values, 'o-');
xlabel('Component Number');
ylabel('Singular Value');
title('Singular Values of the norm_All_s');
%
%% Visualize the first few significant components
%figure;
%plot(V(:,1), 'LineWidth', 2);  % First component
%hold on;
%plot(V(:,2), 'LineWidth', 2);  % Second component
%plot(V(:,3), 'LineWidth', 2);  % Third component
%plot(V(:,4), 'LineWidth', 2);  % Third component
%plot(V(:,5), 'LineWidth', 2);  % Third component
%plot(V(:,6), 'LineWidth', 2);  % Third component
%xlabel('Data Points');
%ylabel('Component Value');
%legend('Component 1', 'Component 2', 'Component 3');
%title('First Three Components of the Frequency Series');
%hold off;
%
%
%figure
%% Superpose the first few significant components (e.g., top 3)
%superposed_component = zeros(size(V(:,1)));  % Same size as one singular vector
%
%% Sum the top 50 components
%for i = 1:4
%    superposed_component = superposed_component + V(:,i);
%end
%
%% Normalize if necessary (optional)
%superposed_component = superposed_component / norm(superposed_component);
%
%% Plot the superposed components
%figure;
%plot(f_s,superposed_component, 'LineWidth', 2);
%xlabel('Data Points');
%ylabel('Amplitude');
%title('Superposition of Top 3 Significant Components');
