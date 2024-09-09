

% Assume signals is a matrix where each row is a signal [15 x N] (15 signals, N points per signal)
N = 3000;  % Length of each signal
signals = randn(15, N);  % Replace with actual signals data

% Step 1: Compute the mean across all signals
mean_signal = mean(signals, 1);

% Step 2: Compute the standard deviation across all signals at each point
std_signal = std(signals, 0, 1);

% Step 3: Compute the standard error of the mean (SEM)
n_signals = size(signals, 1);
sem_signal = std_signal / sqrt(n_signals);

% Step 4: Compute the t-statistic for 95% confidence interval with (n-1) degrees of freedom
confidence_level = 0.95;
t_value = tinv(confidence_level + (1-confidence_level)/2, n_signals - 1);

% Step 5: Compute the margin of error for the 95% confidence interval
margin_of_error = t_value * sem_signal;

% Step 6: Plot the mean signal and the 95% confidence interval
x = 1:N;  % Time or data point axis

figure;
plot(x, mean_signal, 'LineWidth', 2);  % Plot the mean signal
hold on;

% Plot the confidence interval as a shaded area
fill([x, fliplr(x)], [mean_signal + margin_of_error, fliplr(mean_signal - margin_of_error)], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('Data Points');
ylabel('Amplitude');
title('Mean Signal with 95% Confidence Interval');
hold off;