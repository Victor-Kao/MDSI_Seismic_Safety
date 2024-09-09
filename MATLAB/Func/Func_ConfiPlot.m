function [mean_signal,margin_of_error] = Func_ConfiPlot(signals,x,plot_)

%%
% Signals = Size[N,:]

% Step 1: Compute the mean across all signals
mean_signal = mean(signals, 1);

% Step 2: Compute the standard deviation across all signals at each point
std_signal = std(signals, 0, 2);

% Step 3: Compute the standard error of the mean (SEM)
n_signals = size(signals, 1);
sem_signal = std_signal / sqrt(n_signals);

% Step 4: Compute the t-statistic for 95% confidence interval with (n-1) degrees of freedom
confidence_level = 0.95;
t_value = tinv(confidence_level + (1-confidence_level)/2, n_signals - 1);

% Step 5: Compute the margin of error for the 95% confidence interval
margin_of_error = t_value * sem_signal;


if size(x) ~= size(mean_signal)
    x = transpose(x);
end
% Step 6: Plot the mean signal and the 95% confidence interval
%x = 1:N;  % Time or data point axis

if plot_
    figure;   
    % Plot the confidence interval as a shaded area
    fill([x, fliplr(x)], [mean_signal + margin_of_error, fliplr(mean_signal - margin_of_error)], ...
        'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on 

    plot(x, mean_signal,'k', 'LineWidth', 2);  % Plot the mean signal
    
    xlabel('Data Points');
    ylabel('Amplitude');
    title('Mean Signal with 95% Confidence Interval');
    hold off 
end

end

