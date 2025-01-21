% Taylor Series Approximation for Pi with Convergence Plot
clear; clc;

% Define range of terms
max_terms = 100000; % Maximum number of terms
pi_approximations = zeros(1, max_terms); % Preallocate array for approximations

% Compute approximations for each number of terms
for N = 1:max_terms
    % Calculate the approximation of Pi using N terms
    pi_approx = 4 * sum((-1).^(0:N-1) ./ (2*(0:N-1) + 1));
    pi_approximations(N) = pi_approx;
end

% Plot the results
figure;
plot(1:max_terms, pi_approximations, 'b-', 'LineWidth', 1);
hold on;
yline(pi, 'r--', 'LineWidth', 1.5, 'DisplayName', 'True Value of \pi');
xlabel('Number of Terms');
ylabel('Approximation of \pi');
title('Convergence of Pi Approximation using Taylor Series');
legend('Approximation', 'True Value of \pi');
grid on;


