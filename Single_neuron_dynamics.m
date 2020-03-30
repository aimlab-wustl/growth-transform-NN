% Simulate single-neuron dynamics
maxiter = 100;

figure; hold on
set(gcf, 'Color', 'w');

% Tonic spiking
b = zeros(1, maxiter);
b(:, 0.1*maxiter+1 : end) = 0.5;
s = gtnn(1, b, 0, 0, 0);
subplot(2, 2, 1); hold on
plot(s, 'k', 'LineWidth', 2);
plot(b-1.5, 'b', 'LineWidth', 2);
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
ylim([-1.7 1.5]);
title('Tonic spiking');


% Bursting
b = zeros(1, maxiter);
b(:, 0.1*maxiter+1 : end) = 0.5;
s = gtnn(1, b, 0, 0, 1);
subplot(2, 2, 2); hold on
plot(s, 'k', 'LineWidth', 2);
plot(b-1.5, 'b', 'LineWidth', 2);
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
ylim([-1.7 1.5]);
title('Bursting');


% Spike-frequency adaptation
b = zeros(1, maxiter);
b(:, 0.1*maxiter+1 : end) = 0.5;
s = gtnn(1, b, 1, 0, 0);
subplot(2, 2, 3); hold on
plot(s, 'k', 'LineWidth', 2);
plot(b-1.5, 'b', 'LineWidth', 2);
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
ylim([-1.7 1.5]);
title('Spike-frequency adaptation');


% Integration
pulseWidth = 5; pulse = 0.2;
b = - 0.1*ones(1, maxiter);
b(:, 10:10+pulseWidth) = pulse; b(:, 40:40+pulseWidth) = pulse;
b(:, 70:70+pulseWidth) = pulse; b(:, 80:80+pulseWidth) = pulse;
s = gtnn(1, b, 0, 0, 0);
subplot(2, 2, 4); hold on
plot(s, 'k', 'LineWidth', 2);
plot(b-1.5, 'b', 'LineWidth', 2);
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
ylim([-1.8 1.5]);
title('Integration');