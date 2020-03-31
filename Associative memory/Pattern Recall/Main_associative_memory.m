% Plots recall accuracy and spike count vs number of patterns stored for
% 1-200 patterns for an associative network of 100 neurons
% Corresponds to Figures 11(a) and (b) of the paper:
% 'A Spiking Neuron and Population Model based on the Growth Transform Dynamical System', bioRxiv (2019), 523944

Ntimes = 10;
Pmax = 200;
svec = [12, 5, 4, 9, 36, 74, 100, 14, 45, 69, 21, 44, 11, 32, 8, 3, 1, 7, 25, 11];
accuracy_mat = zeros(Pmax, Ntimes);
num_spikes_mat = zeros(Pmax, Ntimes);
for P = 1:Pmax
    for i = 1:Ntimes
        [accuracy, num_spikes] = fn_associative_memory(P, svec(i));
        accuracy_mat(P, i) = accuracy;
        num_spikes_mat(P, i) = num_spikes;
    end
    P
end

accuracy_mat2 = zeros(Pmax, Ntimes);
num_spikes_mat2 = zeros(Pmax, Ntimes);
for P = 1:Pmax
    for i = 1:Ntimes
        [accuracy, num_spikes] = fn_associative_memory_globaladapt(P, svec(i));
        accuracy_mat2(P, i) = accuracy;
        num_spikes_mat2(P, i) = num_spikes;
    end
    P
end


figure; hold on
stdshade(accuracy_mat', 0.5, 'b');
stdshade(accuracy_mat2', 0.5, 'r');
xlabel('Number of patterns stored');
ylabel('Recall accuracy');
box on; grid on;
hold off;

figure; hold on
stdshade(num_spikes_mat', 0.5, 'b');
stdshade(num_spikes_mat2', 0.5, 'r');
xlabel('Number of patterns stored');
ylabel('Number of spikes');
box on; grid on;
hold off;