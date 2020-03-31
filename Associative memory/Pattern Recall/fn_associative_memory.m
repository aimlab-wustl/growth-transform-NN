function [accuracy, num_spikes] = fn_associative_memory(P, s);

N = 100;   % Total number of neurons in the population
k = 10;    % Number of active neurons
l = 5;     % Overlap with a memory pattern during recall

% Create P memory patterns with k randomly selected active cells for each
rng(s);
B = zeros(N, P);
for num = 1:P
    indices = randperm(N, k);
    B(indices, num) = 1;
end

% Create synaptic weight matrix using Willshaw model
Q = B*B';
Q = Q./max(max(Q));

maxiter = 100;
input0 = - 0.5;
input1 = 1;

% Firing patterns corresponding to original memories
rasters = [];
for num = 1:P
    
    pattern = B(:, num);
    indices = find(pattern == 1);
    
    % Input current stimulus
    I_input = input0*ones(N, 1);
    I_input(indices) = input1;
    I_input = repmat(I_input, 1, maxiter);
    
    % Run growth transform network
    raster = gtnn(I_input, Q);
    rasters = [rasters raster];
    
end


% Recall/retrieval phase
rasters_recall = [];
allindices = zeros(k, P);
chosenindices = zeros(l, P);
deletedindices = zeros(k-l, P);
for num = 1:P
    
    pattern = B(:, num);
    indices = find(pattern == 1);
    allindices(:, num) = indices;
    removeindices = randperm(k, k-l);
    deletedindices(:, num) = indices(removeindices);
    indices(removeindices) = [];
    chosenindices(:, num) = indices;
    
    % Input current stimulus
    I_input = input0*ones(N, 1);
    I_input(indices) = input1;
    I_input = repmat(I_input, 1, maxiter);
    
    % Run growth transform network
    raster = gtnn(I_input, Q);
    rasters_recall = [rasters_recall raster];
    
end


% Check similarity between original and retrieved memories
bin = round(0.1*maxiter);
numbins = maxiter/bin;
ratemat = zeros(N, P*maxiter/bin);
ratemat_recall = zeros(N, P*maxiter/bin);
isimat = zeros(N, P*maxiter/bin);
isimat_recall = zeros(N, P*maxiter/bin);
deltaratemat = zeros(N, P*maxiter/bin);
deltaratemat_recall = zeros(N, P*maxiter/bin);
count = 1;
for j = 1:bin:(P*maxiter)-bin+1
    r = rasters(:, j:(j+bin-1));
    r_recall = rasters_recall(:, j:(j+bin-1));
    ratemat(:, count) = sum(r, 2);
    ratemat_recall(:, count) = sum(r_recall, 2);
    for i = 1:N
        spiketimes = find(r(i,:) == 1);
        spiketimes_recall = find(r_recall(i,:) == 1);
        if isempty(spiketimes)
            isimat(i, count) = 2*bin;
        elseif length(spiketimes) == 1
            isimat(i, count) = bin;
        else
            isimat(i, count) = mean(diff(spiketimes));
        end
        
        if isempty(spiketimes_recall)
            isimat_recall(i, count) = 2*bin;
        elseif length(spiketimes_recall) == 1
            isimat_recall(i, count) = bin;
        else
           isimat_recall(i, count) = mean(diff(spiketimes_recall));
        end
    end
    if count > 1
        deltaratemat(:, count) = ratemat(:, count) - ratemat(:, count-1);
        deltaratemat_recall(:, count) = ratemat_recall(:, count) - ratemat_recall(:, count-1);
    end
    count = count + 1;
end


% Compute mean distance between each pair of original-recall trajectories 
sim_matrix = zeros(P, P);
for i = 1:P       %for each recall pattern
    for j = 1:P   %for each original pattern 
        orig_traj = [ratemat(:, (j-1)*numbins+1:(j*numbins)); isimat(:, (j-1)*numbins+1:(j*numbins)); ...
            deltaratemat(:, (j-1)*numbins+1:(j*numbins))];    
        recall_traj = [ratemat_recall(:, (i-1)*numbins+1:(i*numbins)); isimat_recall(:, (i-1)*numbins+1:(i*numbins)); ...
            deltaratemat_recall(:, (i-1)*numbins+1:(i*numbins))];   
        dist = zeros(size(orig_traj, 2), 1);
        for k = 1:size(orig_traj, 2)
           dist(k, 1) = norm(orig_traj(:, k) - recall_traj(:, k)); 
        end
        sim_matrix(i, j) = mean(dist);
    end
end


% Recall accuracy
accuracy = 0;
for i = 1:P
   [~, ind] = min(sim_matrix(i,:));
   if i==ind
       accuracy = accuracy + 1;
   end
end
accuracy = (accuracy/P)*100;
num_spikes = sum(sum(rasters_recall));