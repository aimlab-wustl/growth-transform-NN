% Simulations showing the effect of pre-synaptic adaptation 
% (excitatory and inhibitory coupling) 
% Corresponds to Figure 5a) of the paper:
% 'A Spiking Neuron and Population Model based on the Growth Transform Dynamical System', bioRxiv (2019), 523944

% Number of neurons in the network
N = 2;

thr = 0.0;
Fac = 12*ones(N,1);
dp = -0.1*ones(N,1);
dpprev = dp;
spikecount = zeros(N,1);
sec = 3*ones(N,1);
dec = 0*ones(N,1);
bcount = 3;
C = ones(N,1); 
exp_ad = 0;
a = 3*ones(N, 1);

% Flags for different modes
burstFlag = [0 0]';
adaptFlag = [1 1]';

% Set synaptic weight matrix
Q = [1 0.5; 0 1];  % excitatory pre-synaptic adaptation
%Q = [1 -0.5; 0 1];  % inhibitory pre-synaptic adaptation

% Initialize iteration variables
iter = 1;
maxiter = 250;

% Store iteration variables
dpcheckf = zeros(N,maxiter);
dpf = zeros(N,maxiter);

% Input current stimulus
start_perturb = 0.1*maxiter;
stop_perturb = 1*maxiter;
b = zeros(N,maxiter);
b(1, start_perturb:stop_perturb) = 0.05;
b(2, 0.5*maxiter:stop_perturb) = 0.05;            


while iter <= maxiter
    
    netI = b(:,iter);
    ind = (dp > thr);
    dp(dp > thr) = thr;
    
    C = (9.9+0.1*log2(1+exp_ad))/10.*C + 1*ind;
    TotInp = 1*Q*C - 0.1*C;
    a_iter = 0.5*(1+0.95*tanh(2*TotInp)).*(adaptFlag > 0) + 0.5*ones(N, 1).*((1 - adaptFlag) > 0);
    
    
    spikecount(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) + 1) + (1-(burstFlag(ind,:) > 0)).*(spikecount(ind,:));
    sec(ind,:) = (burstFlag(ind,:) > 0).*(a(ind,:).*(spikecount(ind,:) > bcount) + 1.5*(spikecount(ind,:) <= bcount)) + 3.*(1-(burstFlag(ind,:)>0));
    dec(ind,:) = (burstFlag(ind,:) > 0).*(0*(spikecount(ind,:) > bcount) + 1*(spikecount(ind,:) <= bcount)) + 0*(1-(burstFlag(ind,:)>0));        
    spikecount(ind,:) = 0*(spikecount(ind,:) > bcount) + spikecount(ind,:).*(spikecount(ind,:) <= bcount);
    a_iter(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) > bcount) + 0.5*(burstFlag(ind,:) > 0).*(spikecount(ind,:) < bcount) + (1-(burstFlag(ind,:)>0));
    
    dpcheckf(:, iter) = dp + ind;
    dpf(:, iter) = dp;
    dpprev = dp; 
    quant = sec.*(ind) - dec.*(dp <= thr);
    
    % G = - delH/delV_i
    G = -Q*dp + netI - quant;
    dp = (G + Fac.*dp)./(dp.*G + Fac);
    dp = a_iter.*dp + (1-a_iter).*dpprev;
    
    iter = iter+1;
    
end

% Response of two neurons
figure; 
set(gcf, 'Color', 'w');
xaxis = 1:maxiter;
subplot(2,1,1); hold on
plot(xaxis, dpcheckf(2,:), 'k', 'LineWidth', 4); 
plot(xaxis, 5*b(2, :)-1, 'b', 'LineWidth', 4);
set(gca,'XTick',[],'YTick',[]); box off;
set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
ylim([-1.2 1.2]);
title('Pre-synaptic neuron');
subplot(2,1,2); hold on
plot(xaxis, dpcheckf(1,:), 'k', 'LineWidth', 4); 
plot(xaxis, 5*b(1, :)-1, 'b', 'LineWidth', 4);
set(gca,'XTick',[],'YTick',[]); box off;
set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
ylim([-1.2 1.2]);
title('Post-synaptic neuron');
hold off