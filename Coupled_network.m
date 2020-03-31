function [dpf, dpcheckf] = Coupled_network(Q, netI, burstFlag, adaptFlag, maxiter, inh_fac, a_fac, exp_ad, dp)

N = size(Q, 1);
thr = 0.0;
Fac = 12*ones(N,1);
spikecount = zeros(N,1);
sec = 3*ones(N,1);
dec = 0*ones(N,1);
bcount = 3;
C = ones(N,1); 
a = 3*ones(N, 1);
I = eye(N,N);

% Initialize iteration variables
iter = 1;


% Store iteration variables
dpf = zeros(N,maxiter);
dpcheckf = zeros(N,maxiter);


%% Run growth-transform updates

while iter <= maxiter
    
    ind = (dp > thr);
    
    C = (9.9+0.1*log2(1+exp_ad))/10.*C + 1*ind;
    TotInp = 1*Q*C - inh_fac.*C;
    a_iter = a_fac.*(1+0.95*tanh(2*TotInp)).*(adaptFlag > 0) + 0.5*ones(N, 1).*((1 - adaptFlag) > 0);
    
    %a_iter(ind,:) = 1;
    spikecount(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) + 1) + (1-(burstFlag(ind,:) > 0)).*(spikecount(ind,:));
    sec(ind,:) = (burstFlag(ind,:) > 0).*(a(ind,:).*(spikecount(ind,:) > bcount) + 1.5*(spikecount(ind,:) <= bcount)) + 3.*(1-(burstFlag(ind,:)>0));
    dec(ind,:) = (burstFlag(ind,:) > 0).*(0*(spikecount(ind,:) > bcount) + 1*(spikecount(ind,:) <= bcount)) + 0*(1-(burstFlag(ind,:)>0));        
    spikecount(ind,:) = 0*(spikecount(ind,:) > bcount) + spikecount(ind,:).*(spikecount(ind,:) <= bcount);
    a_iter(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) > bcount) + 0.5*(burstFlag(ind,:) > 0).*(spikecount(ind,:) < bcount) + (1-(burstFlag(ind,:)>0));
    %a_iter = 0.5*ones(N, 1);
    
    dpf(:, iter) = dp;
    dpcheckf(:, iter) = dp + 0.5*(ind);
    dpprev = dp; 
    quant = sec.*(ind) - dec.*(dp <= thr);
    G = 1*netI - (Q+I)*dp - 1*dp - quant; % Gradient - self-inhibitory  
    dp = (G + Fac.*dp)./(dp.*G + Fac);
    dp = a_iter.*dp + (1-a_iter).*dpprev;
    
    
    iter = iter+1;
    
end

% % Response of one neuron 
% figure; 
% xaxis = 1:maxiter;
% subplot(2,1,1);  plot(xaxis,dpcheckf(1,:),'k','LineWidth',4); 
% set(gca,'XTick',[],'YTick',[]); box off;
% subplot(2,1,2);  plot(xaxis,dpcheckf(2,:),'k','LineWidth',4);   
% set(gca,'XTick',[],'YTick',[]); box off;
% hold off


% % Network energy
% figure; 
% Enavg = tsmovavg(Encheckf, 's', W);
% plot(xaxis, Enavg(index, :), 'k', 'LineWidth', 4);