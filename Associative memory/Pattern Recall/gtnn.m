function indf = gtnn(I_input, Q)

[N, maxiter] = size(I_input);

% Growth transform network model parameters
thr = 0.0;
Fac = 12*ones(N,1);
dp = -0.5*ones(N,1);
spikecount = zeros(N,1);
sec = 3*ones(N,1);
dec = 0*ones(N,1);
bcount = 3;
C = ones(N,1); 
exp_ad = 0;
a = 3*ones(N, 1);

% Flags for different modes
burstFlag = zeros(N, 1);
adaptFlag = ones(N, 1);


% Initialize iteration variables
iter = 1;

% Store iteration variables
dpcheckf = zeros(N,maxiter);
dpf = zeros(N,maxiter);
indf = zeros(N, maxiter);

    
while iter <= maxiter
    
    netI = I_input(:,iter);
    ind = (dp > thr);
    indf(:, iter) = ind;
    dp(dp > thr) = thr;
    
    C = (9.9+0.1*log2(1+exp_ad))/10.*C + 1*ind;
    TotInp = 1*Q*C - 0.1*C;
    a_iter = 0.5*(1+0.95*tanh(2*TotInp)).*(adaptFlag > 0) + 0.5*ones(N, 1).*((1 - adaptFlag) > 0);
    
    a_iter(ind,:) = 1;
    spikecount(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) + 1) + (1-(burstFlag(ind,:) > 0)).*(spikecount(ind,:));
    sec(ind,:) = (burstFlag(ind,:) > 0).*(a(ind,:).*(spikecount(ind,:) > bcount) + 1.5*(spikecount(ind,:) <= bcount)) + 3.*(1-(burstFlag(ind,:)>0));
    dec(ind,:) = (burstFlag(ind,:) > 0).*(0*(spikecount(ind,:) > bcount) + 1*(spikecount(ind,:) <= bcount)) + 0*(1-(burstFlag(ind,:)>0));
    spikecount(ind,:) = 0*(spikecount(ind,:) > bcount) + spikecount(ind,:).*(spikecount(ind,:) <= bcount);
    
    dpcheckf(:, iter) = dp + 0.5*(ind);
    dpf(:, iter) = dp;
    dpprev = dp;
    quant = sec.*(ind) - dec.*(dp <= thr);
    
    % G is actually the negative of delH/delV_i)
    G = -Q*dp + netI.*(dp <= thr) - quant;
    dp = (G + Fac.*dp)./(dp.*G + Fac);
    dp = a_iter.*dp + (1-a_iter).*dpprev;
    
    iter = iter+1;
    
end
    

