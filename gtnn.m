function dpcheckf = gtnn(Q, b, adaptFlag1, adaptFlag2, burstFlag)

[N, maxiter] = size(b);
thr = 0.0;
Fac = 12*ones(N,1);
dp = -0.1*ones(N,1);
dpprev = dp;
spikecount = zeros(N,1);
sec = 10*ones(N,1);
dec = 0*ones(N,1);
bcount = 3;
C = ones(N,1); 
a = 3;

% Initialize iteration variables
iter = 1;

% Store iteration variables
dpcheckf = zeros(N,maxiter);
dpf = zeros(N,maxiter);
indf = zeros(N, maxiter);


while iter <= maxiter
    
    netI = b(:,iter);
    ind = (dp > thr);
    indf(:, iter) = ind;
    dp(dp > thr) = thr;
    
    C = 9.9/10.*C + 2*ind;
    TotInp = Q*C;
    a_iter = (0.5*(1+0.95*tanh(-0.2*C))).*(adaptFlag1 > 0) + 0.5*(1+0.95*tanh(2*TotInp)).*(adaptFlag2 > 0);
    if adaptFlag1 == 0 && adaptFlag2 == 0
       a_iter = 0.5*ones(N, 1); 
    end
    
    
    spikecount(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) + 1) + (1-(burstFlag(ind,:) > 0)).*(spikecount(ind,:));
    sec(ind,:) = (burstFlag(ind,:) > 0).*(a(ind,:).*(spikecount(ind,:) > bcount) + 1.5*(spikecount(ind,:) <= bcount)) + 3.*(1-(burstFlag(ind,:)>0));
    dec(ind,:) = (burstFlag(ind,:) > 0).*(0*(spikecount(ind,:) > bcount) + 1*(spikecount(ind,:) <= bcount)) + 0*(1-(burstFlag(ind,:)>0));        
    spikecount(ind,:) = 0*(spikecount(ind,:) > bcount) + spikecount(ind,:).*(spikecount(ind,:) <= bcount);
    
    
    dpcheckf(:, iter) = dp + 1*(ind);
    dpf(:, iter) = dp;
    dpprev = dp; 
    quant = sec.*(ind) - dec.*(dp <= thr);

    
    % G = - delH/delV_i)
    G = -Q*(dp + 10*(ind)) + netI.*(dp <= thr) - quant;
    dp = (G + Fac.*dp)./(dp.*G + Fac);
    dp = a_iter.*dp + (1-a_iter).*dpprev;
    
    iter = iter+1;
    
end