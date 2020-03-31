function indf = gtnn_globaladapt(I_input, Q)

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
flag = 1;
iter = 1;
Cr = 0;
H = 0;
Hdot = 0;
Wr = round(0.1*maxiter);   %Window to calculate energy
Crnl = 0;


% Store iteration variables
dpcheckf = zeros(N,maxiter);
dpf = zeros(N,maxiter);
indf = zeros(N, maxiter);
Encheckf = zeros(1, maxiter);

    
while iter <= maxiter
    
    netI = I_input(:,iter);
    ind = (dp > thr);
    indf(:, iter) = ind;
    dp(dp > thr) = thr;
    
    if rem(iter, Wr) == 0 && flag ==1
        H = [H sum(Encheckf(:, iter-Wr+1:iter))];
        Hdot = [Hdot H(:, end) - H(:, end-1)];
        if abs(H(end))>100 && abs(Hdot(end)) <= 10
            Cr = Cr + 0.1;
            flag = 0;
        end
        Crnl = tanh(Cr);
    end
    
    C = (9.9+0.1*log2(1+exp_ad))/10.*C + 1*ind;
    TotInp = (1/(N*N)).*(1*Q*C) - 0.1*C - Cr;
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
    
    En = 0.5*dp'*Q*dp - netI'*dp + 3*sum(dp.*(dp >0));
    Encheckf(1, iter) = En;
    
    iter = iter+1;
    
end
    

