function Q = createNet(nNeuron,l1N, netType)
Q = zeros(nNeuron);
tmp = full(0.1*sprandn(l1N, l1N,0.02));
idx0 = find(sum(tmp)==0);
for c1 = 1:length(idx0)
    idx = randperm(size(tmp,1),1);
    while idx == idx0(c1)
        disp(c1)
        idx = randperm(size(tmp,1),1);
    end
    tmp(idx,idx0(c1)) = 0.1*randn;
end
Q(1:l1N,1:l1N) = tmp;

tmp = full(0.1*sprandn(nNeuron-l1N, nNeuron-l1N,0.01));

Q(1+l1N:end,1+l1N:end) = tmp;

nConn = 20*l1N;
s = randi(l1N,nConn,1);
t = randi(nNeuron-l1N,nConn,1);

for c1 = 1:nConn
    Q(t(c1)+l1N,s(c1)) = 0.1*randn;
end
switch netType
    case 1
        
    case 2
        Q(3:l1N,1:2) = 0.1*randn(l1N-2,2);
end
Q(Q<0) = 2*Q(Q<0);
for nf = 1:nNeuron
    Q(nf,nf) = 0;
end