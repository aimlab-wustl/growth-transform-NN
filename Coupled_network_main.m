% Plots response trajectories for a 2-neuron network 
% in the absence and presence of pre-synaptic adaptation
% Corresponds to Figure 5a) of the paper:
% 'A Spiking Neuron and Population Model based on the Growth Transform Dynamical System', bioRxiv (2019), 523944

numPlots = 2;
N = 2;
maxiter = 1000;
xaxis = 1:maxiter;

Q = [0 0.5; -0.3 0];

% Input current stimulus
I_input = [0.3 0.04]';

burstFlag = [0 0; 0 0];
adaptFlag = [1 1; 0 0];
inh_fac = [0.1 0.1; 0.1 0.1];
a_fac = [0.5 0.5; 0.5 0.5];
exp_ad_fac = [0 0; 0 0];
lim = 1;
dp0 = 0.8*[-lim -lim; -lim -lim];

cc = [0 0 0; 0.5 0.5 0.5];

figure(1); hold on
[X, Y] = meshgrid(-lim:0.001:lim, -lim:0.001:lim);
Z = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X,2)
        dp = [X(i,j) Y(i,j)]';
        Z(i, j) = dp'*(Q+eye(N,N))*dp - I_input'*dp + 3*sum(dp.*(dp >0));
    end
end
colormap gray;
cmap = colormap;
contourf(X, Y, Z, 20, 'edgecolor', 'none', 'HandleVisibility', 'Off'); colormap(flipud(cmap));
ind = Z==0;



for i = 1:numPlots
    
   [dpf, dpcheckf] = Coupled_network(Q, I_input, burstFlag(i,:)', adaptFlag(i,:)', maxiter, inh_fac(i,:)', a_fac(i,:)', exp_ad_fac(i,:)', dp0(i,:)');
   figure; hold on
   for n = 1:N
      plot(xaxis,dpcheckf(n,:)-1.5*n,'k','LineWidth',2);  
      set(gca, 'XTick', [], 'YTick', []);
      set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1]);
      set(gcf, 'Color', [1 1 1]);
      if i==1
          title('With pre-synaptic adaptation');
      else 
          title('Without pre-synaptic adaptation');
      end
   end
   hold off
   figure(1); plot(dpf(1,:), dpf(2,:), 'Color', cc(i,:), 'LineWidth',2); 
   
end

legend('With adaptation', 'Without adaptation');
line([0, 0], [0, -lim], 'Color', 'r', 'LineStyle', '--', 'HandleVisibility', 'Off');
hold on
line([-lim, 0], [0, 0], 'Color', 'r', 'LineStyle', '--', 'HandleVisibility', 'Off');
hold off