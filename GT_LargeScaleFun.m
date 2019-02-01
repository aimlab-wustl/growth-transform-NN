function [y, yInit, spk] = GT_LargeScaleFun(Q,I_input,nIter,yInit)
%GROWTHTRANSFORMNEURONS Growth Transform Neuron Model
% 	This function creates a GUI for simulating large network of Neuron
% 	Models based on
% growth transforms (Ref. 1, 2). More details about the model and its
% applications can be found in Ref. 3. In brief, the Growth Transform (GT)
% Neuron Network generates spiking activity while minimizing an objective
% function. The GUI simulates a network of neurons with inhibitory and
% excitatory connections. Inputs to the function are:
% 
% Q - connectivity matrix square matrix of size = number of neurons with 
% diagnoal elements = 0
% I_input - input current 
% nIter - number of iterations 
% yInit - the initial state

% Copyright (c) [2018] Washington University  in St. Louis Created by:
% [Darshit Mehta, Ahana Gangopadhyay, Kenji Aono, Shantanu Chakrabartty] 1.
% Gangopadhyay, A., and Chakrabartty, S. (2017). Spiking, bursting, and
% population dynamics in a network of growth transform neurons. 2.
% Gangopadhyay, A., Chatterjee, O., and Chakrabartty, S. (2017). Extended
% polynomial growth transforms for design and training of generalized
% support vector machines. IEEE Transactions on Neural Networks and
% Learning Systems 3.  Gangopadhyay, A., Aono, K.  Mehta, D., and
% Chakrabartty, S. (in Review). A Coupled Network of Growth Transform
% Neurons for Spike-Encoded Auditory Feature Extraction
% 
% Washington University hereby grants to you a non-transferable,
% non-exclusive, royalty-free, non-commercial, research license to use and
% copy the computer code provided here (the “Software”).  You agree to
% include this license and the above copyright notice in all copies of the
% Software.  The Software may not be distributed, shared, or transferred to
% any third party.  This license does not grant any rights or licenses to
% any other patents, copyrights, or other forms of intellectual property
% owned or controlled by Washington University.  If interested in obtaining
% a commercial license, please contact Washington University's Office of
% Technology Management (otm@dom.wustl.edu).
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS
% PROVIDED “AS IS”, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED,
% INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR
% ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT,
% COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  IN NO EVENT SHALL THE
% CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY
% DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN
% ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS
% AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH
% PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES. YOU ALSO AGREE THAT
% THIS SOFTWARE WILL NOT BE USED FOR CLINICAL PURPOSES.

nNeuron = size(Q,1);
epsilon = 0.1;
% Set synaptic weight matrix
% Uncoupled case

% Spike parameters
w = 5;
gd = 0.5;
u = gd + epsilon;
l = gd - epsilon;
% P(:,1) = l + yInit/2-0.001;
% P(:,2) = u - yInit/2+0.001;
thr = 0;





% Initialize iteration variables

% Convergence hyperparameters
C = 0*ones(nNeuron,1);  %Regularization hyper-parameter
Fac = 10;
a = ones(nNeuron,2);
%spikecount = zeros(nNeuron,1);
exp_ad = 9.99*ones(nNeuron,1);
y = zeros(nNeuron,nIter);
y(:,1) = yInit;
dp = zeros(nNeuron,1);
dpprev = dp;

for c1 = 2:nIter
    
%     I = I_input-0.005+0.002*randn(size(I_input));
% %     biter = 2*u-I;
%     biter = 2*(1+epsilon)-I;
%     bias = [0.5*biter -0.5*biter];
%     Pprev = P;
%     C = C.*exp(-exp_ad.*(dp>thr)/100).*exp((1-C)/100); % Adaptation
% %     quant = w*(0.5<P & P<(0.5+epsilon)) - sec.*(P<0.5);
%     quant = (-(P<=l)+w*sign(P-(1/2)).*(l<P & P<u)+(P>=u));
%     G = - (bias+P+Q*[P(:,1) P(:,2)-2*epsilon]) - diag(C)*quant; % Gradient
% %     Fac = -min(G(:))+0.1;
%     
%     if any(G+Fac)<0
%         fprintf('No growth\n')
%     end
%     P = P.*(G + Fac); % Growth Transform
%     P = P./(sum(P,2)*ones(1,2)); % Normalization
%     a_iter = a; % Step size for axonal delay
%     a_iter(dp>thr,:) = 1; %Step size = 1, if neuron is spiking
%     P = a_iter.*P+(1-a_iter).*Pprev;
%     dp = P(:,1)-P(:,2);
%     y(:,c1) = P(:,1)-P(:,2)+2*epsilon;
    
    
     % External stimuli current
     I = I_input-0.007+0.005*randn(size(I_input));  
     %netI = I_input+ac_amp.*sin(2*pi*freq*iter/1000)+0.3*pulseFlag; % Net input current
 
     ind = (dp > thr);        
     C = exp_ad/10.*C + 1*ind;                               
     TotInp = 1*Q*C - 0.1*C;
     a_iter = 0.5*(1+0.95*tanh(2*TotInp));

     a_iter(ind,:) = 1;        
     %spikecount(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) + 1) + (1-(burstFlag(ind,:) > 0)).*(spikecount(ind,:));
     %sec(ind,:) = (burstFlag(ind,:) > 0).*(a(ind,:).*(spikecount(ind,:) > 10) + 1.5*(spikecount(ind,:) <= 10)) + a(ind,:).*(1-(burstFlag(ind,:)>0));
     %dec(ind,:) = (burstFlag(ind,:) > 0).*(0*(spikecount(ind,:) > 10) + 1*(spikecount(ind,:) <= 10)) + 0*(1-(burstFlag(ind,:)>0));        
     %spikecount(ind,:) = 0*(spikecount(ind,:) > 10) + spikecount(ind,:).*(spikecount(ind,:) <= 10);

        
%        y = [y(:,2:end), dp+ 0.5*(ind)];
     dpprev = dp; 
%        quant = sec.*(ind) - dec.*(dp <= thr);
     quant = 3*(ind) - 0*(dp <= thr);        
     G = 1*I - Q*dp - 1*dp - quant; % Gradient - self-inhibitory  
     dp = (G + Fac.*dp)./(dp.*G + Fac);
     if any(abs(dp))> 1
        fprintf('No growth\n')
     end
     dp = a_iter.*dp + (1-a_iter).*dpprev;
     y(:,c1) = dp + 0.5*(ind); 
        
%        I_hist = [I_hist(:,2:end), netI];
%
%        iter = mod(iter,100000)+1;
%        burstIter = burstIter+1;  
    
end
yInit = y(:,end);
spk = y(:,2:end)>0;