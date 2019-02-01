%GROWTHTRANSFORMNEURONS Growth Transform Neuron Model
% 	This function creates a GUI for simulating Neuron Models based on 
% growth transforms (Ref. 1, 2). More details about the model and its
% applications can be found in Ref. 3. In brief, the Growth Transform (GT)
% Neuron Network generates spiking activity while minimizing an objective
% function. The GT Neron model is capable of emulating a number of dynamics
% which are observed in biological systems, such as tonal spiking in
% presence of constant input, excitation by sinsuidoilly varying input
% stimuls, leaky integrator, and as bursting neurons. All of these dynamics
% can be observed in the GUI by clicking the relevant buttons. Users have
% the ability to change the input current (DC and AC components), neuronal
% adaptation, bursting behavior, and the axonal delay. This can be done for
% an individual neuron or for all neurons simultaneously. Clicking reset
% resets all inputs to their default values.
% The GUI can also simulate a network of neurons with inhibitory and
% excitatory connections. There is an option for generating sparse random
% connection matrix in the GUI. The user can also input their own
% connection matrix 'Q' by defining it in the workspace and then importing
% it into the GUI using the Connectivity Matrix Dropdown Menu. The matrix
% should be called Q and it should be a square matrix of size = number of
% neurons. This GUI works best for 1-40 neurons. Swapping lines 255, 261
% and 273 for 254, 260 and 272 will simulate the network with a different
% spike mapping function.

% Copyright (c) [2018] Washington University  in St. Louis Created by:
% [Darshit Mehta, Ahana Gangopadhyay, Kenji Aono, Shantanu Chakrabartty]

% Citations for this tool are: 
% 1. Gangopadhyay, A., Mehta, D. and Chakrabartty, S. (2019). A Spiking
% Growth Transform Neuron and Population Model,. BioArxiv.

% 2. Gangopadhyay, A., and Chakrabartty, S. (2017). Spiking, bursting, and
% population dynamics in a network of growth transform neurons. IEEE Trans.
% Neural Network and Learning Systems.

% 3.  Gangopadhyay, A., Aono, K.  Mehta, D., and Chakrabartty, S. (2018). 
% A Coupled Network of Growth Transform Neurons for Spike-Encoded Auditory 
% Feature Extraction, BioArxiv.
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


function GT_NeuronModel
% function GT_NeuronModel
% Input number of neurons
prompt={'Enter number of neurons between 1 and 40'};
name = 'Number of neurons';
defaultans = {'5'};
answer = inputdlg(prompt,name,[1 40],defaultans);
temp = str2double(answer);
if isempty(answer)
    disp("Error in fetching number of neurons\n");
    return
elseif ~isnan(temp) && temp>0 && temp<41
    nNeuron = round(temp);
else
    disp('Number of neurons needs to be integer between 1 and 100')
    return
end
% Parameters
% Start from a random or fixed initial point


% Set synaptic weight matrix
% Uncoupled case
Q = zeros(nNeuron,nNeuron);
I = eye(nNeuron,nNeuron);
% Spike parameters
w = 5;
epsilon = 0.1;
gd = 0.5;
u = gd + epsilon;
l = gd - epsilon;
thr = 0;
% Initialize iteration variables
L = 1000;
P = [l*ones(nNeuron,1)-eps u*ones(nNeuron,1)+eps];
% Convergence hyperparameters
C = ones(nNeuron,1);  %Regularization hyper-parameter
Fac = 12;
a = 3*ones(nNeuron,1);
I_input = 0*ones(nNeuron,1);

nSpeed = 1;
ac_amp = zeros(nNeuron,1);
freq = 5*ones(nNeuron,1);
exp_ad = 0*ones(nNeuron,1);
y = zeros(nNeuron,L);
I_hist = zeros(nNeuron,L);
dp = zeros(nNeuron,1);

burstFlag = zeros(nNeuron,1);
burstIter = 50*ones(nNeuron,1);
pulseFlag = zeros(nNeuron,1);
nSelect = 1;
nDisp = 1;
pauseFlag = 0;
nStr = cell(1,nNeuron+1);
for n1 = 1:nNeuron
    nStr{n1} = num2str(n1);
end
nStr{nNeuron+1}='All';
% Figure
figNumber = figure(1);
clf;
set(figNumber,'NumberTitle','off',...
    'Name','Growth Transform Neuron Model',...
    'Units','normalized','toolbar','figure',...
    'Position',[0.05 0.1 0.9 0.8]);

% Control buttons

% Input fields
nInputFields = 7;
ycoord = linspace(1,0,nInputFields+1);
yHt = 1/nInputFields - 0.05;
pnl_input = uipanel('Title','Input Fields','FontSize',10,...
    'BackgroundColor','white',...
    'Position',[.31 .1 .35 .3]);
inputFieldStrings = {'Input current = ','AC Amplitude = ','AC Frequency = ','Adaptation = ','After-Poralization: ', 'Burst: '};
inputFieldDefaults = {'0','0','5','0','3',''};

inFieldDisp = cell(1,6);
for iter = 1:numel(inFieldDisp)
    uicontrol('Parent',pnl_input,'Style','text', 'Units','normalized', ...
        'BackgroundColor','white','Position',[0.2 ycoord(iter+2)+0.03 0.3 yHt],'string',inputFieldStrings{iter},...
        'FontSize',9, 'HorizontalAlignment', 'right');
    
    inFieldDisp{iter} = uicontrol('Parent',pnl_input,'Style','text', 'Units','normalized', ...
        'BackgroundColor','white','Position',[0.51 ycoord(iter+2)+0.03 0.1 yHt],'string',inputFieldDefaults{iter},...
        'FontSize',9, 'HorizontalAlignment', 'left');
end
inField{1} = uicontrol('Parent',pnl_input,'Style','slider',...
    'Min',-0.9,'Max',0.9, 'SliderStep',[0.05 0.2], 'Units','normalized', ...
    'Position',[0.62 ycoord(3)+0.04 0.3 yHt-0.01],...
    'tag','input','Callback',@changepars);


inField{2} = uicontrol('Parent',pnl_input,'Style','slider',...
    'Min',0,'Max',0.2, 'SliderStep',[0.05 0.2],'Units','normalized', ...
    'Position',[0.62 ycoord(4)+0.06 0.3 yHt-0.01],...
    'Value',ac_amp(1),'tag','ac_amp','Callback',@changepars);


inField{3} = uicontrol('Parent',pnl_input,'Style','slider',...
    'Min',0,'Max',10, 'SliderStep',[0.1 0.2],'Units','normalized',...
    'Position',[0.62 ycoord(5)+0.04 0.3 yHt-0.01],...
    'Value',freq(1),'tag','freq','Callback',@changepars);


inField{4} = uicontrol('Parent',pnl_input,'Style','slider',...
    'Min',0,'Max',0.999, 'SliderStep',[0.01 0.1], 'Units','normalized', ...
    'Position',[0.62 ycoord(6)+0.04 0.3 yHt-0.01],...
    'Value',exp_ad(1),'tag','exp_ad','Callback',@changepars);


inField{5} = uicontrol('Parent',pnl_input,'Style','slider',...
    'Min',0,'Max',5, 'SliderStep',[0.05 0.2], 'Units','normalized', ...
    'Position',[0.62 ycoord(7)+0.04 0.3 yHt-0.01],...
    'Value',a(1),'tag','alpha','Callback',@changepars);

inField{6} = uicontrol('Parent',pnl_input,'Style','checkbox',...
    'Min',0,'Max',1, 'Units','normalized', ...
    'Position',[0.51 ycoord(8)+0.04 0.05 yHt-0.01],...
    'BackgroundColor','white','tag','burst','Callback',@changepars);

uicontrol('Parent',pnl_input, 'Units','normalized', ...
        'Position',[0.62 ycoord(8)+0.02 0.1 yHt+0.01],...
        'string','Pulse','tag','Pulse','Callback',@changepars);


if nNeuron>1
    uicontrol('Parent',pnl_input,'Style','text', 'Units','normalized', ...
        'BackgroundColor','white','Position',[0.2 ycoord(2)+0.03 0.3 yHt],'string','Select neuron',...
        'HorizontalAlignment', 'right','FontSize',9);
    uicontrol('Parent',pnl_input,'Style', 'popup', 'Units','normalized',...
        'String', nStr,...
        'Position', [0.51 ycoord(2)+0.035 0.2 yHt],...
        'tag','neurons','Callback', @changepars);
end

spiking_mode_tags = {'Tonic Spiking','Tonal Excitation','Bursting','Integration', 'Adaptation'};
for iter = 1:length(spiking_mode_tags)
    uicontrol('Parent',pnl_input, 'Units','normalized', ...
        'Position',[0.02 0.95-0.17*iter 0.18 0.15],...
        'string',spiking_mode_tags{iter},'tag',spiking_mode_tags{iter},'Callback',@spiking_modes);
end

% Simulation speed
uicontrol('Style','text', 'Units','normalized', ...
    'Position',[0.7  0.96 0.2 0.03],'string','Simulation speed','FontSize',12);
uicontrol('Style', 'slider',...
    'Min',1,'Max',10,'Value',nSpeed, 'Units','normalized', ...
    'Position',[0.7  0.92 0.2 0.03],'SliderStep',[0.1 0.2],...
    'tag','speed','Callback',@changepars);

% Pause button
uicontrol('Style', 'togglebutton','String','Pause/Resume',...
    'Min',0,'Max',1,'Value',0, 'Units','normalized', ...
    'Position',[0.8 0.02 0.1 0.05],...
    'tag','pauseflag','Callback',@changepars);

% Reset button
uicontrol('Style', 'togglebutton','String','Reset',...
    'Min',0,'Max',1,'Value',0, 'Units','normalized', ...
    'Position',[0.65 0.02 0.1 0.05],...
    'tag','reset','Callback',@changepars);

h1 = axes('Position',[0.05 0.5 0.9 0.4]);
hold on

nAx = cell(1,nNeuron);
for n1 =1:nNeuron
    nAx{n1} = plot(h1,1:L,y(n1,:)+n1);
end



axis manual
axis([0 1000 0 nNeuron+1])
title('Membrane potential')
xlabel('Time (ms)');
ylabel('Neuron Index');
set(gca,'ytick',1:nNeuron)

h2 = axes('Position',[0.7 0.15 0.25 0.25]);
%I_Ax = plot(h2,1:L,I_hist(1,:));
E_av = zeros(1,1000);
I_Ax = plot(h2,1:1000,E_av);
axis manual
axis([0 1000 -1/nNeuron 1/nNeuron])
title('Network Energy')
xlabel('Time (ms)');
ylabel('Energy (a.u.)');
grid on;

colorMap = repmat(linspace(0,0.7,30)',1,3);
colorMap = [colorMap;[1 1 1];colorMap(end:-1:1,:)];
colorMap(1:30,1) = 1;
colorMap(32:61,3) = 1;
if nNeuron>1
    h3=axes('Position',[0.05 0.1 0.25 0.25]);
    conn_im = imagesc(h3,Q+I);
    colormap(colorMap)
    set(gca,'xtick',1:nNeuron,'ytick',1:nNeuron)
    ylabel('Post-synaptic')
    xlabel('Pre-synaptic')
    title('Connectivity Matrix')
    uicontrol('Style', 'popup', 'Units','normalized',...
        'String', {'Identity','Random Sparse','Random Non-sparse','From workspace'},...
        'Position', [0.08  0.22 0.2 0.2],...
        'tag','Q_mat','Callback', @changeQ);
    caxis([-1 1])
    caxis manual
    colorbar    
end
% Update and plot
thr = 0.0;
Fac = 12*ones(nNeuron,1);
iter = 1;
decay = 0.999;
dp = zeros(nNeuron,1);
dpprev = dp;
spikecount = zeros(nNeuron,1);
sec = 3*ones(nNeuron,1);
dec = 0*ones(nNeuron,1);
bcount = 3;
W = 900;
E_hist = zeros(1,W);

while ishandle(figNumber)
    for c1 = 1:nSpeed        
        
        % External stimuli current
        netI = I_input+ac_amp.*sin(2*pi*freq*iter/1000)+0.3*pulseFlag; % Net input current
        pulseFlag = zeros(nNeuron,1);

        ind = (dp > thr);
        
        C = (9.9+0.1*log2(1+exp_ad))/10.*C + 1*ind;               
                
        TotInp = 1*Q*C - 0.1*C;
        a_iter = 0.5*(1+0.95*tanh(2*TotInp));

%        a_iter(ind,:) = (burstFlag(ind,:) > 0).*((a_iter(ind,:).*(spikecount(ind,:) < 10)) + (spikecount(ind,:) >= 10)) + ...
%            (1-(burstFlag(ind,:)>0));
        a_iter(ind,:) = 1;        
        spikecount(ind,:) = (burstFlag(ind,:) > 0).*(spikecount(ind,:) + 1) + (1-(burstFlag(ind,:) > 0)).*(spikecount(ind,:));
        sec(ind,:) = (burstFlag(ind,:) > 0).*(a(ind,:).*(spikecount(ind,:) > bcount) + 1.5*(spikecount(ind,:) <= bcount)) + a(ind,:).*(1-(burstFlag(ind,:)>0));
        dec(ind,:) = (burstFlag(ind,:) > 0).*(0*(spikecount(ind,:) > bcount) + 1*(spikecount(ind,:) <= bcount)) + 0*(1-(burstFlag(ind,:)>0));        
        spikecount(ind,:) = 0*(spikecount(ind,:) > bcount) + spikecount(ind,:).*(spikecount(ind,:) <= bcount);

        
        y = [y(:,2:end), dp+ 0.5*(ind)];
        dpprev = dp; 
        quant = sec.*(ind) - dec.*(dp <= thr);
        G = 1*netI - Q*dp - 1*dp - quant; % Gradient - self-inhibitory  
        dp = (G + Fac.*dp)./(dp.*G + Fac);
        dp = a_iter.*dp + (1-a_iter).*dpprev;
        
        % Estimate the energy
        En = 0.5*dp'*(Q+I)*dp - netI'*dp + sum((sec.*dp).*(dp > 0));
        E_hist = [E_hist(2:end), En];        
        E_av = [E_av(2:end), sum(E_hist)/W];
               
%        I_hist = [I_hist(:,2:end), netI];

        iter = mod(iter,100000)+1;
        burstIter = burstIter+1;
    end
    for n1 =1:nNeuron
        set(nAx{n1},'ydata',y(n1,:)+n1) % Update the membrane potential plot
    end
    
%    set(I_Ax,'ydata',I_hist(nDisp,:)) % Update the input current plot
    set(I_Ax,'ydata',E_av) % Update the input current plot
    drawnow
    while pauseFlag && ishandle(figNumber)
        drawnow
        pause(0.1)
    end
end
% Functions
    function changepars(source, ~)
        t = source.Tag;
        switch t
            case 'input'
                nv=source.Value;
                if isempty(nv)
                    nv=I_input(nSelect,1);
                end
                I_input(nSelect,1)=nv;
                set(inFieldDisp{1},'string',num2str(nv,'%.2f'));
            case 'speed'
                nv=round(source.Value);
                if isempty(nv)
                    nv=nSpeed;
                end
                nSpeed=nv;
                set(source,'string',num2str(nv));
            case 'ac_amp'
                nv=source.Value;
                if isempty(nv)
                    nv=ac_amp(nSelect,1);
                end
                ac_amp(nSelect,1)=nv;
                set(inFieldDisp{2},'string',num2str(nv,'%.2f'));
            case 'exp_ad'
                nv=source.Value;
                if isempty(nv)
                    nv=exp_ad(nSelect,1);
                end
                exp_ad(nSelect,1)=nv;
                set(inFieldDisp{4},'string',num2str(nv,'%.2f'));
            case 'freq'
                nv=source.Value;
                if isempty(nv)
                    nv=freq(nSelect,1);
                end
                freq(nSelect,1)=nv;
                set(inFieldDisp{3},'string',num2str(nv,'%.2f'));
            case 'alpha'
                nv=source.Value;
                if isempty(nv)
                    nv=a(nSelect,1);
                end
                a(nSelect,:)=nv;
            case 'burst'
                nv=source.Value;
                burstFlag(nSelect,:) = nv;
                burstIter(nSelect,:) = 50;
            case 'Pulse'
                pulseFlag(nSelect,1) = 1;
            case 'neurons'
                val = source.Value;
                if val<=nNeuron
                    nSelect = val;
                    nDisp = val;
                else
                    nSelect = 1:nNeuron;
                end
            case 'pauseflag'
                pauseFlag = source.Value;
            case 'reset'
                I_input(:,1) = 0;
                ac_amp(:,1) = 0;
                freq(:,1) = 5;
                exp_ad(:,1) = 0;
                burstFlag(:,1) = 0;
                a(:,1)=3;
        end
        displayParams;
    end
    function spiking_modes(source,~)
        t = source.Tag;
        switch t
            case 'Tonic Spiking'
                I_input(nSelect,1) = 0.1;
                ac_amp(nSelect,1) = 0;
                freq(nSelect,1) = 5;
                exp_ad(nSelect,1) = 0;
                burstFlag(nSelect,1) = 0;
            case 'Tonal Excitation'
                I_input(nSelect,1) = 0;
                ac_amp(nSelect,1) = 0.1;
                freq(nSelect,1) = 5;
                exp_ad(nSelect,1) = 0;
                burstFlag(nSelect,1) = 0;
            case 'Bursting'
                burstFlag(nSelect,1) = 1;
                I_input(nSelect,1) = 0.1;
                ac_amp(nSelect,1) = 0;
                freq(nSelect,1) = 5;
                exp_ad(nSelect,1) = 0;
            case 'Integration'
                I_input(nSelect,1) = -0.02;
                ac_amp(nSelect,1) = 0;
                freq(nSelect,1) = 5;
                exp_ad(nSelect,1) = 9;
                burstFlag(nSelect,1) = 0;
            case 'Adaptation'
                I_input(nSelect,1) = 0.2;
                ac_amp(nSelect,1) = 0;
                freq(nSelect,1) = 5;
                exp_ad(nSelect,1) = 0.99;
                burstFlag(nSelect,1) = 0;
                
        end
        displayParams;
    end
    function displayParams        
        set(inFieldDisp{1},'string',num2str(I_input(nDisp)));
        set(inFieldDisp{2},'string',num2str(ac_amp(nDisp)));
        set(inFieldDisp{3},'string',num2str(freq(nDisp)));
        set(inFieldDisp{4},'string',num2str(exp_ad(nDisp)));
        set(inFieldDisp{5},'string',num2str(a(nDisp)));        
        set(inField{1},'Value',I_input(nDisp));
        set(inField{2},'Value',ac_amp(nDisp));
        set(inField{3},'Value',freq(nDisp));
        set(inField{4},'Value',exp_ad(nDisp));
        set(inField{5},'Value',a(nDisp));
        set(inField{6},'Value',burstFlag(nDisp,1));
    end
    function changeQ(source,~)
        val = source.Value;
        switch val
            case 1
                Q = eye(nNeuron, nNeuron);
            case 2
                Q = full(0.1*sprandn(nNeuron, nNeuron,0.1));
            case 3
                Q = 1*(rand(nNeuron, nNeuron)-0.5);
                for nf = 1:nNeuron
                    Q(nf,nf) = 0;
                end
            case 4
                nv = evalin('base','Q');
                disp(nv)
                if size(nv,1) == size(nv,2) && size(nv,2) == nNeuron
                    Q = nv;
                else
                    disp('Q is not a square matrix of size = number of neurons')
                    disp(size(nv,1))
                    disp(size(nv,2))
                end
                for nf = 1:nNeuron
                    Q(nf,nf) = 0;
                end
        end
        conn_im.CData = Q;
    end
end