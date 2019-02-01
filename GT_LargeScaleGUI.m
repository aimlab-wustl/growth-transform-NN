function GT_LargeScaleGUI(Q, Xc, Yc,l1Neuron)
% Creates a GUI for simulating a large network of neurons. In this case half
% of the neurons in the first layer get random inputs during stimulation.
% Second layer gets inputs (excitatory/inhibitory) from the first layer.
% There is no feedback to the first layer. Rasters and PCA trajectories are
% plotted for visualization.

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

C =[0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560 ;0.4660    0.6740    0.1880];
% colors for plots
nNeuron = size(Q,1);
arrowCord = [max(Xc(1:l1Neuron)) + 0.1, min(Xc(l1Neuron+1:end))];
figNumber = figure;
set(figNumber,'NumberTitle','off',...
    'Name','Growth Transform Neuron Network Model',...
    'Units','normalized','toolbar','figure',...
    'Position',[0.05 0.1 0.9 0.8]);
h1 = axes('Position',[0.05 0.5 0.9 0.4],'Color',[0 0 0]);
% colormap(h1,'gray')
colormap(h1,[zeros(10,1) linspace(0,1,10)' zeros(10,1)])
caxis([-3 5])
hold on
hGraph = scatter(h1,Xc,Yc,[],zeros(nNeuron,1),'filled');
% text(arrowCord(1)+0.5,0,'\rightarrow','FontSize',14,'FontWeight','Bold',...
%     'Color',[1 1 1])
text(mean(Xc(1:l1Neuron)),max(Yc(1:l1Neuron))+2,'Level 1 neurons',...
    'HorizontalAlignment','center','FontSize',14,'FontWeight','Bold',...
    'Color',[0 1 0])
text(mean(Xc(l1Neuron+1:end)),max(Yc(1:l1Neuron))+2,'Level 2 neurons',...
    'HorizontalAlignment','center','FontSize',14,'FontWeight','Bold',...
    'Color',[0 1 0])
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);


colorMap = repmat(linspace(0,0.7,30)',1,3);
colorMap = [colorMap;[1 1 1];colorMap(end:-1:1,:)];
colorMap(1:30,1) = 1;
colorMap(32:61,3) = 1;
hConn=axes('Position',[0.05 0.15 0.25 0.25]);
IM = 5*imgaussfilt(Q',2);
conn_im = imagesc(hConn,IM);
colormap(hConn,colorMap)
set(gca,'xtick',[],'ytick',[])
xlabel('Post-synaptic')
ylabel('Pre-synaptic')
title('Connectivity Matrix')
caxis([-1 1])
caxis manual
colorbar

hPCA=axes('Position',[0.65 0.15 0.25 0.25]);
title('PCA trajectories')
view(3)
hRaster=axes('Position',[0.35 0.15 0.25 0.25]);
title('Raster plot')
xlabel('Time (ms)')
% Pause button
runFlag = 0;
runButton = uicontrol('Style', 'pushbutton','String','Run',...
    'Min',0,'Max',1,'Value',0, 'Units','normalized', ...
    'Position',[0.55 0.02 0.1 0.05],...
    'tag','runFlag','Callback',@changepars);
soundButton = uicontrol('Style', 'pushbutton','String','Sonify',...
    'Min',0,'Max',1,'Value',0, 'Units','normalized', ...
    'Position',[0.85 0.02 0.1 0.05],...
    'tag','runFlag','Callback',@playSound);
uicontrol('Style', 'pushbutton','String','PCA Plot',...
    'Min',0,'Max',1,'Value',0, 'Units','normalized', ...
    'Position',[0.7 0.02 0.1 0.05],...
    'tag','runFlag','Callback',@pcaButton);
%%
nIter = 10;
% rng(1234,'twister');

iInput0 = zeros(nNeuron,1)-0.01;
iInput0(l1Neuron+1:end) = -0.01;

nTimes = 30;
tStim = [10 20];
nStims = 3;
D = zeros(nNeuron,nIter*nTimes,nStims);
SPK = zeros(nNeuron,nIter*nTimes-1,nStims);
axes(h1)
o1 = 0;
SonWt = zeros(2,nNeuron-l1Neuron);
SonWt(1,1:2:end) = 1;
SonWt(2,2:2:end) = 1;
% SonWt(1,1:50) = abs(randn(1,50));
% SonWt(2,51:100) = abs(randn(1,50));
while ishandle(figNumber)
    if runFlag
        o1 = mod(o1,3)+1;
        runbutton.Enable = 'off';
        axes(h1)
        Y = [];
        S = [];
        yInit = zeros(nNeuron,1);
        iInput = iInput0;
        inputIdx = randperm(l1Neuron,l1Neuron/2);
        iInput(inputIdx,:) = 0.008*randn(l1Neuron/2,1);
        txt1 = text(min(Xc(1:l1Neuron))-0.5,mean(Yc(1:l1Neuron)),'',...
            'HorizontalAlignment','right','FontSize',14,'FontWeight','Bold');
        for t1 = 1:nTimes
            if t1>tStim(1) && t1<tStim(2)
                I = iInput;
                txt1.String = 'Stimulus';
                txt1.Color = C(o1,:);
            else
                I = iInput0;
                txt1.String = '';
            end
            [yTarget,  yInit, spk] = GT_LargeScaleFun(Q,I,nIter+1,yInit);
            Y = cat(2,Y,yTarget(:,2:end));
            S = cat(2,S,spk);
            hGraph.CData = sum(spk,2);
            drawnow
            pause(0.1)
        end
        runFlag = 0;
        runButton.Value = 0;
        runbutton.Enable = 'on';
        D(:,:,o1) = Y;
        SPK(:,:,o1) = o1*(logical(S(:,1:end-1))|logical(S(:,2:end)));
        PCA_data = D(l1Neuron+1:end,tStim(1)*nIter:tStim(2)*nIter,:);
        Son_data = Y(l1Neuron+1:end,tStim(1)*nIter:tStim(2)*nIter,:);
        title('PCA of level 2 neurons')
        
        axes(hRaster)
        SPK_plot = [SPK(:,:,1) SPK(:,:,2) SPK(:,:,3)];
        imagesc(SPK_plot)
 %       colormap(hRaster,[1 1 1; C(1:max(SPK_plot(:)),:)])
        title('Raster plot')
        xlabel('Time (ms)')
%         SPK_plot(SPK,C,hRaster)
        
    end
    pause(0.1)
end
% Functions
    function changepars(source, ~)
        t = source.Tag;
        switch t
            case 'runFlag'
                runFlag = source.Value;
                
        end
    end
    function playSound(~,~)
        soundVec = SonWt*Son_data;
        sound(resample(10*soundVec',10,1),1000)
    end
    function pcaButton(~,~)
        PCA_plot(PCA_data,C,3,hPCA)
    end
end
function PCA_plot(data,C,nFilt,figNumber)
nCh = size(data,1);
X = reshape(data,nCh,[]);

X = X';

k = size(data,3);
[~,score,latent] = pca(X);
D2 = round(10000*latent/sum(latent))/100;
X_PCA2 = score(:,1:3);
X_PCA2 = filtfilt(ones(1,nFilt),nFilt,X_PCA2);
data_PCA = permute(reshape(X_PCA2',3,size(data,2)/1,[]),[2 1 3]);
axes(figNumber)
cla(figNumber)
hold on
view(3)
for o1 = 1:k
    plot3(data_PCA(:,1,o1),data_PCA(:,2,o1),data_PCA(:,3,o1),...
        'Color',C(o1,:),'linewidth',1.5,'Marker','o','MarkerFaceColor',C(o1,:),'MarkerEdgeColor',C(o1,:));
end
grid on
box off
xlabel(['PC1 (' num2str(D2(1)) '%)'],'FontName','Arial','FontSize',10,'FontWeight','Bold')
ylabel(['PC2 (' num2str(D2(2)) '%)'],'FontName','Arial','FontSize',10,'FontWeight','Bold')
zlabel(['PC3 (' num2str(D2(3)) '%)'],'FontName','Arial','FontSize',10,'FontWeight','Bold')
end
function SPK_plot(S,C,figNumber)
axes(figNumber)
cla(figNumber)
hold on
nT = size(S,2);
nN = size(S,1);
for o1 = 1:3
    tmp = S(:,:,o1);
    xc = nT*(o1-1);
    for n1 = 1:nN
        for t1 = 1:nT
            if tmp(n1,t1)
                plot(xc+[t1 t1],[n1-0.4 n1+0.4],'Color',C(o1,:))
            end
        end
    end
end
end