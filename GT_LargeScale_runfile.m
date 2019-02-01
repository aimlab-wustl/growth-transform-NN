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


% Creting a network with 300 neurons in the input layer and 700 neurons in
% the second layer
nNeuron = 1000;
l1Neuron = 300;
Q = createNet(nNeuron,l1Neuron,1);
%Graph functions are used to generate X and Y coordinates from the
%connectivity matrix. The graphs can be examined by commenting out close
%fig commands.
G1 = digraph(Q(1:l1Neuron,1:l1Neuron));
fig1 = figure;
hGraph = plot(G1,'Layout','force');
Xc1 = hGraph.XData;
Yc1 = hGraph.YData;
close(fig1)
G1 = digraph(Q(1+l1Neuron:end,1+l1Neuron:end));
fig2 = figure;
hGraph = plot(G1,'Layout','force');
Xc2 = hGraph.XData;
Yc2 = hGraph.YData;
close(fig2)
X = [Xc1 Xc2+-min(Xc2)+max(Xc1)+2];
Y = [Yc1 Yc2];
GT_LargeScaleGUI(Q,X,Y,l1Neuron)
