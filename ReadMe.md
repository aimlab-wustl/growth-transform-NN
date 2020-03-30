# Growth Transform neural network

This repository contains the codes for the paper 'A Spiking Neuron and Population Model based on the Growth Transform Dynamical System', bioRxiv (2019), 523944. Please see the PDF pre-print of the paper that describes the Growth Transform neuron model.

The file descriptions are given below:

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) **Single-neuron_dynamics.m:** Generates plots for single-neuron response characteristics, e.g. tonic spiking, bursting, spike-frequency adaptation, etc. 

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) **GT_NeuronModel.m:** This function creates a GUI for simulating Neuron Models based on growth transforms [1, 2].
The GT Neron model is capable of emulating a number of dynamics which are observed in 
biological systems, such as tonal spiking in presence of constant input, excitation by sinusoidally varying input stimuli, leaky 
integrator, and as bursting neurons. All of these dynamics can be observed in the GUI by clicking the relevant buttons. Users 
have the ability to change the input current (DC and AC components), neuronal adaptation, bursting behavior, and the hyperpolarization
parameter. This can be done for an individual neuron or for all neurons simultaneously. Clicking reset resets all inputs to their 
default values.
The GUI can also simulate a network of neurons with inhibitory and excitatory connections. There is an option for 
generating sparse random connection matrix in the GUI. The user can also input their own connection matrix 'Q' by defining it in 
the workspace and then importing it into the GUI using the Connectivity Matrix Dropdown Menu. The matrix should be called Q and 
it should be a square matrix of size = number of neurons. This GUI works best for 1-40 neurons.

![#1589F0](https://placehold.it/15/1589F0/000000?text=+) **GT_LargeScale_runfile.m:** Runs a large network of neurons and creates visualization for the spiking activity and population trajectories.



Copyright (c) [2018] Washington University  in St. Louis
Created by: [Darshit Mehta, Ahana Gangopadhyay, Kenji Aono, Shantanu Chakrabartty]

Citations for this tool are: 

1. Gangopadhyay, A., and Chakrabartty, S. Spiking, bursting, and population dynamics in a network of growth transform neurons. IEEE Transactions on Neural Networks and Learning Systems, 2017.
2. Gangopadhyay, A., Chatterjee, O., and Chakrabartty, S. Extended polynomial growth transforms for design and training of generalized support vector machines. IEEE Transactions on Neural Networks and Learning Systems, 2017.
3. Gangopadhyay, A., Aono, K. Mehta, D., and Chakrabartty, S. A Coupled Network of Growth Transform Neurons for Spike-Encoded Auditory Feature Extraction. Preprint at bioRxiv, 2018
4. Gangopadhyay, A., Mehta, D., and Chakrabartty, S. A Spiking Neuron and Population Model based on the Growth Transform Dynamical System. Preprint at bioRxiv, 2019



Washington University hereby grants to you a non-transferable, non-exclusive, royalty-free, non-commercial, research license to use and copy the computer 
code provided here.  You agree to include this license and the above copyright notice in all copies of the Software.  The Software may not be distributed, 
shared, or transferred to any third party.  This license does not grant any rights or licenses to any other patents, copyrights, or other 
forms of intellectual property owned or controlled by Washington University.  If interested in obtaining a commercial license, please contact 
Washington University's Office of Technology Management (otm@dom.wustl.edu).

YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND 
IS PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR 
FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  IN NO EVENT SHALL 
THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR 
IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF 
SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES. YOU ALSO AGREE THAT THIS SOFTWARE WILL NOT BE USED FOR CLINICAL PURPOSES.
