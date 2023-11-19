# Mini_Scripts

## Description
  Source code for detecting signals from a noisy trace. The method utilizes a feedforward neural network trained with the linear summation of noise and signal on a binary classification task, with targets *signal* and *noise*. Then, using a sliding window technique thrugh the trace of interest, a confidence value (CV) curve is produced where each point, x, on the CV curve gives the confidence level of the x-th window being *signal*. Locations that are deemed by the CV curve as *signal* are then smoothed out from the trace of interest, and the process is repeated until no new signals are found.
  
  The program was initially designed to detect minature excitatory postsynaptic currents (mEPSCs) from a noisy electrophisiological recording which adopts the shape of a pearson distribution. However, this method can be generalized to detect signals of any shape by training the neural network. More information can be found in the following paper (link to come) 

## Toolbox
  Here are the MATLAB toolboxes required for running the program
  
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox

  Instructions on installing toolboxes can be found [here](https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html)

## Usage
  We have provided three MATLAB scripts:
  
  1. `ROC.m` generates the ROC curve given a trained neural network and a noise trace. We have provided sampled trained network and noise trace in the **pre_trained_networks** and **noise_files** folders. The script is pre-set to use a trained network from the folder, you can change this by changing the path in the script to a trained network of your choice. When executing the script, a window is prompted for you to select a noise trace, you can use the default one we have provided in the folder
     
  2. `detect_peaks.m` is for detecting signals from a trace of your choice. When executing the sript, a window is prompted for you to select a trace you want to detect the peaks from. It is preset to use a trained network and noise trace provided. You can train your own network using `train_net.m`.
  
  3. `train_net.m` trains a neural network given a noise trace and a canonical shape of the signal of your choice. The signal must be a n-by-1 array where n is the number of data points that make up the signal peak.

  4. `create_noise.m` creates the noise trace from a signal file of your choice. When executing the script figures containing sections of your trace are prompted, select the region you believe contains **noise only** by left clicking on the start of the region then right click on the end of the region. Continue until you have adequate length of noise trace then press enter to finish
