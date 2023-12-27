%% load the canonical peak shape

% set the correct working directory (same directory this file is located in) 
wd = matlab.desktop.editor.getActive;
wd = fileparts(wd.Filename);
cd(wd);
addpath 'helper_functions';

% prompt window to get the peak shape file. Should be a dx1 array stored in
% a .mat file. Check out 'FastMiniBo.mat' for the format of the signal
[peakFileName, peakFilePath] = uigetfile('*.mat');
load([peakFilePath, peakFileName]);
signal = mini;     % set signal to be the variable name that stores the canonical peak shape

load([peakFilePath, peakFileName]);


%% load the noise file

% prompt window to get the noise file
[noiseFileName, noiseFilePath] = uigetfile('*.mat');
load([noiseFilePath, noiseFileName]);
ASnoise = ASnoise * 2.5;                            % hyperparameter, scale the noise 

%% create training data

n = 10000;                          % training set contains n positives and n negatives
SL = size(signal, 1);               % sweeping length. number of data points for one peak
amp = 5;                            % hyperparameter, mean peak amplitude of the signals in the training set     

% [x] = buildTrainSet(n, SL, amp, 1, signal);     % create n signals with varying width and amplitude
[x] = MakeMiniMatFS(n, SL, amp, 1, signal);

% [X, Y]=MakeTrainMat(-x,reshape(ASnoise, [1 numel(ASnoise) ]), n, SL);% mini mean=StNset
[X, Y]=MakeTrainMat(-x,reshape(ASnoise, [1 numel(ASnoise) ]), n, SL);% mini mean=StNset

net0 = patternnet([200 100 100]);   % initialize neural net  
netN = configure(net0,X,Y);
[trainedNet,~] = train(netN,X,Y);

%%
uisave({'trainedNet'});

