%% load the neural net

wd = matlab.desktop.editor.getActive;
wd = fileparts(wd.Filename);
cd(wd);
% addpath 'helper_functions';

load([wd,'/pre_trained_networks/trained_net.mat']);     % load trainedNetN

%% select the noise file to be used for testing
                                
noiseScale = 2.5;                                                                                     
[fileName,filePath] = uigetfile('*.mat');
load([filePath, fileName]);
ASnoise = ASnoise * noiseScale;     
noiseSD = std(ASnoise);

%% create the ROC
% set the parameters
D = 1000;   % number of synthetic peaks in the test trace
mode1 = 1;  % 1 if peaks with constant amplitude. 2 if pearson distributed. 3 if uniform between average-0.5 to average+0.5. 
mode2 = 2;  % 1 if evenly spaced minis 2 if variable spacing
SL = 300;               % sweeping length. number of data points for one peak
mpw = 3;                % minimum peak width
mpd = 4;                % minimum peak distance 
traceLength = 700000;   % length of trace to create ROC curve
s = 3;                  % scale for ROC bench mini (to create ROC for minis with avg amplitude of s, 2s, 3s, ...)
numAmp = 1;             % number of amplitudes to test
timeBeAf = 10*2;        % time window about peak to be considered a true positive
% thresholdList = [0.1,0.2,0.4,0.6,0.7,0.75,0.8,0.85,0.9,0.95,0.99,0.999]; % thresholds to create ROC
% smoothingList = [1,3,10,20,30,40,50,60,70,80];                           % smoothings to create ROC
smoothingList = [5];
thresholdList = [0.975];

% pad the noise trace to a longer trace
T = round((.5*traceLength/numel(ASnoise)),0);
ASnoise2 = ASnoise;
for n = 1:T
    ASnoise2 = [ASnoise2;ASnoise(100:end)];
end

% initialize cell arays
TPR = cell(1,numAmp);
FDR = cell(1,numAmp);
FNR = cell(1,numAmp);
TPLOC = cell(1,numAmp);
FPLOC = cell(1,numAmp);
FNLOC = cell(1,numAmp);
SIGNALTRACE = cell(1,numAmp);   % trace with no noise
NOISETRACE = cell(1,numAmp);    % noise trace
TRUELOC = cell(1,numAmp);       % peak locations
FF = cell(1,numAmp);            % store figures from each iteration

for M = 1:numAmp
    Ms = s * (M);
    disp('running genROC loop')

    [TPRt,FDRt,FNRt,TPLoc,FPLoc,FNLoc,signalTrace,noiseTrace,trueLoc,FFt]=genROC5b(smoothingList,thresholdList,Ms,mpw,mpd,SL,ASnoise2,trainedNetN,D,mode1,mode2,timeBeAf); % 1 for CV smoothing

    TPR{M} = TPRt;
    FDR{M} = FDRt;
    FNR{M} = FNRt;
    TPLOC{M} = TPLoc;
    FPLOC{M} = FPLoc;
    FNLOC{M} = FNLoc;
    SIGNALTRACE{M} = signalTrace;
    NOISETRACE{M} = noiseTrace;
    TRUELOC{M} = trueLoc;
    FF{M} = FFt;
end

%% save results

uisave({'TPR', 'FDR', 'FNR', 'TPLOC', 'FPLOC', 'FNLOC', 'SIGNALTRACE', 'NOISETRACE', 'TRUELOC', 'FF'});
