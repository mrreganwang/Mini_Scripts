%% 1. prompt an window to select the file to be tested 

% set the correct working directory (same directory this file is located in) 
wd = matlab.desktop.editor.getActive;
wd = fileparts(wd.Filename);
cd(wd);
% addpath 'helper_functions';
[fileName,filePath] = uigetfile('*.abf');

load([wd,'/pre_trained_networks/trained_net.mat']);     % load trainedNetN

% load the data 
trace = abfload([filePath, fileName]);  % select the data you want to run
trace = trace - movmean(trace, 10000);  % take out slow data deviations

%% 2. load the noise file. 

load([wd, '/noise_files/2022_02_28_0000_noise.mat']);
ASnoise = ASnoise * 2.5;        % scale noise by 2.5
noiseSD = std(ASnoise);

%% 3. find peaks

% INITIALIZE
dataLength = size(trace, 1);       % total length of the trace 
smoothedTrace = trace;             % smoothed to obtain the new CV every iteration
originalTrace = smoothedTrace;  % for plotting purposes
TN = trainedNetN;               % trained neural network
pks = 1;                        % store the found peaks in each iteration                                    
Pk_list = [];                   % store the found peaks across all iterations 
POS = 0;                        % parameter for plotting figures
SL = 300;                       % sweeping length. number of data points for one peak
mpd = 4;                        % minimum peak distance   
mpw = 3;                        % minimum peak width
Nn = ASnoise;                   % value to set the found peaks to after each iteration
i = 1;                          % keep track of the number of iterations
numPeakTable = [];              % num peaks found in each iteration and total number of peaks found

while ~isempty(pks)
 
    % GENERATE CV
    CV = TN(MakeCircMatData(SL, smoothedTrace));
    CV = CV(1,:);
    CV(CV < 0) = 0;
    CV = movmean(CV, 10);
 
    % FIND CV PEAKS
    [~,pks,~,~]=findpeaks...
        (CV,'MinPeakDistance', mpd,'MinPeakProminence',.9,'MinPeakWidth',mpw);
    Pk_list = [Pk_list, pks]; % add peaks found to list
 
    % SUBTRACT DETECTED MINI FROM TRACE
    if ~isempty(pks)
        for n=1:numel(pks)
            inval = pks(n)-50+70:pks(n)+70+100;
            smoothedTrace(inval) = 0;   % zero out peak
            smoothedTrace(inval+100) = Nn(1:numel(inval)); % set trace after mini peak to noise values
        end
    end
 
    % PLOT CV ORIGINAL AND REMADE TRACE
    f = figure; 
    f.OuterPosition = [10, 900 - POS, 2200, 500];
    POS = POS + 300;
    plot(originalTrace(71:end));
    hold on; 
    plot(smoothedTrace(71:end));
    plot(10 * CV(1:end));
    scatter(Pk_list, originalTrace(Pk_list+71), 'filled');
    disp([num2str(numel(pks)), ' peaks are found in iteration ', num2str(i)]);
    i = i + 1;
    numPeakTable = [numPeakTable; [numel(pks), numel(Pk_list)]];
end

%% 4. save the result

uisave({'numPeakTable', 'trace', 'Pk_list'});
