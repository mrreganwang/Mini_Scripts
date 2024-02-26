%% 1. prompt an window to select the file to be tested 

% set the correct working directory (same directory this file is located in) 
wd = matlab.desktop.editor.getActive;
wd = fileparts(wd.Filename);
cd(wd);
addpath 'helper_functions';
[fileName,filePath] = uigetfile('*.abf');

% load the data 
[trace, si, h] = abfload([filePath, fileName]);  % select the data you want to run 
sampling_freq = 1 / si * 10^6;
trace = trace - movmean(trace, sampling_freq);  % take out slow data deviations

if sampling_freq == 10000
    load([wd,'/pre_trained_networks/trained_net_10k.mat']);     % load trainedNetN
elseif sampling_freq == 2000
    load([wd,'/pre_trained_networks/trained_net_2k.mat']);     % load trainedNetN
end


% create the noise file. choose random interval of size 10 from the trace. 
% concatenate 1000 times to achieve full noise trace
randomS = floor(rand() * (size(trace, 1) - 11));
randInterval = trace(randomS: randomS + 9);
randInterval = randInterval - mean(randInterval);

ASnoise = repmat(randInterval', 1, 1000);
noiseSD = std(ASnoise);

%% 2. find peaks

% INITIALIZE
dataLength = size(trace, 1);       % total length of the trace 
smoothedTrace = trace;             % smoothed to obtain the new CV every iteration
originalTrace = smoothedTrace;  % for plotting purposes
TN = trainedNet;               % trained neural network
pks = 1;                        % store the found peaks in each iteration                                    
Pk_list = [];                   % store the found peaks across all iterations 
POS = 0;                        % parameter for plotting figures
SL = floor(300 * (sampling_freq / 10000)); % sweeping length. number of data points for one peak
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
    CV = movmean(CV, 50);
 
    % FIND CV PEAKS
    [~,pks,~,~]=findpeaks...
        (CV,'MinPeakDistance', mpd,'MinPeakProminence',0.2,'MinPeakWidth',mpw);
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
    disp([num2str(numel(pks)), ' peak(s) found in iteration ', num2str(i)]);
    disp(pks);
    i = i + 1;
    if i > 11
        break;
    end
    numPeakTable = [numPeakTable; [numel(pks), numel(Pk_list)]];
end

%% 3. save the result

uisave({'numPeakTable', 'trace', 'Pk_list'});
