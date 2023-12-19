%% create the noise file. 

% set the correct working directory (same directory this file is located in) 
wd = matlab.desktop.editor.getActive;
wd = fileparts(wd.Filename);
cd(wd);
% addpath 'helper_functions';
[fileName,filePath] = uigetfile('*.abf');

% load the data 
AS = abfload([filePath, fileName]);  % select the data you want to run
AS = AS - movmean(AS, 10000);  % take out slow data deviations

% create the noise trace. When each figure window is prompted, select the 
% region you believe contains only noise by left clicking on the start of
% the region then right click on the end of the region. Continue until you 
% have adequate length of noise trace then press enter to finish
noiseFileName = [fileName, '_noise.mat'];
[ASnoise] = MakeNoiseFromFile(AS);
ASnoise = ASnoise - movmean(ASnoise, 10000);% take out slow data deviations

%% save the noise trace created
save([wd, '/data_to_test/', noiseFileName], "ASnoise");