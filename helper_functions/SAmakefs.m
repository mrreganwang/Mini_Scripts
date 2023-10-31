%% Import Steph's Data Get File; take out pulse; make into matrix of traces
function [fsmat,filename]=SAmake_fs(gain1,F,I,d)
filename=F{I};
source_files=[d,filename);
fileid=fopen(source_files);
fs=fread(fileid, 'int16');%
fs(1:512)=[];
for n=1:(length(fs)/14848)
    tmpmat(n,:)=fs((n-1)*14848+1:n*14848-10); % last 10 vals junk
end
fsmat=tmpmat(:,1+118+320*2:end); % take out pulse
fsmat=fsmat-mean(fsmat,2); % zero sweeps
Ns=(length(fs)/14848); % N sweeps
end