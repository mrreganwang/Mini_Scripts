function [mnTr, mnTa]=RMmytrainset1(mininew,circMat, D, StNset)
% StN -> set signal to noise of mini/noise
%D = # minis generated, D=#noise traces
% mini has peak at pk; length=V; 
n1=mean(std(circMat,1,1));

load fastmini mini
V=160; 
pk=60; 
strt=53;
MM=1;
StN=.1;
while StN<StNset && StN>0 % no idea why i put this in
MM=MM+.1;
StNs=.1*MM;
StN=(10/n1)*StNs;
end

%make D noisy aligned mini of variable amplitude (skewed amplitude distribution)
clear tmp tmp1 miniTr noiseTr mnTr mnTa
tmp1=circMat;%tmp1=tmp1*0;
miniTr=zeros(D,160);

for n=1:D
    tmp=round(rand*(size(tmp1,1)-1)+1);
    miniTr(n,:)=tmp1(tmp,:)+(mininew(n,:))*StNs; %add real noise scale mini by noise
end
miniTr=miniTr'; % mini sweeps as columns

%generate D noise sweeps
noiseTr=zeros(D,160);
for n=1:D
    tmp2=round(rand*(size(tmp1,1)-1)+1);
    noiseTr(n,:)=tmp1(tmp2,:); %add real noise StN*
end
noiseTr=noiseTr'; 

% cat mini and noise sweeps into training matrix
mnTr=cat(2,miniTr,noiseTr); 

% target matrix
mnTa=[ones(1,D),zeros(1,D);zeros(1,D),ones(1,D)];

% % Get tmp5=std(mnTr(:,D+1:end),1,1); 
% StNTR=(mean(AMP)*StNs)/mean(std(mnTr(:,D+1:end),1,1),2); %measured S/N for training set
% StNR=(mean(AMP)*StNs)/n1; %measured S/N for Real
% AMPval=mean(AMP);
% 'here'
end