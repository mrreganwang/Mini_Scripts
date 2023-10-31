function [mininew, AMPval, AMP]=MakeMiniWN(circMat, D, StNset)
% StN -> set signal to noise of mini/noise
%D = # minis generated, D=#noise traces
% mini has peak at pk; length=V; 
n1=mean(std(circMat,1,1));

load fastmini mini
V=160; 
pk=60; 
strt=53;
% MM=1;
% StN=.1;
% while StN<StNset && StN>0
% MM=MM+.1;
% StNs=.1*MM;
% StN=(10/n1)*StNs;
% end

%make D scaled and stretched clean minis
AMP=zeros(1,D);
Fused=zeros(1,D);
mininew=zeros(160,D);
for n=1:D 
    FF=.7+.3*rand; %expands mini by 1/FF
    minits=timeseries(mini,1:V);
    dataRange1=FF*(1:(1/FF)*V); dataRange1(1)=1;
    minirs=resample(minits,dataRange1);%resample mini at FF intervals
    dataRange=floor((1/FF)*pk)-pk+1:V+floor((1/FF)*pk)-pk;%keep mini peak at pk
    minitmp=minirs.Data(dataRange); 
    AMP(n)=(pearsrnd(10,3,.5,3,1,1).^2)/10; %ADD 3.5?
    Fused(n)=FF;% keep track of mini stretch
    mininew(:,n)=AMP(n)*minitmp; %scale mini by rand from skewed dist
end
mininew=mininew';
AMPval=mean(AMP);
% 
% %make D noisy aligned mini of variable amplitude (skewed amplitude distribution)
% tmp1=circMat;
% miniTr=zeros(D,160);
% clear tmp
% for n=1:D
%     tmp=round(rand*(size(tmp1,1)-1)+1);
%     miniTr(n,:)=tmp1(tmp,:)+(mininew(n,:)*StNs); %add real noise 
% end
% miniTr=miniTr'; % mini sweeps as columns
% 
% %generate D noise sweeps
% noiseTr=zeros(D,160);
% for n=1:D
%     tmp2=round(rand*(size(tmp1,1)-1)+1);
%     noiseTr(n,:)=tmp1(tmp2,:); %add real noise StN*
% end
% noiseTr=noiseTr'; 
% 
% % cat mini and noise sweeps into training matrix
% mnTr=cat(2,miniTr,noiseTr); 
% 
% % target matrix
% mnTa=[ones(1,D),zeros(1,D);zeros(1,D),ones(1,D)];
% 
% % Get tmp5=std(mnTr(:,D+1:end),1,1); 
% StNTR=(mean(AMP)*StNs)/mean(std(mnTr(:,D+1:end),1,1),2); %measured S/N for training set
% StNR=(mean(AMP)*StNs)/n1; %measured S/N for Real
% AMPval=mean(AMP);
end