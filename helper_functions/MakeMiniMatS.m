function [mininew]=MakeMiniMatS(D,SL,ftl)
% StN -> set signal to noise of mini/noise
%D = # minis generated, D=#noise traces
% mini has peak at pk; length=V; 
% n1=std(allsweeps);
load(ftl, 'mini')
% mini=smini2;
DD=D;
V=SL; 
[~,pk]=max(mini(:)); 
strt=53;
AMP=zeros(1,DD);
Fused=zeros(1,DD);
mininew=zeros(300,DD);
for n=1:DD 
    FF=.6+.4*rand; %expands mini by 1/FF
    minits=timeseries(mini,1:V);
    dataRange1=FF*(1:(1/FF)*V); dataRange1(1)=1;
    minirs=resample(minits,dataRange1);%resample mini at FF intervals
    dataRange=floor((1/FF)*pk)-pk+1:V+floor((1/FF)*pk)-pk;%keep mini peak at pk
    minitmp=minirs.Data(dataRange); 
    AMP(n)=(pearsrnd(10,3,.5,3,1,1).^2)/100; %ADD 3.5? makes mini with mean 1
    Fused(n)=FF;% keep track of mini stretch
    mininew(:,n)=AMP(n)*minitmp; %scale mini by rand from skewed dist
end
clear tmp tmp2
% tmp=zeros(V,D);%tmp2=zeros(20000/DD);
tmp2=repmat(mininew,1,D/DD);
mininew=tmp2;
% mininew=reshape(tmp2,[D,size(tmp2,1)]);
% 
% mininew=mininew';
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