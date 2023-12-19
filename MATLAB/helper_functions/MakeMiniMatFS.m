function [mininewF, AMPval, AMP]=MakeMiniMatFS(D,SL,K,mode, signalShape)
% StN -> set signal to noise of mini/noise
%D = # minis generated, D=#noise traces
% mini has peak at pk; length=V; 
% n1=std(allsweeps);
tStart = tic;
% load FastMiniBo mini
% miniF=mini;
% load SlowMiniBo mini
% miniS=mini;
miniF = signalShape;
V=SL; 
[~,pk]=max(miniF(:));
% pk=60; 
strt=53;
% MM=1;
% StN=.1;
% while StN<StNset && StN>0
% MM=MM+.1;
% StNs=.1*MM;
% StN=(10/n1)*StNs;
% end

%make D scaled and stretched clean minis
mini=miniF;
% D=D/2;
AMP=zeros(1,D);
Fused=zeros(1,D);
mininew=zeros(SL,D);
for n=1:D 
    FF=.6+.4*rand; %expands mini by 1/FF
    minits=timeseries(mini,1:V);
    dataRange1=FF*(1:(1/FF)*V); dataRange1(1)=1;

    minirs=resample(minits,dataRange1);%resample mini at FF intervals
    size(dataRange1)
    size(minirs.Data)
    dataRange=floor((1/FF)*pk)-pk+1:V+floor((1/FF)*pk)-pk;%keep mini peak at pk
    pk
    FF
    floor((1/FF)*pk)-pk+1
    V+floor((1/FF)*pk)-pk
    minitmp=minirs.Data(dataRange); 

%     AMP(n)=3 + K*(pearsrnd(10,3,.5,3,1,1).^2)/100; %ADD 3.5? makes mini with mean 1
% AMP(n) = K*((pearsrnd(10,3,.5,3,1,1).^2)/100);
    if mode == 1
        AMP(n) = K;
    elseif mode == 2
        AMP(n)=K*(pearsrnd(10,3,.5,3,1,1).^2)/100; %ADD 3.5? makes mini with mean 1
    else 
        AMP(n)=unifrnd(K-0.5,K+0.5);
    end
    Fused(n)=FF;% keep track of mini stretch
    mininew(:,n)=AMP(n)*minitmp; %scale mini by rand from skewed dist
end
mininewF=mininew';%figure;for n=1:10;plot(mininewF(n,:));hold on; end
AMPval=mean(AMP);

% mini=miniS;
% AMP=zeros(1,D);
% Fused=zeros(1,D);
% mininew=zeros(SL,D);
% for n=1:D 
%     FF=.6+.4*rand; %expands mini by 1/FF
%     minits=timeseries(mini,1:V);
%     dataRange1=FF*(1:(1/FF)*V); dataRange1(1)=1;
%     minirs=resample(minits,dataRange1);%resample mini at FF intervals
%     dataRange=floor((1/FF)*pk)-pk+1:V+floor((1/FF)*pk)-pk;%keep mini peak at pk
%     minitmp=minirs.Data(dataRange); 
%     if mode == 1
%         AMP(n) = K;
%     elseif mode == 2 
%         AMP(n)=K*(pearsrnd(10,3,.5,3,1,1).^2)/100; %ADD 3.5? makes mini with mean 1
%     else
%         AMP(n)=unifrnd(K-0.5,K+0.5);
%     end
%     Fused(n)=FF;% keep track of mini stretch
%     mininew(:,n)=AMP(n)*minitmp; %scale mini by rand from skewed dist
% end
% mininewS=mininew';%figure;for n=1:10;plot(mininewF(n,:),'k');hold on; end
% mininew=[mininewF;mininewS];%
% % figure;plot(mean(mininew(1:100,:),1));hold on; plot(mean(mininew(end-100:end,:),1));
% AMPval=mean(AMP);
tEnd = toc(tStart);
['MakeMiniMatFS took: ', num2str(tEnd), ' seconds']
toc
end