function [mininew]=MakeMiniMatK(D,Kset)
%                     MakeMiniMatK(D,A,mini);

load fastmini mini

mini=reshape(mini,[length(mini), 1]);
V=160; 
pk=60; 
strt=53;
%make D scaled and stretched clean minis
DD=D;
AMP=zeros(1,DD);
Fused=zeros(1,DD);
mininew=zeros(160,DD);
tic
for n=1:D
    FF=.7+.3*rand; %expands mini by 1/FF
    minits=timeseries(mini,1:V);
    dataRange1=FF*(1:(1/FF)*V); dataRange1(1)=1;
    minirs=resample(minits,dataRange1);%resample mini at FF intervals
    dataRange=floor((1/FF)*pk)-pk+1:V+floor((1/FF)*pk)-pk;%keep mini peak at pk
    minitmp=minirs.Data(dataRange); 
    AMP(n)=Kset;%(pearsrnd(10,3,.5,3,1,1).^2)/100; %ADD 3.5? makes mini with mean 1
    Fused(n)=FF;% keep track of mini stretch
    mininew(:,n)=AMP(n)*minitmp; %scale mini by rand from skewed dist
end
% toc
% clear tmp tmp2
% % tmp=zeros(V,D);%tmp2=zeros(20000/DD);
% tmp2=mininew;
% tic
% for n=1:D/1000-1
% tmp2=[tmp2,mininew];
% end
% tmp2=repmat(mininew,1,D/DD);
toc
% mininew=reshape(tmp2,[D,size(tmp2,1)]);

end