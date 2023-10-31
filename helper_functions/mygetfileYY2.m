%% Get File
function [circMatallsweepsdd,allsweeps_NBQX_added,AMPval,Fusedval,SNstdev]=mygetfileYY3(allsweeps_NBQX)
%%
load fastmini mini
V=160; 
pk=60; 
D=floor((size(allsweeps_NBQX,1)/160))*size(allsweeps_NBQX,2);
D1=floor((size(allsweeps_NBQX,1)/160));

%make D scaled and stretched clean minis for validation
AMPval=zeros(1,D);
Fusedval=zeros(1,D);
mininew1=zeros(160,D);
for n=1:D
    F=.3+.7*rand; %expands mini by 1/F
    minits=timeseries(mini,1:V);
    dataRange1=F*(1:(1/F)*V); dataRange1(1)=1;
    minirs=resample(minits,dataRange1);%resample mini at F intervals
    dataRange=floor((1/F)*pk)-pk+1:V+floor((1/F)*pk)-pk;%keep mini peak at pk
    minitmp=minirs.Data(dataRange); 
    AMPval(n)=(pearsrnd(5,3,.5,3,1,1).^2)/10; %ADD 3.5?
    Fusedval(n)=F;% keep track of mini stretch
    mininew1(:,n)=AMPval(n)*minitmp; %scale mini by rand from skewed dist
end

%make D noisy aligned mini of variable amplitude (skewed amplitude distribution)
tmp1=allsweeps_NBQX;
SNstdev=zeros((D1-1),size(allsweeps_NBQX,2));

for n=1:size(allsweeps_NBQX,2)
    for nn=0:(D1-2) %The size of this seems to matter based on the size of the initial array
    SNstdev(nn+1,n)=std(tmp1(161+(nn*160):320+(nn*160),n));
    tmp1(161+(nn*160):320+(nn*160),n)=tmp1(161+(nn*160):320+(nn*160),n)+mininew1(:,((nn+1)+((n-1)*D1)));
    end
end
allsweeps_NBQX_added=tmp1;
SNstdev=SNstdev(:)';
%%
circMatallsweepsdd=[];%zeros(160, (size(allsweeps,2)*round(run1num)));
circMat=zeros(160,(size(tmp1,1)-159));
for nnn=0:(size(tmp1,2)-1)
    for n=1:(size(tmp1,1)-159)
        circMat(:,n)= tmp1((n:n+159), (nnn+1));
    end
circMatallsweepsdd(1:160, (nnn*size(circMat,2)+1):(nnn*size(circMat,2)+size(circMat,2)))=circMat;
end

end