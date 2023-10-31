%% Clear Previous Data
clear 
clc
%% Import Data MJM
gain1=20; %100 for Yifan 50 or 20 for MJM
selection1='Get Before NBQX File';
[circMatallsweeps]=mygetfilemjm(gain1, selection1);
selection2='Get After NBQX File';
[circMatallsweeps_NBQX]=mygetfilemjm(gain1, selection2);
clear gain1 selection1 selection2

%% Import Data SA
selection1='Get Before NBQX File';
[circMatallsweeps, allsweeps]=mygetfileSA(selection1);
selection2='Get After NBQX File';
[circMatallsweeps_NBQX, allsweeps_NBQX]=mygetfileSA(selection2);
clear selection1 selection2

%% Import Data YY
gain1=100;
selection1='Get Before NBQX File';
[circMatallsweeps]=mygetfileYY(gain1, selection1);
selection2='Get After NBQX File';
[circMatallsweeps_NBQX,allsweeps_NBQX]=mygetfileYYnbqx(gain1, selection2);
clear gain1 selection1 selection2

%% generate training set (or only minis)
N=4;
StNset=N;
D=20000; %number of training set minis
[mnTr, mnTa, AMPval]=mytrainset1(circMatallsweeps,D,StNset);
n3=mean(std(circMatallsweeps_NBQX));
StNNBQX=mean(AMPval)/n3; %measured S/N for NBQX

%Create Test Data
[circMatdd,AMPval,fusedval,SNstdev]=mygetfileYY2(allsweeps_NBQX);

% make net and train
net = patternnet([200 100 100]);%cascadeforwardnet%t | fitnet | network | feedforwardnet  patternnet
net.trainFcn = 'trainscg'; 
tmp11=mnTr;
tmp12=mnTa;
net = configure(net,tmp11,tmp12);
[trainedNet,net] = train(net,tmp11,tmp12);
nbf=net.best_vperf; %this gives crossentropy measure of performance

net1 = patternnet([200 100 100]);%cascadeforwardnet%t | fitnet | network | feedforwardnet  patternnet
net1.trainFcn = 'trainscg'; 
net1 = configure(net1,tmp11,tmp12);
[trainedNet1,net1] = train(net1,tmp11,tmp12);
nbf1=net1.best_vperf; %this gives crossentropy measure of performance

net2 = patternnet([200 100 100]);%cascadeforwardnet%t | fitnet | network | feedforwardnet  patternnet
net2.trainFcn = 'trainscg'; 
tmp11=mnTr;
tmp12=mnTa;
net2 = configure(net2,tmp11,tmp12);
[trainedNet2,net2] = train(net2,tmp11,tmp12);
nbf2=net2.best_vperf; %this gives crossentropy measure of performance

[Iminall, prombest, tpbestall, fpbestall, CVrawdd]=myfindbestnetwork(circMatdd, trainedNet, trainedNet1, trainedNet2, AMPval);

trainedNetbest={trainedNet,trainedNet1,trainedNet2};
trainedNetbest=trainedNetbest{Iminall};

%% calculate confidence values CV for real data and NBQX Data
test_circ=trainedNetbest(circMatallsweeps);
CVraw=movmean(test_circ(1,:),1);
CVraw(CVraw<0)=0;% confidence values of real data

test_circN=trainedNetbest(circMatallsweeps_NBQX);
CVrawN=movmean(test_circN(1,:),1);
CVrawN(CVrawN<0)=0;% confidence values of NBQX

test_circdd=trainedNetbest(circMatallsweepsdd);
CVrawdd=movmean(test_circdd(1,:),1);
CVrawdd(CVrawdd<0)=0;% confidence values of test data

% Use Peak Finder on CV to Find Mini Locations
thrProm=prombest;
[pks1,locs1,width1,prominence1]=findpeaks(CVraw,'MinPeakDistance',20,'MinPeakProminence',thrProm,'MinPeakWidth',8);
[pks1N,locs1N,width1N,prominence1N]=findpeaks(CVrawN,'MinPeakDistance',20,'MinPeakProminence',thrProm,'MinPeakWidth',8);

%% Find True Positives and False Positives in Synthetic Data
[pksdd,locsdd,widthdd,prominencedd]=findpeaks(CVrawdd,'MinPeakDistance',20,'MinPeakProminence',0,'MinPeakWidth',0);
[tlocs,flocs, trueamppos]=myfindminicvdd(locsdd,size(AMPval,2));
AMPvaltrue=AMPval(truetramppos);
fusedvaltrue=fusedval(trueamppos);
prominenceddreal=prominencedd(tlocs);
prominenceddfake=prominencedd(flocs);
widthddreal=widthdd(tlocs);
widthddfake=widthdd(flocs);
TPdd=zeros(1,999999);
FPdd=zeros(1,999999);
for nn=1:9999
TP1=size(prominenceddreal(prominenceddreal>((nn)/10000)),2)/size(prominenceddreal,2);
FP1=sum(prominenceddfake>(((nn)/10000)))/(size(tlocs,1)-size(prominenceddreal,2));
TPdd(nn)=TP1;
FPdd(nn)=FP1;
end

validation=TPdd(prombest*10000);
falsepositive=FPdd(prombest*10000);
nbqxvreal=sum(prominence1N)/sum(prominence1);
%% Matrix of Identified sweeps
M=circMatallsweeps(:, locs1(prominence1>thrProm));
MN=circMatallsweeps_NBQX(:, locs1N(prominence1N>thrProm));
Mdd= circMatallsweepsdd(:,locsdd(shit(:,2)));
Mddf= circMatdd(:,locsdd(flocs & prominencedd>thrProm));

% calculate amplitudes removes those with excessive error in fitting
[AMPr, risetime, decay]=myampr4(M, width1(prominence1>thrProm));
if 0==isempty(pks1N)
[AMPrN, risetimeN, decayN]=myampr4(MN, width1N(prominence1N>thrProm));
AMPrallN(N,:)=AMPrN;
else
AMPrallN(N,:)=0;
end
[AMPrdd]=myampr4(Mdd, widthdd(shit(:,2)));
[AMPrddf]=myampr4(Mddf, widthdd(~tfrealminidd & prominencedd>thrProm));

AMPrall(N,:)=AMPr;
AMPralldd(N,:)=AMPrdd;
AMPrallddf(N,:)=AMPrddf;

sizeamp=mean(AMPr);
stdevamp=std(AMPr)/sqrt(size(AMPr,1));
rlSN=sizeamp/mean(std(circMatallsweeps_NBQX,1,1));