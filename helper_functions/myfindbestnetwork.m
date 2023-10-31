function [Iminall, prombest, tpbestall, fpbestall, CVrawdd]=myfindbestnetwork(circMatdd, trainedNet, trainedNet1, trainedNet2, AMPval)
thrProm=.1;

test_circdd=trainedNet(circMatdd);
CVrawdd1=movmean(test_circdd(1,:),1);
CVrawdd1(CVrawdd1<0)=0;% confidence values of circular interleved synthetic minis and noise of Training Data
[~,locsdd,~,prominencedd]=findpeaks(CVrawdd1,'MinPeakDistance',20,'MinPeakProminence',thrProm);
[tlocs,flocs]=myfindminicvdd2(locsdd,size(AMPval,2));
prominenceddreal=prominencedd(tlocs);
prominenceddfake=prominencedd(flocs);
TPdd1=zeros(1,9999);
FPdd1=zeros(1,9999);
for nn=1:9999
TP1=size(prominenceddreal(prominenceddreal>((nn)/10000)),2)/size(prominenceddreal,2);
FP1=sum(prominenceddfake>(((nn)/10000)))/(size(tlocs,1)-size(prominenceddreal,2));
TPdd1(nn)=TP1;
FPdd1(nn)=FP1;
end
[min1,I1]=min(sqrt(((1-TPdd1).^2)+((FPdd1).^2)));
tpbest1=TPdd1(I1);
fpbest1=FPdd1(I1);

test_circdd1=trainedNet1(circMatdd);
CVrawdd2=movmean(test_circdd1(1,:),1);
CVrawdd2(CVrawdd2<0)=0;% confidence values of circular interleved synthetic minis and noise of Training Data
[~,locsdd1,~,prominencedd1]=findpeaks(CVrawdd2,'MinPeakDistance',20,'MinPeakProminence',thrProm);
[tlocs1,flocs1]=myfindminicvdd2(locsdd1,size(AMPval,2));
prominenceddreal1=prominencedd1(tlocs1);
prominenceddfake1=prominencedd1(flocs1);
TPdd2=zeros(1,9999);
FPdd2=zeros(1,9999);
for nn=1:9999
TP1=size(prominenceddreal1(prominenceddreal1>((nn)/10000)),2)/size(prominenceddreal1,2);
FP1=sum(prominenceddfake1>(((nn)/10000)))/(size(tlocs1,1)-size(prominenceddreal1,2));
TPdd2(nn)=TP1;
FPdd2(nn)=FP1;
end
[min2,I2]=min(sqrt(((1-TPdd2).^2)+((FPdd2).^2)));
tpbest2=TPdd2(I2);
fpbest2=FPdd2(I2);

test_circdd2=trainedNet2(circMatdd);
CVrawdd3=movmean(test_circdd2(1,:),1);
CVrawdd3(CVrawdd3<0)=0;% confidence values of circular interleved synthetic minis and noise of Training Data
[~,locsdd2,~,prominencedd2]=findpeaks(CVrawdd3,'MinPeakDistance',20,'MinPeakProminence',thrProm);
[tlocs2,flocs2]=myfindminicvdd2(locsdd2,size(AMPval,2));
prominenceddreal2=prominencedd2(tlocs2);
prominenceddfake2=prominencedd2(flocs2);
TPdd3=zeros(1,9999);
FPdd3=zeros(1,9999);
for nn=1:9999
TP1=size(prominenceddreal2(prominenceddreal2>((nn)/10000)),2)/size(prominenceddreal2,2);
FP1=sum(prominenceddfake2>(((nn)/10000)))/(size(tlocs2,1)-size(prominenceddreal2,2));
TPdd3(nn)=TP1;
FPdd3(nn)=FP1;
end
[min3,I3]=min(sqrt((((1-TPdd3).^2)+((FPdd3).^2))));
tpbest3=TPdd3(I3);
fpbest3=FPdd2(I3);

[~,Iminall]=min([min1,min2,min3]);
% [~,Itpall]=max([tpbest1,tpbest2,tpbest3]);
% [~,Ifpall]=min([fpbest1,fpbest2,fpbest3]);

Ibest=[I1,I2,I3];
prombest=Ibest(Iminall)/10000;

tpbestall=[tpbest1,tpbest2,tpbest3];
tpbestall=tpbestall(Iminall);

fpbestall=[fpbest1,fpbest2,fpbest3];
fpbestall=fpbestall(Iminall);

CVrawdd=[CVrawdd1; CVrawdd2; CVrawdd3];
CVrawdd=CVrawdd(Iminall,:);

end
