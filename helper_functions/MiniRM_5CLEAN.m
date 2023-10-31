%% (relatively) clean mini program 

% add paths
addpath 'C:\Users\roberto\Documents\MATLAB\function needed for code'
addpath 'C:\Users\roberto\Documents\MATLAB\MiniRM'

%% CHOOSE data folder (either look for it with ui or past below)
datF='C:\Users\roberto\Desktop\DiscussionsW_formerPD\Mini Program Development\Bo\May20 2022\';
datF='C:\Users\roberto\Desktop\DiscussionsW_formerPD\Mini Program Development\Alfonso Data\From Server';
datF='C:\Users\roberto\Desktop\DiscussionsW_formerPD\Mini Program Development\Helmut Data\';
noiseF='C:\Users\roberto\Desktop\DiscussionsW_formerPD\Mini Program Development\Bo\Noise';
fn=uigetfile([datF,'\*.*']);% get filename of data

%% 1b   open data from Stephanie Alfonso (SA) NEED TO FIX

% [SAmat,SAdat]=SAmake_fs(fn,[datF,'\']);%for SA data[SAmat,SAdat]=SAmake_fs(fn,[datF,'\']);%for SA data
% [SAdat]=SAmake_fs(fn,[datF,'\']);%for SA data
% figure;plot(AS(1:10000))
% AS=SAdat/100; % if 100 is gain
%% 1c    open .abf file

AS=abfload([datF,fn]);ASo=AS;AS=AS-movmean(AS,10000);% take out slow data deviations
figure;plot(AS(1:100000))
%% 1ca for Helmut's DATA
AS=ASo;
AS = lowpass(AS,4000,20000);AS=downsample(AS,10);
%% 1cb get  noise file
fnN=uigetfile([datF,'\*.*']);% get filename of noise
load([datF,fnN]);
ASo=AS;ASnoise=ASnoise-movmean(ASnoise,10000);% take out slow data deviations
numel(ASnoise)/numel(AS)
std(ASnoise)
%%  1cc or make noise file
% or make noise file from sections 
% for each window, find region with no minis; 
% left click at start and right click at end of no-mini section
% hit return to go back to main program
[ASnoise]=MakeNoiseFromFile(AS);
std(ASnoise)
tmpn=[datF,fn];tmpn=[tmpn(1:end-4),'Noise2'];
save(tmpn, "ASnoise");
%% LOAD GOOD SAVED NET???
load(['C:\Users\roberto\Documents\MATLAB\MiniRM\','goodnetN.mat'])
D=5000;SL=300;
%% make and save minis from file (fmini or smini) AND SAVE
MakeMiniFromFile(AS)

%% 2a. get fast and slow mini filenames and save
% clear mininewF mininewS  fmfn smfn
% fmfn=input('fastmini filename','s');
% smfn=input('slowmini filename','s');
% % %SlowMiniBo FasMiniBo
% fn=input('filename','s');
% save(fn,'mini')% slowmini2 smini2
%% 2b make D fast and slow minis for NN (and 3. training matrix and 4. NN)
SL=300;K=3;
D=5000;tic;['make minis'] %number of training set minis
[mininewFS]=MakeMiniMatFS(D,SL,K); toc%Pearsondist meanAmp=K
% [mininewS]=MakeMiniMatS(D,SL,smfn); toc%Pearsondist meanAmp=K
% mininewFSo=mininewFS;
%% 3a.Make Training Set use original data or noisefile as noise (1 class mini)
% TRAINING MINI MUST PEAK AT -1 pA
clear mnTr mnTa
%Pearson distribution for training
[mnTr, mnTa]=MakeTrainMat(-mininewFS,reshape(AS,[1 numel(AS) ]),D,SL);% mini mean=StNset
[mnTrN, mnTaN]=MakeTrainMat(-mininewFS,reshape(ASnoise,[1 numel(ASnoise) ]),D,SL);% mini mean=StNset

% plot minis made for training
% figure;plot(mean(mnTrN(:,1:D/2),2));hold on;plot(mean(mnTrN(:,D/2+1:D),2))
% figure;plot(mean(mnTrN(:,1:D/2),2));hold on;plot(mean(mnTrN(:,D/2+1:D),2))
% plot(mean(mnTrN(:,D+1:3*D/2),2));hold on;plot(mean(mnTrN(:,3*D/2:2*D),2))
% figure;for n=1:D/20:D; plot(mnTrN(:,n),'k');hold on;end
%% 4a. make 2 nets one using all file; one using noisefile  (1 class mini)
'make nets'
clear net trainedNet net net0 trainedNetN netN
% tic; 
net0=patternnet([200 100 100]);% other options: cascadeforwardnet%t | fitnet | network | feedforwardnet  patternnet
net=configure(net0,mnTr,mnTa);
netN=configure(net0,mnTrN,mnTaN);
tic;[trainedNet,~] = train(net,mnTr,mnTa);toc
tic;[trainedNetN,~] = train(netN,mnTrN,mnTaN);toc
['time to make net']
% number of iterations in training (in generating some NN; for reference): 62 85 64 35
%% 5. generate ROC

% mini amplitudes tested are N*s;D=5000
clear TP FP LOCS TPt FPt locsB FDRt TP FP FDR NNList%variables in this block
SL=300; minpeakwidth=3;mpd=40;T=round((.5*numel(AS)/numel(ASnoise)),0);% SL=mini sweep length`; CV smoothing; peakfinder argument
ASnoise2=ASnoise;for n=1:T;ASnoise2=[ASnoise2;ASnoise(100:end)];end
% NNList=[.01,.05,.1, .2, .7,.75,.8,.85,.9,.99,.999,.9999,.99999,.999999,.9999999]; % CV thresholds
% NNList=[.001,.01,.02,.03,.1, .2,.9,.999,.99999];
NNList=[.03,.04,.05,.06,.075,.1,.2,.3,.4,.5,.8,.9,.925,.95,.975,.99,.995,.999];
MMList=[1,3,5,10,30,40,50,60,75,90];
s=2;% scale for ROC bench mini
% load FastMiniBo.mat mini % this is basic mini shape loaded from saved file
    clear T1 T2 TP FP LOCS FDR TMS
for M=1:5
    Ms=s*(M);%*std(ASnoise)*2;
clear T1 T2
    'doing genROC2 loop'
    tic
%     [TTP12,FFDR,FFP12,tms,locsB]
    [TPt,FDRt,FPt,tms,tms2,locsBroc]=genROC3(MMList,NNList,Ms,minpeakwidth,mpd,SL,ASnoise2,trainedNetN,D,T); % 1 for CV smoothing
%     TMS(M)=tms;TP(:,M)=TPt; FP(:,M)=FPt;LOCS{M}=locsBroc;FDR(:,M)=FDRt; toc
    TMS2{M}=tms2;TMS{M}=tms;TP{M}=TPt; FP{M}=FPt;
    LOCS{M}=locsBroc;FDR{M}=FDRt;toc
end
%% 5b heat map indicating ((1-TP)^2+FDR^2)^.5 for (smooth,thresh)
ff=figure; ff.Position(1:4)=[50 100 1000 300] 
MMList;NNList;clear t tt B AccuracyHeatmap tmpx tmpy
for n=1:5;for m=1:numel(MMList);
            t(n,m,:)=round(FDR{1,n}{1,m},3);
            tt(n,m,:)=round((TP{1,n}{1,m}),3);
B(n,m,:)=round((((1-tt(n,m,:)).^2+t(n,m,:)).^2).^.5,2);
end;end
C=[0 2];
for n=1:5;ttt=squeeze(B(n,:,:));   subplot(1,5,n);
    imagesc(-log10((ttt)),C);[mv,ml]=min(ttt);[mv1,ml1]=min(mv);
minl=[ml(ml1),ml1];minv=[tt(n,ml(ml1),ml1),t(n,ml(ml1),ml1)];
    title([{num2str(s*(n)),'pA '},{' minL,V=',num2str([minl,minv])}]);
ylabel('CVsmoothing');xlabel('threshold');
AccuracyHeatmap(:,:,n)=-log10(ttt);
tmpA=squeeze(AccuracyHeatmap(:,:,n));
[tmpx(n,:),tmpy(n,:)]=max(tmpA(:));
AccuracyMaxVal(:,n)=tmpx(n,:);
[t1,t2]=ind2sub([size(tmpA)],tmpy(n,:));
A1(n,1)=[t1];A2(n,2)=[t2];
end
AccuracyMaxLoc=[A1(:,1),A2(:,2)];
NNList;
MMList;

%% plot ROC for chosen MMval and threshold
ff=figure; ff.Position(1:4)=[150 100 400 400] 
MMList;NNList; M=5;N=15;
% for n=1:5;for m=numel(MMList)
%             t(n,m,:)=round(FDR{1,n}{1,m},3);
%             tt(n,m,:)=round((TP{1,n}{1,m}),3);
% B(n,m,:)=round((((1-tt(n,m,:)).^2+t(n,m,:)).^2).^.5,2);
% end;end
fdr=squeeze(t);tp=squeeze(tt);
Rfdr=squeeze(fdr(:,M,:));Rtp=squeeze(tp(:,M,:));
cth=5;for n=1:5; scatter(Rfdr(n,:),Rtp(n,:));hold on;end%figure;
sth=1;
for n=1:numel(NNList);
    plot(Rfdr(sth:end,n),Rtp(sth:end,n),'k')
text(Rfdr(sth,n),Rtp(sth,n),num2str(NNList(n)))
text(.1*n+.1,.1,[num2str(round(numel(LOCS{1,1}{1,numel(NNList)-n+1})*2000/numel(AS),2))]);end
% text(0,.1,'FPfreq (Hz)= ');xlim([0 1]);ylim([0 1])
legend(num2str(s),num2str(2*s),num2str(3*s),num2str(4*s),num2str(5*s))
ylabel('True Positive Rate');xlabel('False Discovery Rate');
%% make example
K=8;N=11;NNList(N);D=5000;
tmp=zeros(D,SL);mpd=30;A=4;
% [tmp, AMPval, AMP]=MakeMiniMatKFS(D,A,SL); %Pearsondist mean=1
[tmp]=MakeMiniMatFS(D,SL,K);
mininewForBench=-tmp;%AS=ASnoise;
% Aso=AS;AS=ASnoise;figure;plot(mininewForBench(5000,:))
% ['doing alternating minis and zeros  ',num2str(A)]
MNtrace=reshape(mininewForBench',[1 numel(mininewForBench)]);
MNpZtrace=MakeMNpZ(MNtrace,SL);% alternating minis and zeros
minL=min(length(MNpZtrace),length(ASnoise2));ASmod=reshape(ASnoise2,[1 length(ASnoise2)]);
Bench_AS=ASmod(1:minL)+MNpZtrace(1:minL); % MNpZtrace minis are mean 1

[circDat]=MakeCircMatData2(SL,Bench_AS);
['doing trainedNet ',num2str(A)]
test_circ=trainedNetN(circDat);
CV=movmean(test_circ(1,:),MMList(M)); %can smooth if CV noisy
CV(CV<0)=0;% confidence values of real data
% figure;plot(MNtrace(1:2000));hold on;plot(MNtrace(end-2980:end))
[~,locsM,widthB,prominenceB]=findpeaks...
    (CV,'MinPeakDistance',mpd,'MinPeakProminence',NNList(N),'MinPeakWidth',minpeakwidth);
ff=figure; ff.Position(1:4)=[150 100 400 400] ;
plot(Bench_AS(70:end)-1);hold on;plot(-MNpZtrace(70:numel(CV))+3,'LineWidth',2)
plot(10*CV(1:end)+3,'LineWidth',1);line([1 numel(CV)],3+10*[NNList(N) NNList(N)])
% plot(ASnoise(70:end)+3)
scatter(locsM,-MNpZtrace(locsM+70)+3,'filled')
% scatter(locsM,3+10*NNList(N),'filled','r')
%% figs for paper
%% TP vs FDR for diff smoothing
f=figure;f.Position(1:4)=([500 100 1000 500]);
Mo=3;No=15;
%  ff=figure; ff.Position(1:4)=[150 100 400 400];
hold on;
for n=1:numel(MMList);MMlistString{n}=num2str(MMList(n));end
for M=1:5
    for j=1:18
        for i=1:10
            FDRmat(i,j)=(FDR{1,M}{1,i}(1,j));
            TPmat(i,j)=(TP{1,M}{1,i}(1,j));
        end
    end
    %
    subplot(1,5,M);ax = gca;
    
    for m=1:10
        hold on;
        plot(FDRmat(m,:),TPmat(m,:),'LineWidth',2)
        legend(MMlistString);
    end
scatter(FDRmat(Mo,No),TPmat(Mo,No),'filled')
    xlim([0 1]);ylim([0 1]);grid on;ax.GridAlpha=1;ax.LineWidth = 1.5;
    title([num2str(s*M),' pA']);ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on' ;
end

%% histogram b 
figure;plot(ASnoise,'k')
locsM;numel(CV)/600;
 b = mod(locsM,SL*2);%
 ff=figure; ff.Position(1:4)=[150 100 400 400]; histogram(b,0:3:600,'FaceColor','k')
 ff=figure; ff.Position(1:4)=[150 100 400 400]; histogram(MNpZtrace(locsM+70),100,'FaceColor','k')

%% 5c plot ROC and false discover rate from above values (needs work)
figure(1); hold on;xlim([0 1]); ylim([0 1])
set(gca, 'color', 'none')
clear Leg1 ; MMM=0;
for MM=1:size(FP,2)
    MMM=MMM+1;
    figure(1);plot(FP(:,MM),TP(:,MM),'LineWidth',1.5);hold on;
    Leg1{MMM}=([num2str(MM*.5),' pA']);
end
for NN=2:numel(NNList)
    figure(1);plot(FP(NN,:),TP(NN,:),'k');hold on;
end
plot([0,1],[0,1])
Leg1{MMM+1}='x=y';strNNList=string(NNList);
legend(Leg1,'Location','southeast');
text(FP(:,1),TP(:,1)-.03,strNNList)
% title([fn(1:end-4),'  std=',num2str(std(AS)),' miniN=', num2str(D),]);
xlabel('False Positive Rate');ylabel('True Positive Rate'); grid on

% Plot TP vs False Discovery Rate
figure; hold on;xlim([0 1]); ylim([0 1])
set(gca, 'color', 'none')
clear Leg2; MMM=0;
for MM=1:size(FP,2)
    MMM=MMM+1;
    plot(FDR(:,MM),TP(:,MM),'LineWidth',1.5);hold on
    Leg2{MMM}=([num2str(MM*.5),' pA']);
end
for NN=1:numel(NNList)
    plot(FDR(NN,:),TP(NN,:),'k')
end
plot([0,1],[0,1])
Leg2{MMM+1}='x=y';
legend(Leg2,'Location','southeast')
text(FP(:,1),TP(:,1)-.03,strNNList)
% title([fn(1:end-4),'  std=',num2str(std(AS)),' miniN=', num2str(D),]);
xlabel('False Discovery Rate');ylabel('True Positive Rate'); grid on
NNList; strNNList=string(NNList);
% for n=1:numel(NNList);strNNList(n)=num2str(n);end
% text(FDR(:,3),TP(:,3)-.03,strNNList)
% text(.2,.1,num2str([1:numel(NNList)]))
% text(.2,.2,'choose CVthresh n from here')

%% 7. detect minis (needs work; use thrProm, MMval from ROC)
 N=15;MMval=MMList(M);thrProm=NNList(N); minpeakwidth=3;
%   MMval=MMlist(4);thrProm=.8; minpeakwidth=3;
[circDatMat]=MakeCircMatData2(SL,AS);
clear test_circ
% test_circ=trainedNet(circDatMat);
test_circN=trainedNetN(circDatMat);
clear CV locsB widthB prominenceB TP1 FP1
CV=movmean(test_circN(1,:),MMval); %can smooth if CV noisy

CV(CV<0)=0;% confidence values of real data

[~,locsB1,widthB1,prominenceB1]=findpeaks(CV,'MinPeakDistance',40,'MinPeakProminence',thrProm,'MinPeakWidth',minpeakwidth);

figure;hold on;plot(CV(1:end)*10);% plot(CVN3(1:end)*10)
ASs=movmean(AS,3);
plot(ASs(OS:end));line([0 length(AS)],[thrProm*10 10*thrProm])

clear L;for n=1:numel(locsB1);[V,LL]=min(ASs(locsB1(n)+OS-5:locsB1(n)+OS+10));L(n)=LL(1);locsB1M(n)=locsB1(n)-5+L(n);end
figure; plot(CV(1:end)*10);hold on;plot(ASs(OS:end));line([0 length(AS)],[NNList(N)*10 NNList(N)*10])
% scatter(locsB,10*ones(1,numel(locsB)),'.')
scatter(locsB1M,ASs(locsB1M+OS-1),'filled')
NNList;L;
figure(2);histogram(ASs(locsB1M+OS-1),-35:1:5,'FaceColor','k')
numel(locsB1M)/(length(ASs)/2000)
[2673 529.5]
    %% 8. compute mini amplitudes (needs fixing)
    clear basE1 peakl1 CalcA1 CalcA2 CalcA3 BE1
Bench=0;OS=70;AS=ASo;   locs=locsM;
TTanalyze=Bench_AS;numel(Bench_AS)/300;%AS;
    for n=1:numel(locs)
        [BasE1(n),peakl1(n),CalcA1(n)]=findMiniAmp2(locs(n),TTanalyze,CV,NNList(N),Bench,SL,OS);
   TTanalyzeM=movmean(TTanalyze,5);
        [~,tloc]=min(TTanalyzeM(OS+locs(n)-5:OS+locs(n)+5));
Asimp(n)=mean(TTanalyzeM(OS+locs(n)-5+tloc-2:OS+locs(n)-5+tloc-2));
    end
figure;hold on;histogram(tmpp,50);
text(10,10,[{'mean   std'},{num2str(round([nanmean(CalcA1) nanstd(CalcA1)],1))}])
title('CalcA1 4 pA standard Helmut data');%xlim([ 0 15]);
%histogram(CalcA1,20);
% numel(locsBn)*2000/numel(AS)
tmpp=K*(pearsrnd(10,3,.5,3,1,2500).^2)/100;mean(tmpp)
%%
figure;plot(TTanalyze);hold on;plot(MNpZtrace(1:numel(TTanalyze)))
for n=1:numel(locs);text(locs(n)+OS,-10,num2str(round(Asimp(n),1)));end
%% 9. make traces for example figure for detecting mini of amplitue A
A=3;D=10000;%  SET MINI AMPLITUDE TO DISPLAY
[mininewForBench]=MakeMiniMatK(D,A); %Pearsondist mean=1
clear MNtrace
MNtrace=reshape(mininewForBench',[1 numel(mininewForBench)]);
MNpZtrace=MakeMNpZ(MNtrace,SL);% alternating minis and zeros
minL=min(length(MNpZtrace),length(AS));ASmod=reshape(AS,[1 length(AS)]);
Bench_AS=ASmod(1:minL)+MNpZtrace(1:minL); % MNpZtrace minis are mean 1
['doing MakeCircMatData ',num2str(A)]
[circDatMatBench]=MakeCircMatData2(SL,Bench_AS);
clear test_circBench
['doing trainedNet ',num2str(A)]
test_circBench=trainedNetN(circDatMatBench);
clear CVrawB locsB widthB prominenceB TP1 FP1
CV=movmean(test_circBench(1,:),1); %can smooth if CV noisy
CV(CV<0)=0;% confidence values of real data
[~,locsB,widthB,prominenceB]=findpeaks(CV,'MinPeakDistance',20,'MinPeakProminence',thrProm,'MinPeakWidth',minpeakwidth);
CP=[1:320:numel(CV)];numel(locsB)/numel(CP)

%% 10b detect and compute mini amplitudes (version 2)

clear Bench circDatMat test_circ
Bench=1;SL=160;MMval=1;minpeakwidth=8;%only count benchmark loci
thrProm=NNList(N);
if Bench~=1 % find minis in data
    [circDatMat]=MakeCircMatData2(SL,AS);
    clear test_circ
    test_circ=trainedNet(circDatMat);
    clear CV locsB widthB prominenceB TP1 FP1
    CV=movmean(test_circ(1,:),MMval); %can smooth if CV noisy
    CV(CV<0)=0;% confidence values of real data
    [~,locsB,widthB,prominenceB]=findpeaks(CV,'MinPeakDistance',20,'MinPeakProminence',thrProm,'MinPeakWidth',minpeakwidth);
    clear basE1 peakl1 CalcA1 CalcA2 CalcA3 BE1
    locsB=LL;
    for n=1:numel(locsB)
        [BasE1(n),peakl1(n),CalcA1(n)]=findMiniAmp2(locsB(n),AS,CV,thrProm,Bench);
    end
end
clear basE1 peakl1  CalcA2 CalcA3
if Bench==1 % if minis are all amplitude A
    for n=1:numel(locsB)
        [BasE2(n),peakl2(n),CalcA2(n)]=findMiniAmp2(locsB(n),TpM,CV,thrProm,Bench);
        [BasE3(n),peakl3(n),CalcA3(n)]=findMiniAmp2(locsB(n),MpZ,CV,thrProm,Bench);
    end;
end
%% 11 plot mean and amp histo of minis
L=locsB;BE=BasE2;TtP=TpM;%INDICATE WHICH TRACE TO CONSIDER
HTP=CalcA2;
clear ASm ASm1 ASm2 hFP hFDR hTP
if Bench==1;itmp=(logical((mod(L,320)<5)+(mod(L,320)>315)));% find TP
    L=locsB.*((mod(locsB,320)<=5)+(mod(locsB,320)>=315))...
        +((mod(locsB,320)>5).*(mod(locsB,320)<315))*100;Lm=L(itmp);
else Lm=L;
end
for n=1:numel(Lm)
    % ASm1(n,:)=ASo(L(n):L(n)+160)-mean(ASo(L(n)-BE(n)-20:L(n)-BE(n)));
    ASm2(n,:)=TtP(Lm(n):Lm(n)+159);%-mean(TtP(Lm(n)-BE(n)-20:Lm(n)-BE(n)));
end
figure; plot(mean(ASm2,1),'LineWidth',2);hold on;%plot(mean(ASm2,1));
title(['mean of all ',num2str(numel(HTP)),'  minis'])
if Bench~=1
    figure;h=histogram(HTP,[0:std(AS)/2:5*std(AS)]);hold on;%histogram(CalcA1,100); % histogram(CalcA2,100);
    vv=h.Values; be=h.BinEdges;
    for n=1:numel(vv); hFDR(n)=(FDR(12,n)).*vv(n);hTP(n)=(1/TP(12,n)).*vv(n);end
    hold on;scatter(be(3:end)-(std(AS)/4),hFDR,10,'red','filled')
    hold on;scatter(be(3:end)-(std(AS)/4),hTP,10,'blue','filled')
    title([fn(1:end-4),'  CVt=',num2str(thrProm)]);
    xlabel('amplitude (pA)');ylabel('count')
    text(0,max(vv),{['miniN= ',num2str(numel(HTP))];...
        ,['  freq= ',num2str(round(10000*numel(HTP)/(numel(AS)),1)),' Hz']})
end
if Bench==1; figure
    h=histogram(HTP,[0:A/10:3*A]);title([fn(1:end-4),'  testA=',num2str(A)]); mv=max(h.Values);
    xlabel('amplitude (pA)');ylabel('count')
    tmp=HTP(itmp);tmp(tmp==0)=[];text(0,numel(tmp)^.6,['std= ',num2str(round(std(tmp),2))])
    text(0,numel(tmp)^.75,['mean= ',num2str(round(mean(tmp),2))])
end
hold on; x = [0:.1:3]; y = .8*mv*normpdf(x,A,std(tmp)); plot(x,y,'LineWidth',2)
figure; histogram(CV,100);title(['CV std= ',num2str(std(CV))])
% FN='2022-02-28-0001 02Filter.ABF';
% std(CalcA2)
%%
[~,locsM,widthB,prominenceB]=findpeaks(CV,'MinPeakDistance',20,'MinPeakProminence',.000001,'MinPeakWidth',minpeakwidth);

%% 12 plot minis individually to check
figure
BE=BasE1;
for n=1:numel(L)
    plot(AS(L(n)-100:L(n)+160)-mean(AS(L(n)-BE(n)-10+60:L(n)-BE(n)+60)));hold on%
    plot(1*CV(L(n)-100-60:L(n)+160-60))
    bas=mean(AS(60+L(n)-BE(n)-10:60+L(n)-BE(n))); %use L=CV rather than peak of noisy trace
    PK=mean(AS(L(n)-3+60:L(n)+3+60));
    text(10,0,[num2str(CalcA1(n)),'  100*CV= ',num2str(100*CV(L(n))),...
        '  miniN=',num2str(n),'  BV=',num2str(-mean(AS(60+L(n)-BE(n)-10:60+L(n)-BE(n)))),...
        '  bas=',num2str(bas),'  PK=',num2str(PK)])
    line([160-BasE1(n) 160-BasE1(n)],...
        [0 -CalcA1(n)],'LineWidth',2.5)
    line([100+60 100+60],[PK-bas 1])
    % if 1(n)>0; ylim([-CalcA1(n)-1 1 ]);end
    pause; clf
end
%% 13 plot minis by group
[NminiA,AvMinA]=PlotMBG(AS,L,CalcA1,peakl1);
[NminiA]
[AvMinA]
%% 14 plot mini of amplitude A at location L
figure;
plot((TpM(60:M)+5),'Color','k','LineWidth',1);hold on;
plot(MpZ(60:M),'LineWidth',1);
% plot(2*CV);
% scatter(LL,TpM(LL+60)+5,'black','Marker','o','MarkerFaceColor','black','SizeVariable',1); %ID minis
% scatter(LL-BasE1,CV(LL-BasE1)*2,'black','Marker','.','SizeVariable',1); %ID minis
% scatter(LL-BasE1,MpZ(LL+60-BasE1),'green','Marker','.','SizeVariable',1); %ID minis
line([LL ;LL],[MpZ(LL+60) ;TpM(LL+60)+5],'LineStyle','--','Color','red','LineWidth',1.5)
% for n=1:numel(L)
%
% end
%% make mini
[Mmini]=MakeMiniFromFile(AS);

mini=smini;
fn=input('filename','s');
% vn=input('smini or fmini?');

save(fn,'mini')% slowmini2 smini2


figure;hold on;plot(mini)

NNList

%, ZMinorGrid — Minor grid lines