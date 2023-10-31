function [mnTrM, mnTaM]=MakeTrainMatM(TestAmp,mini1,AS, D, StNset)
% StN -> set signal to noise of mini/noise
%D = # minis generated, D=#noise traces
% mini has peak at pk; length=V; 
AS=reshape(AS,[1 numel(AS)]);
%make D noisy aligned mini of variable amplitude (skewed amplitude distribution)
clear tmp tmp1 miniTr noiseTr mnTrM mnTaM
tmp1=AS;%tmp1=tmp1*0;
miniTr=zeros(D,160);
tic
for n=1:D
    tmp=round(rand*(length(AS)-160-1)+1);
    minitmp=circshift(mini1(n,:),(randi(20)-10),2);
    miniTr(n,:)=tmp1(tmp+1:tmp+160)+minitmp*(TestAmp/2); %add real noise scale mini by noise
end
miniTr1=miniTr;
toc;1
tic
for n=1:D
    tmp=round(rand*(length(AS)-160-1)+1);
    minitmp=circshift(mini1(n,:),(randi(20)-10),2);
    miniTr(n,:)=tmp1(tmp+1:tmp+160)+minitmp*TestAmp; %add real noise scale mini by noise
end
miniTr2=miniTr;
toc;2
tic
for n=1:D
    tmp=round(rand*(length(AS)-160-1)+1);
    minitmp=circshift(mini1(n,:),(randi(20)-10),2);
    miniTr(n,:)=tmp1(tmp+1:tmp+160)+minitmp*1.5*TestAmp; %add real noise scale mini by noise
end
miniTr3=miniTr;

miniTr=[miniTr1;miniTr2;miniTr3];
miniTr=miniTr'; % mini sweeps as columns
toc;3
%generate D noise sweeps
noiseTr=zeros(D,160);
for n=1:D
    tmp2=round(rand*(length(AS)-161)+1);
    noiseTr(n,:)=tmp1(tmp2+1:tmp2+160); %add real noise StN*
end
noiseTr=noiseTr'; 

% cat mini and noise sweeps into training matrix
mnTrM=[miniTr,noiseTr];%cat(2,miniTr,noiseTr); 
mnTrM=mnTrM';
% target matrix
mnTaM=[ones(1,D),zeros(1,D),zeros(1,D),zeros(1,D);...
    zeros(1,D),ones(1,D),zeros(1,D),zeros(1,D);...
    zeros(1,D),zeros(1,D),ones(1,D),zeros(1,D);...
    zeros(1,D),zeros(1,D),zeros(1,D),ones(1,D)];
% % Get tmp5=std(mnTr(:,D+1:end),1,1); 
% StNTR=(mean(AMP)*StNs)/mean(std(mnTr(:,D+1:end),1,1),2); %measured S/N for training set
% StNR=(mean(AMP)*StNs)/n1; %measured S/N for Real
% AMPval=mean(AMP);
% 'here'
end