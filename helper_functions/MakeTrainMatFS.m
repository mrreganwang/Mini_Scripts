function [mnTr, mnTa]=MakeTrainMatFS(mininewF,mininewS,AS, D, SL)
clear tmp tmp1 miniTr noiseTr mnTr mnTa
% D=D/2; 
mininew=mininewF';
%make D noisy aligned mini of variable amplitude (skewed amplitude distribution)
tmp1=AS;  %tmp1=tmp1*0;
miniTr=zeros(D,SL);
for n=1:D
    tmp=round(rand*(length(AS)-SL-1)+1);
    miniTr(n,:)=tmp1(tmp+1:tmp+SL)+(mininew(n,:))*10; %add real noise scale mini by noise
end
miniTr=miniTr'; % mini sweeps as columns
%generate D noise sweeps
miniTrF=miniTr;clear miniTr

mininew=mininewS';
miniTr=zeros(D,SL);
for n=1:D
    tmp=round(rand*(length(AS)-SL-1)+1);
    miniTr(n,:)=tmp1(tmp+1:tmp+SL)+(mininew(n,:))*10; %add real noise scale mini by noise
end
miniTr=miniTr';
miniTrS=miniTr;clear miniTr

noiseTr=zeros(D,SL);
for n=1:2*D
    tmp2=round(rand*(length(AS)-SL-1)+1);
    noiseTr(n,:)=tmp1(tmp2+1:tmp2+SL); %add real noise StN*
end
noiseTr=noiseTr'; 
% cat mini and noise sweeps into training matrix
mnTr=cat(2,miniTrS,miniTrF,noiseTr); 

% target matrix
mnTa=[ones(1,2*D),zeros(1,2*D),;...
    zeros(1,2*D),ones(1,2*D)];

end