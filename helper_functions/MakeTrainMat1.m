function [mnTr, mnTa]=MakeTrainMat1(mininew,AS, D, SL, StNset)
% mininew=E;%                        (E,AS1,D,SL,StNset)                                 
% (E,AS,D,SL,StNset)
% % D=D/2; 
% mininew=mininew';
% %make D noisy aligned mini of variable amplitude (skewed amplitude distribution)
% clear tmp tmp1 miniTr noiseTr mnTr mnTa
% tmp1=AS;  %tmp1=tmp1*0;
miniTr=zeros(D,SL);
tmp1=AS1;tmp=tmp1;
for n=1:D
    tmp=round(rand*(length(AS)-SL-1)+1);
    miniTr(n,:)=tmp1(tmp+1:tmp+SL)+(mininew(n,:))*StNset; %add real noise scale mini by noise
end
miniTr=miniTr'; % mini sweeps as columns

%generate D noise sweeps
noiseTr=zeros(D,SL);
for n=1:D
    tmp2=round(rand*(length(AS)-SL-1)+1);
    noiseTr(n,:)=tmp1(tmp2+1:tmp2+SL); %add real noise StN*
end
noiseTr=noiseTr'; 

% cat mini and noise sweeps into training matrix
mnTr=cat(2,miniTr,noiseTr); 

% target matrix
mnTa=[ones(1,D),zeros(1,D);zeros(1,D),ones(1,D)];

end