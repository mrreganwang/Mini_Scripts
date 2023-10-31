%% Import Steph's Data Get File
function [circMat]=...
    MakeCircMat(D,SL,allsweeps)
%     source_dir=uigetdir(path,selection1);
%     filename=uigetfile([source_dir,'\*.*'],selection1);
%     source_files=fullfile(source_dir, filename);
    
%     source_dir=d; %'C:\Users\Chenyu\Desktop\yyminifolder\mini files';
%     filename=F{I};
%     source_files=fullfile(source_dir, filename);
%     fileid=fopen(source_files); 
%     allsweeps=[];
%     allsweeps=fread(fileid, 'int16');% figure; plot(arawsweeparray,'.')
%     allsweeps=allsweeps/gain1;
%     
%     sweeplength_SA=SL;%15085*2; %30170;
clear swpST circMat tmpMat tmp2Mat
Nswps=(length(allsweeps)/SL);
CswpsPswp=floor(D/Nswps);
% CswpsPswp=round(length(allsweeps)/(Nswps*10))-160;
% NN=CswpsPswp*Nswps;%div by 1000 to ftest
tmpMat=zeros(D,160);
% tmpMat=zeros(Nswps*ceil(CswpsPswp/10),160);
    for n=1:Nswps
%         tic;
        for nn=1: CswpsPswp
            swpST=(n-1)*CswpsPswp;%swpST=ceil(swpST/10);
            swpST=swpST+nn;
            RswpST=nn+(n-1)*(SL);
            tmpMat(swpST,:)=allsweeps(RswpST:RswpST+159)...
                -mean(allsweeps(RswpST:RswpST+159));%add offset if there is a pulse
        end
%         [n,toc]
    end
% tmp2Mat=[tmp2Mat;tmpMat];
% end
circMat=tmpMat;
% 
% circMatallsweeps=zeros(160, length(allsweeps)/sweeplength_SA);
% circMat=zeros(160,(sweeplength_SA/2));
% for nn=0:(size(allsweeps,2))-1
%     clear n 
%     for n=1:(length(allsweeps)-160)
%         circMat(:,n)= allsweeps((n:n+159), (nn+1));
%     end
% circMatallsweeps(1:160, (nn*size(circMat,2)+1):(nn*size(circMat,2)+size(circMat,2)))=circMat;
% end

end