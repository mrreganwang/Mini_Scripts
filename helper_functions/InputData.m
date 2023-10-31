%% Import Steph's Data Get File
function [allsweeps]=...
    InputData(gain1,F,I,d,SL)
%     source_dir=uigetdir(path,selection1);
%     filename=uigetfile([source_dir,'\*.*'],selection1);
%     source_files=fullfile(source_dir, filename);
    
    source_dir=d; %'C:\Users\Chenyu\Desktop\yyminifolder\mini files';
    filename=F{I};
    source_files=fullfile(source_dir, filename);
    fileid=fopen(source_files); 
% Tt = readtable(source_files);
% Tdat=table2array(Tt);
    allsweeps=[];
    allsweeps=fread(fileid, 'int16');% figure; plot(arawsweeparray,'.')
    allsweeps=allsweeps/gain1;
    
    sweeplength_SA=SL;%15085*2; %30170;
%     if addmin==1
%         clear tmp
%         tmp=traceWmini(1:length(allsweeps));
%        allsweeps=allsweeps+tmp'; 
%         
%     end
%     allsweeps=allsweepsmm3;
% clear swpST circMat tmpMat tmp2Mat
% Nswps=(length(allsweeps)/sweeplength_SA);
% CswpsPswp=floor(D/Nswps);
% tmpMat=zeros(length(allsweeps),160);
% for n=1:Nswps
%     for nn=1:SL-160
%         tmpN=nn+(n-1)*(SL-160);
%         tmpNdat=nn+(n-1)*(SL-160);
%     tmpMat(tmpN,:)=(allsweeps(tmpNdat:tmpNdat+159));
%     end
% end
% circDatMat=tmpMat;
% 
% % CswpsPswp=round(length(allsweeps)/(Nswps*10))-160;
% % NN=CswpsPswp*Nswps;%div by 1000 to ftest
% tmpMat=zeros(D,160);
% % tmpMat=zeros(Nswps*ceil(CswpsPswp/10),160);
%     for n=1:Nswps
% %         tic;
%         for nn=1: CswpsPswp
%             swpST=(n-1)*CswpsPswp;%swpST=ceil(swpST/10);
%             swpST=swpST+nn;
%             RswpST=nn+(n-1)*(SL);
%             tmpMat(swpST,:)=allsweeps(RswpST:RswpST+159)...
%                 -mean(allsweeps(RswpST:RswpST+159));%add offset if there is a pulse
%         end
% %         [n,toc]
%     end
% % tmp2Mat=[tmp2Mat;tmpMat];
% % end

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