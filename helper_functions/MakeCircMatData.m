%% Import Steph's Data Get File
function [circDatMat]=...
    MakeCircMatData2(SL,AS)

clear tmp1 circMat tmpMat tmp2Mat

tmp1=mod(length(AS),SL);AS(end-tmp1+1:end)=[];

tmpMat2=zeros(SL,numel(AS));
for n=1:numel(AS)-SL
tmpMat2(:,n)=AS(1+(n-1):SL+(n-1));
end
circDatMat=tmpMat2;
% 
% for n=1:floor(Nswps)-1
% %     tic;
%     for nn=1: CswpsPswp
% %         swpST=(n-1)*(CswpsPswp+160);%swpST=ceil(swpST/10);
% %         swpST=swpST+nn;
%         RswpST=nn+(n-1)*(SL);
%         tmpMat(RswpST,:)=AS(RswpST:RswpST+159)...
%             -mean(AS(RswpST:RswpST+159));%add offset if there is a pulse
%     end
% %     [n,toc]
% end
% % tmp2Mat=[tmp2Mat;tmpMat];
% % end
% circDatMat=tmpMat;
% % 
% % circMatallsweeps=zeros(160, length(allsweeps)/sweeplength_SA);
% % circMat=zeros(160,(sweeplength_SA/2));
% % for nn=0:(size(allsweeps,2))-1
% %     clear n 
% %     for n=1:(length(allsweeps)-160)
% %         circMat(:,n)= allsweeps((n:n+159), (nn+1));
% %     end
% % circMatallsweeps(1:160, (nn*size(circMat,2)+1):(nn*size(circMat,2)+size(circMat,2)))=circMat;
% % end

end