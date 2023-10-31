function [tlocs, flocs, trueamppos]=myfindminicvdd2(locsdd, expectedmininumber)
tlocs=zeros(expectedmininumber,1);
trueamppos=zeros(expectedmininumber,1);
for nnn=0:expectedmininumber-1
    [Mdog, Idog]=min(abs((161+(960*nnn))-locsdd));
     tlocs(nnn+1)=Idog;
     trueamppos(nnn+1)=Mdog;
end
rdog=1:size(locsdd,2);
flocs=setdiff(rdog,tlocs);
trueamppos=trueamppos<50;
