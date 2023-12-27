function [out]=...
    MakeMNpZ(MNtrace,SL)
odd=[1:2:numel(MNtrace)];
even=[0:2:numel(MNtrace)];
% SL=3*SL;
Z=zeros(1,SL);
for n=1:numel(MNtrace)/(SL)
  tmptrace(even(n)*SL+1:even(n)*SL+SL)=MNtrace((n-1)*SL+1:n*SL);
  tmptrace(odd(n)*SL+1:odd(n)*SL+SL)=Z; 
end
tmp=[tmptrace(1:numel(MNtrace)/2)];
    tmp=[tmp,tmptrace(numel(MNtrace)+1:(3/2)*numel(MNtrace))];
% tmptrace(numel(MNtrace)+1:end)=[];
out=tmp;
end  
% figure;hold on; plot(tmp(1:2000))