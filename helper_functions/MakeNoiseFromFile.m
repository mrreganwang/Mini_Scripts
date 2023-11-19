function [ASnoise]=MakeNoiseFromFile(AS)

% Make noise trace from non-mini data regions
% figure(1);plot(AS)

t=0;n=0;tmp=0;clear x y
while t==0
    n=n+1;SS=[(n-1)*1000+1:n*(1000)];
    figure(1);plot(AS(SS))
    [t1,t2]=getpts(figure(1));
    % figure(2);plot(AS(t1(1):t1(2)))
    if isempty(t1); break;end
    x(n)=t1(1);y(n)=t1(2);
    tmp=t1(2)-t1(1)+tmp;
    SSS=[(n-1)*1000+1+round(x(n)):(n-1)*(1000)+1+round(y(n))];
    figure(2);plot(AS(SSS));title(num2str(tmp))
end
% anneal bits
ASnoise=[];
for n=1:numel(x)-1
    SS=[(n-1)*1000+1+round(x(n)):(n-1)*(1000)+round(y(n))];
    tmp=AS(SS);tmp=tmp-mean(tmp);
    ASnoise=[ASnoise;tmp];
end
figure;plot(ASnoise);hold on;plot(AS)

end
