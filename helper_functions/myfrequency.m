function [interminidistance3, interminidistance5, interminidistance8]=myfrequency(amplitudethreshold3, amplitudethreshold5, amplitudethreshold8, frequencylocations3, frequencylocations5, frequencylocations8, locs1)
interminidistance3=zeros(1,(size(locs1(amplitudethreshold3),2)-1));    
for n=1:(size(locs1(amplitudethreshold3),2)-1)
    interminidistance3(n)=(frequencylocations3(n+1)-frequencylocations3(n))/5; %dividing by 5 to get mS
end  

interminidistance5=zeros(1,(size(locs1(amplitudethreshold5),2)-1));    
for n=1:(size(locs1(amplitudethreshold5),2)-1)
    interminidistance5(n)=(frequencylocations5(n+1)-frequencylocations5(n))/5; %dividing by 5 to get mS
end  

interminidistance8=zeros(1,(size(locs1(amplitudethreshold8),2)-1));    
for n=1:(size(locs1(amplitudethreshold8),2)-1)
    interminidistance8(n)=(frequencylocations8(n+1)-frequencylocations8(n))/5; %dividing by 5 to get mS
end

end