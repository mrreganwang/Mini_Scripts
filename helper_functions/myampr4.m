%Check sum of squares for Amplitude
function [AMPr, risetime, decay, terror, berror, derror]=myampr4(Mtemp1, width)
AMPr=zeros(size(Mtemp1,2),1);
risetime=zeros(size(Mtemp1,2),1);
decay=zeros(size(Mtemp1,2),1);
terror=zeros(size(Mtemp1,2),1);
berror=zeros(size(Mtemp1,2),1);
derror=zeros(size(Mtemp1,2),1);
% Mtemp1=smooth(Mtemp1);

for N=(1:size(Mtemp1,2))%finds max min and baseline of mini
w=round(width(N));
   if w<57
    [minM, minI]=min(Mtemp1(55:65, N)); 
    avgmini=(minM+Mtemp1((minI+54)-1, N)+Mtemp1((minI+54)+1, N))/3;
    meanM=mean(Mtemp1((60-(w+3)):(60-(w-3)),N));%baseline based on width
    topdiff=abs(avgmini-meanM);
    
    %find current at which 90% of topdiff occurs and 10% of topdiff
    tenpct=((topdiff*.9)+Mtemp1((54+minI),N));
    ninetypct=((topdiff*.1)+Mtemp1((54+minI),N));
    [mdifft,tenpctpt]=min(abs(Mtemp1(50:60,N)-tenpct)); %could use meanM to subtract may want to change 50 to toppoint to constrain
    [mdiffb,ninetypctpt]=min(abs(Mtemp1(50:(54+minI),N)-ninetypct));
    
        %rise time if statements to get at interperlated pts 
        if Mtemp1(tenpctpt+49+1,N)<tenpct && Mtemp1(49+tenpctpt,N)>tenpct
           truetenpctpt=(49+tenpctpt)+((Mtemp1(49+tenpctpt,N)-tenpct)/((Mtemp1(49+tenpctpt,N)-Mtemp1(49+tenpctpt+1,N))));
           terrortemp=0;
        elseif Mtemp1(tenpctpt+49-1,N)>tenpct && Mtemp1(49+tenpctpt,N)<tenpct
           truetenpctpt=(49+tenpctpt-1)+((Mtemp1(49+tenpctpt-1,N)-tenpct)/((Mtemp1(49+tenpctpt-1,N)-Mtemp1(49+tenpctpt,N))));
           terrortemp=0; 
        else
           truetenpctpt=tenpctpt+49; 
           terrortemp=mdifft;
        end   

        if Mtemp1(ninetypctpt+49+1,N)<ninetypct && Mtemp1(49+ninetypctpt,N)>ninetypct
           trueninetypctpt=(49+ninetypctpt)+((Mtemp1(49+ninetypctpt,N)-ninetypct)/((Mtemp1(49+ninetypctpt,N)-Mtemp1(49+ninetypctpt+1,N))));
           berrortemp=0; 
        elseif Mtemp1(ninetypctpt+49-1,N)>ninetypct && Mtemp1(49+ninetypctpt,N)<ninetypct
           trueninetypctpt=(49+ninetypctpt-1)+((Mtemp1(49+ninetypctpt-1,N)-ninetypct)/((Mtemp1(49+ninetypctpt-1,N)-Mtemp1(49+ninetypctpt,N))));
           berrortemp=0; 
        else
           trueninetypctpt=ninetypctpt+49; 
           berrortemp=mdiffb;
        end   
          
    risetimetemp=(trueninetypctpt-truetenpctpt)/5; %ms
    
        %Find decay
        [merrord,toptenpctd]=min(abs(Mtemp1(54+minI:100,N)-tenpct));
        if Mtemp1(toptenpctd+54+minI-1,N)<tenpct && Mtemp1(54+minI+toptenpctd,N)>tenpct
           trueninetypctptdecay=(49+toptenpctd)+((Mtemp1(49+toptenpctd,N)-tenpct)/((Mtemp1(49+toptenpctd,N)-Mtemp1(49+toptenpctd+1,N))));
           derrortemp=0; 
        elseif Mtemp1(toptenpctd+49+1,N)>tenpct && Mtemp1(54+minI+toptenpctd,N)<tenpct
           trueninetypctptdecay=(49+toptenpctd-1)+((Mtemp1(49+toptenpctd-1,N)-tenpct)/((Mtemp1(49+toptenpctd-1,N)-Mtemp1(49+toptenpctd,N))));
           derrortemp=0; 
        else
           trueninetypctptdecay=toptenpctd+54+minI; 
           derrortemp=merrord;
        end   
    
    decaytemp=(trueninetypctptdecay-(54+minI))/5; %ms
   
   else
    [minM, minI]=min(Mtemp1(55:65, N)); 
    avgmini=(minM+Mtemp1((minI+54)-1, N)+Mtemp1((minI+54)+1, N))/3;
    meanM=mean(Mtemp1((1):(6),N));%baseline based on width
    topdiff=abs(avgmini-meanM);

    tenpct=((topdiff*.9)+Mtemp1((54+minI),N));
    ninetypct=((topdiff*.1)+Mtemp1((54+minI),N));
    [mdifft,tenpctpt]=min(abs(Mtemp1(50:60,N)-tenpct)); %could use meanM to subtract may want to change 50 to toppoint to constrain
    [mdiffb,ninetypctpt]=min(abs(Mtemp1(50:(54+minI),N)-ninetypct));
    
    %rise time
    risetimetemp=(ninetypctpt-tenpctpt)/5; %ms
    terrortemp=mdifft;
    berrortemp=mdiffb;
    %Find decay
    [merrord,toptenpctd]=min(abs(Mtemp1(54+minI:100,N)-tenpct));
    decaytemp=(toptenpctd-(54+minI))/5; %ms
    derrortemp=merrord;
   end
   
AMPr(N)=topdiff;
risetime(N)=risetimetemp;
decay(N)=decaytemp;
terror(N)=terrortemp;
berror(N)=berrortemp;
derror(N)=derrortemp;
 end

end
