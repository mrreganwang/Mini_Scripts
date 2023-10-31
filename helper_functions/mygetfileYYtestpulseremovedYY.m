%% Import Steph's Data Get File
function [circMatallsweeps, allsweeps,appendedsweeps,filename]=mygetfileYYtestpulseremovedYY(gain1,F,I,d)
%     source_dir=uigetdir(path,selection1);
%     filename=uigetfile([source_dir,'\*.*'],selection1);
%     source_files=fullfile(source_dir, filename);
    
    source_dir=d; %'C:\Users\Chenyu\Desktop\yyminifolder\mini files';
    filename=F{I};
    source_files=fullfile(source_dir, filename);
    fileid=fopen(source_files); 
    arawsweeparray=fread(fileid, 'int16');% figure; plot(arawsweeparray,'.')
%  14563 = calculated sweeplength
    intersweeplength_SA=14563;%15085*2; %30170;
    sweeplength_SA=14563;%15085*2; %30170;
    totallength_SA=length(arawsweeparray)-1024;
    run1num_SA=(((totallength_SA))/intersweeplength_SA); 
    run1num_SA=round(run1num_SA,0);
    gain=gain1; %change this depending on gain gain1=100;

%startsweep=(starttime*5)/intersweeplength_SA;
%endsweep=(endtime*5)/intersweeplength_SA;
    
% Matrix of Sweeps
    startidxsweepsfull_SA=zeros(run1num_SA,1); %Get starting positions of every sweep
for i=1:run1num_SA
    startidxsweepsfull_SA(i)=(1024+((i-1)*intersweeplength_SA));
end

    allsweeps=[]; %loop creates a matrix of every sweep: columns are sweeps, rows are time 
    appendedsweeps=[];
for l=1:numel(startidxsweepsfull_SA)-1
    q=arawsweeparray(((startidxsweepsfull_SA(l)):(startidxsweepsfull_SA(l)+(sweeplength_SA))),1);
    q=q/gain;
    q=movmean(q,2);% Downsample from 10K to 5K
    q=downsample(q,2);
    allsweeps=cat(2,allsweeps,q);
    appendedsweeps=cat(1, appendedsweeps, q);
end

circMatallsweeps=zeros(160, ((sweeplength_SA/2)*(run1num_SA-1)));
circMat=zeros(160,(sweeplength_SA/2));
for nn=0:(size(allsweeps,2))-1
    clear n 
    for n=1:(length(allsweeps)-160)
        circMat(:,n)= allsweeps((n:n+159), (nn+1));
    end
circMatallsweeps(1:160, (nn*size(circMat,2)+1):(nn*size(circMat,2)+size(circMat,2)))=circMat;
end

end