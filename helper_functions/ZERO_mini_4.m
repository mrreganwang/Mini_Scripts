
% FOR GENERATING ROC, HEAT MAP, MINI DETECTION IN TRACE AND AMPLITUDE DISTRIBUTION

                % INITIALIZE (1)
AAAList=[4];
TP=[];FP=[];FN=[];TPr=[];FPr=[];FNr=[];FDr=[];CVfirst={};
TN=trainedNetList{1,1};masterPK_list={};D=1000;SL=300;mode1=2;mode2=2;
mpd=1;minpeakwidth=1;
            %MAKE NOISE from loaded ASnoise
T = round((.5*numel(AS)/numel(ASnoise)),0);
ASnoise2 = ASnoise;
for n = 1:T
    ASnoise2 = [ASnoise2;ASnoise(100:end)];
end
%%
tic
for AAA=1:numel(AAAList)

                % MAKE Mpz SECTION
    tmp=zeros(D,SL);AS=AS-movmean(AS,10000);%ASo=AS;AS=ASnoise;
% disp(['mini of amp: ', num2str(AAAList(AAA))]);
    [tmp,Amean,AMP]=MakeMiniMatFS(D,SL,AAAList(AAA),mode1); %1,K;2,Pearsondist
    mininewForBench=-tmp;trueLoc=[];%
    if mode2 == 1
        MNtrace=reshape(mininewForBench',[1 numel(mininewForBench)]);
        MNpZtrace=MakeMNpZ(MNtrace,SL);% alternating minis and zeros
        trueLoc = [71:600:71+600*(D/2-1)];
    else
        %%%%%%%%%%%%%% for random spacing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        spacing = normrnd(500,200,[1,D-1]);
        spacing = round(spacing);spacing=abs(spacing);
        MNpZtrace(1:SL) = tmp(1,:);
        trueLoc(1) = 71;
        for i = 2:1:D
            MNpZtrace = MNpZtrace(1:trueLoc(i-1)+SL-71);
            MNpZtrace = cat(2,MNpZtrace,zeros(1,SL-71+spacing(i-1)));
            MNpZtrace(trueLoc(i-1)+spacing(i-1)-70:trueLoc(i-1)+spacing(i-1)+SL-71) = MNpZtrace(trueLoc(i-1)+spacing(i-1)-70:trueLoc(i-1)+spacing(i-1)+SL-71) + tmp(i,:);
            trueLoc(i) = trueLoc(i-1) + spacing(i-1);
        end

        MNpZtrace = -MNpZtrace;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    TTrueLoc{AAA}=trueLoc;

                % DETECT MINI SECTION
                    % INITIALIZE (2)
    Nn=ASnoise2;MpZ=MNpZtrace;
    endL=min(length(Nn),length(MpZ));
    Nn=reshape(Nn,[1,numel(Nn)]);MpZ=reshape(MpZ,[1,numel(MpZ)]);
    tmp=[];tmp=Nn(1:endL)+MpZ(1:endL); tmp1=tmp;%make test trace = tmp = tmp1
    tmptmp=tmp;tmpinval=[];tmpkeep{AAA}=tmp;
    lastTL=numel(find(trueLoc<numel(tmp1))); %last true location

    NNNList=[.8,.85,.875,.9,.925,.95,.99];% 7 entries
    MMMList=[10,20,30,40,50];% 5 entries
    % NNNList=bestN;MMMList=bestM;
    for MM=1:numel(MMMList);
        for NN=1:numel(NNNList)
% toc;[NN,MM]
    % for MM=1%1:numel(MMMList);
    %     for NN=4%1:numel(NNNList)
            Pk_list=[];nn=1;pks=ones(numel(trueLoc)/20+1,1);CVtmp=[];tmptmp=tmp1;CVtime=1;
            while  numel(pks)>lastTL/50
                % GENERATE CV`
                CVtmpa=TN(MakeCircMatData2(SL,tmptmp));
                CVtmp=CVtmpa(1,:);CVtmp(CVtmp<0)=0;
                CVtmp=movmean(CVtmp,MMMList(MM));
                if CVtime==1; CVfirst{AAA,NN,MM}=CVtmp; CVtime=0;end %save CV for first iteration

                % FIND CV PEAKS
                clear pks
                [~,pks,~,~]=findpeaks...
                    (CVtmp.^1,'MinPeakDistance', mpd,'MinPeakProminence',NNNList(NN),'MinPeakWidth',minpeakwidth);

                Pk_list=[Pk_list,pks]; % add pk to list
                % [nn,numel(Pk_list)]

                % SUBTRACT DETECTED MINI FROM TRACE
               
              if numel(pks)lastTL/50


                    for n=1:numel(pks)
                        inval=[pks(n)-50+70:pks(n)+70+100];
                        tmp(inval)=0;   %zero peak
                        tmp(inval+100)=Nn(1:numel(inval)); %set trace after mini peak to noise values
                        tmpinval=[tmpinval,inval,inval(end)+1:inval(end)+100];
                        % [nn,n,length(tmpinval)]
                    end
                ['numel(pks) = ', num2str(numel(pks)), '  nn = ', num2str(nn)]                  
                    tmptmp(unique(tmpinval))=Nn(1:numel(unique(tmpinval)));
                    nn=nn+1;

                end
                tmpinval=[];

            end

                        % COMPARE DETECTED TO TRUE LOCATIONSc
            AA=[];BB=[];%lastTL=[];
            Pk_list1=sort(Pk_list+71); % offset detected locations to compare to true locations

            AA=ismembertol(Pk_list1,trueLoc,20,'DataScale',1); % detected minis within 40 locations of true locations
            BB=ismembertol(trueLoc(1:lastTL),Pk_list1,20,'DataScale',1);%true locations not detected

            TP=sum(AA==1);TPr(AAA,NN,MM)=TP/lastTL; %true pos rate = TP/total minis
            FP=sum(AA==0);FDr(AAA,NN,MM)=FP/numel(Pk_list1);  %false detection rate= % detected minis not TP
            FN=sum(BB==0);FNr(AAA,NN,MM)=FN/lastTL;  %false neg rate = FN/(total minis)
            B(AAA,NN,MM) = round(((1-TPr(AAA,NN,MM)).^2+FDr(AAA,NN,MM).^2).^.5,2);
            masterPK_list{AAA,NN,MM}=Pk_list;pks=ones(numel(trueLoc)/20+1,1);
            [AAA,NN,MM,TPr(AAA,NN,MM),FDr(AAA,NN,MM)]
 toc;tic          
        end; end
% toc;tic
end

%% plot ROC for diff smoothing values
    % NNNList=[0.80, 0.85 ,0.875, 0.90, 0.925, 0.95, 0.99];% 7 entries
    % MMMList=[10,20,30,40,50];% 5 entries
figure;%AAA=1; load stuff2
for AAA=1:numel(AAAList)
    for MM=1:numel(MMMList)
        subplot (3,2,AAA)
        scatter(FDr(AAA,:,MM),TPr(AAA,:,MM),'LineWidth',2)
        hold on;xlim([0 1]);ylim([0 1]);
    end
    xlabel('False Discovery Rate');ylabel('True Positive Rate')
    title(['mini amplitude =',num2str(AAAList(AAA))])
    % title(['smoothing = ',num2str(MMMList(MM))])
grid on
end
%% heat map
figure
for AAA=1:numel(AAAList)
subplot(2,3,AAA)
% heat map
imagesc(-log10(squeeze(B(AAA,:,:))),[0 1]); title({['mini amplitude =',num2str(AAAList(AAA))];['-log10 dpt']})
xlabel('smoothing');ylabel('threshold')
pause(.1)
end

%% plot trace and mini locations compute mean mini amplitude
OS=0;clear testminisV testminisL %to ID mini peak in trace
AAA=1;NN=1;MM=1; % choose mini amplitude, threshold, smoothing values
tmp1=tmpkeep{AAA};
for n=1:numel(masterPK_list{AAA,NN,MM})
    [testminisV(n),testminisL(n)]=min(tmp1(masterPK_list{AAA,NN,MM}(n)-20+71-OS:masterPK_list{AAA,NN,MM}(n)+20+71-OS));
tmA(n)=mean(tmp1(masterPK_list{AAA,NN,MM}(n)+71-20-5+testminisL(n):masterPK_list{AAA,NN,MM}(n)+71-20+5+testminisL(n)));
end
  f=figure; f.OuterPosition = [10,1200-500,2200,500 ];%POS=POS+200;
            plot(tmp1(71:end));hold on; %plot(tmp(71:end));
            CVT=CVfirst{AAA,NN,MM};CVT=CVT.^10;
            plot(10*CVT(1:end))
            lastTL=numel(find(TTrueLoc{AAA}<numel(tmp1))); %last true location
            scatter(TTrueLoc{AAA}(1:lastTL)-70,tmp1(TTrueLoc{AAA}(1:lastTL)),'filled','black')
            % scatter(masterPK_list{3,4},tmp1(masterPK_list{3,4}+71),'filled')
            scatter((masterPK_list{AAA,NN,MM}+1*testminisL-20-OS),tmp1(masterPK_list{AAA,NN,MM}+70+1*testminisL-20-OS),'filled','magenta')
          % +testminisL-20-OS

figure;histogram(-tmA,[0:50],'Normalization','probability');hold on
FF=3+10*(pearsrnd(10,3,.5,3,1,100000).^2)/100; %ADD 3.5? makes mini with mean 1
[P,edges] = histcounts(-tmA,[0:50],'Normalization','probability');
[N,edges] = histcounts(FF,[0:50],'Normalization','probability');
scatter(4.5,P(5)/(.73),'filled','black')
scatter(4.5,P(5)*(1-.12),'filled','black')
bar(edges(1:end-2)+.5,N(2:end),'FaceColor','red','FaceAlpha',.5,'EdgeAlpha',.5)

%% save stuff
save stuff2 masterPK_list TPr FDr tmpkeep CVfirst B AAAList MMMList NNNList TTrueLoc

%% get weighted heatmap

%
% f(AAA)=
clear wB
% ff=1;
for MM=1:numel(MMMList);
        for NN=1:numel(NNNList)
            ttmp=ff'.*(-log10(B(:,NN,MM)));
wB(NN,MM)=sum(ttmp,1)./sum(ff);
        end
end

figure;imagesc(squeeze(wB))
%%
figure;
for n=1:numel(AAAList)
    subplot(2,3,n)
    imagesc(squeeze(-log10(B(n,:,:))))
end
%% find 'optimal' smoothin and threshold
heatmapVal = ones(9,8);  
heatmapVal = wB;
% change the assignment to be the heatmap values
totalSum = sum(heatmapVal,"all");
weightedColSum = sum(heatmapVal, 1) / totalSum;
weightedRowSum = sum(heatmapVal, 2) / totalSum;
% bestM = round(sum(weightedColSum .* [1:numel(MMMList)], 2));
% bestN = round(sum(weightedRowSum .* [1:numel(NNNList)]', 1));
bestM = round(sum(weightedColSum .* [MMMList], 2));
bestN = round(sum(weightedRowSum .* [NNNList]', 1),2);
%%  here and below is junk
% get weighted mean from heatmap
% tic;load stuff2 ;toc;B
clear BB BBN BBM
    NNNList=[0.80, 0.85 ,0.875, 0.90, 0.925, 0.95, 0.99];% 7 entries
    MMMList=[10,20,30,40,50];% 5 entries
for MM=1:numel(MMMList);
        for NN=1:numel(NNNList)
BB(NN,MM,AAA)=(-log10(B(AAA,NN,MM))).*NNNList(NN).*MMMList(MM);
        end;
end

for NN=1:numel(NNNList)
    BBN(AAA,NN)=mean(BB(NN,:,AAA),2).*NNNList(NN);
end
for MM=1:numel(MMMList)
    BBM(AAA)=sum(BB(:,MM,AAA),2)/sum(NNNList);
end
%% 
PKK=cell2mat(PKS);
% rebuild tmp from tmp2{};
Pk_list=PKS{1};
for n=2:numel(PKS-1)
    for nnn=1:numel(PKS{n})
        PPP=numel(find(PKS{n-1}<PKS{n}(nnn)));
        Pk_list=[Pk_list,PKS{n}(nnn)+(251*numel(PKS{n-1}))];
    end
end

%%
figure;hold on
plot(mean(mean(FDr,1),3),mean(mean(TPr,1),3))
plot(mean(squeeze(mean(FDr,1)),1),mean(squeeze(mean(TPr,1)),1))

[mean(NNNList) mean(MMMList)]
FF=3+10*(pearsrnd(10,3,.5,3,1,100000).^2)/100; %ADD 3.5? makes mini with mean 1
% F=5*(pearsrnd(10,3,.5,3,1,100000).^2)/100;
figure;f=histogram(-FF,[-50:0],'Normalization','probability')
[P,edges] = histcounts(-tmA,[0:50],'Normalization','probability');
for n=1:numel(AAAList)
    ff(n)=N(find(edges>AAAList(n),1));
end