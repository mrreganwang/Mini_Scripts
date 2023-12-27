function [TPR,FDR,FNR,TPLoc,FPLoc,FNLoc,MNpZtrace2,ASmod2,trueLoc,FF] = genROC5(MMList,NNList,A,minpeakwidth,mpd,SL,AS,trainedNet,D,mode1,mode2,timeBeAf, signal)
    
    AS = AS - movmean(AS, 10000);       % smooth the noise
    disp(['mini of amp: ', num2str(A)]);
    % [tmp, ~, ~] = MakeMiniMatFS(D, SL, A, mode1, signal); % Pearsondist mean=1
    [tmp, ~, ~] = MakeMiniMatFS1(D, SL, A, mode1); % Pearsondist mean=1
    mininewForBench = -tmp;

    if mode2 == 1       % equal spacing
        MNtrace = reshape(mininewForBench', [1 numel(mininewForBench)]);
        MNpZtrace = MakeMNpZ(MNtrace,SL);% alternating minis and zeros
        trueLoc = 71:600:71+600*(D/2-1);
    else                % random spacing    
        spacing = normrnd(500,200,[1,D-1]);
        spacing = round(spacing);
        MNpZtrace(1:SL) = tmp(1,:);
        trueLoc(1) = 71;
        for i = 2:1:D
            MNpZtrace = MNpZtrace(1:trueLoc(i-1)+SL-71);
            MNpZtrace = cat(2,MNpZtrace,zeros(1,SL-71+spacing(i-1)));
            MNpZtrace(trueLoc(i-1)+spacing(i-1)-70:trueLoc(i-1)+spacing(i-1)+SL-71) = MNpZtrace(trueLoc(i-1)+spacing(i-1)-70:trueLoc(i-1)+spacing(i-1)+SL-71) + tmp(i,:);
            trueLoc(i) = trueLoc(i-1) + spacing(i-1);
        end
        MNpZtrace = -MNpZtrace;
    end
    
    ASmod = reshape(AS, [1 length(AS)]);
    minL = min(length(MNpZtrace), length(AS));
    % Bench_AS = ASmod(1:minL) + MNpZtrace(1:minL);   % mini trace + noise
    ASmod2 = ASmod(1:minL);                         % noise trace
    MNpZtrace2 = MNpZtrace(1:minL);                 % mini trace
    TN = trainedNet;                                % rename the trained network
    
    tmp = MNpZtrace2 + ASmod2;
    tmp1 = tmp;
    
    TPR = cell(1,numel(MMList));
    FDR = cell(1,numel(MMList));
    FNR = cell(1,numel(MMList));
    TPLoc = cell(1,numel(MMList));
    FPLoc = cell(1,numel(MMList));
    FNLoc = cell(1,numel(MMList));
    FF = cell(numel(MMList), numel(NNList));    % store the figures from each iteration
    for MM = 1:1:numel(MMList)
        % initialize cell arrays
        TPRt = zeros(1,numel(NNList));
        FDRt = zeros(1,numel(NNList));
        FNRt = zeros(1,numel(NNList));
        TPLoct = cell(1,numel(NNList));
        FPLoct = cell(1,numel(NNList));
        FNLoct = cell(1,numel(NNList));
        for NN = 1:1:numel(NNList)
            [A MM NN]
            pks = 1;
            Pk_list = [];
            POS = 0;
            F = cell(1,100);
            iteration = 0;
            tmp = tmp1;
            while ~isempty(pks) && iteration < 5
                iteration = iteration + 1;
                            % GENERATE CV
                f = figure;
                plot(tmp);
                CVtmpa = TN(MakeCircMatData(SL, tmp));
                CVtmp = CVtmpa(1,:);
                CVtmp(CVtmp < 0) = 0;
                CVtmp = movmean(CVtmp, MMList(MM));
             
                               % FIND CV PEAKS
                 [~,pks,~,~]=findpeaks...
                     (CVtmp,'MinPeakDistance', mpd,'MinPeakProminence',NNList(NN),'MinPeakWidth',minpeakwidth);
                 Pk_list = [Pk_list, pks]; % add pk to list
                 numel(Pk_list)
             
                            % SUBTRACT DETECTED MINI FROM TRACE
                     % can for loop be replaced by faster lines??
                     % is this generally applicable to all signals and noise??
                if ~isempty(pks)
                    for n = 1:numel(pks)
                        inval = pks(n)-50+70:pks(n)+70+100;
                        tmp(inval) = 0;   %zero peak
                        tmp(inval+100) = ASmod2(1:numel(inval)); %set trace after mini peak to noise values
                    end
                end
             
                            % PLOT CV ORIGINAL AND REMADE TRACE
                % f = figure; 
                f.OuterPosition = [10,900-POS,2200,500 ];
                POS = POS + 300;
                plot(tmp1(71:end));
                hold on; 
                plot(tmp(71:end));
                plot(10 * CVtmp(1:end))
                scatter(Pk_list, tmp1(Pk_list+71), 'filled');
                lastTL = numel(find(trueLoc < numel(tmp1))); %last true location
                scatter(trueLoc(1:lastTL)-70,tmp1(trueLoc(1:lastTL)-70),'filled','black')
                F{iteration} = f;
                % pause
                close all
                disp([num2str(numel(pks)), ' peaks are found in iteration ', num2str(iteration)]);
            end
            F = F(1:iteration);
            
            Pk_list1=sort(Pk_list+71); % offset detected locations to compare to true locations
 
            AA = ismembertol(Pk_list1,trueLoc,timeBeAf,'DataScale',1); % are detected minis within 40 locations of true locations?
            lastTL = numel(find(trueLoc<numel(tmp1))); %last true location in trace
            BB = ismembertol(trueLoc(1:lastTL),Pk_list1,timeBeAf,'DataScale',1);%true locations not detected
             
            TP = sum(AA == 1);
            TPr = TP / lastTL; %true pos rate = TP/total minis
            FP = sum(AA == 0);
            FDr = FP / numel(Pk_list1);  %false detection rate= % detected minis not TP
            FN = sum(BB == 0);
            FNr = FN / lastTL;  %false neg rate = FN/(total minis)

            TPLoct{NN} = zeros(1,TP);
            FPLoct{NN} = zeros(1,FP);
            FNLoct{NN} = zeros(1,TP);
            TPcount = 1;
            FPcount = 1;
            FNcount = 1;
            for i = 1:1:numel(AA)
                if AA(i) == 1
                    TPLoct{NN}(1,TPcount) = Pk_list1(i);
                    TPcount = TPcount + 1;
                else 
                    FPLoct{NN}(1,FPcount) = Pk_list1(i);
                    FPcount = FPcount + 1;
                end
            end
            for i = 1:1:numel(BB)
                if BB(i) == 0
                    FNLoct{NN}(1,FNcount) = trueLoc(i);
                    FNcount = FNcount + 1;
                end
            end
            
            TPRt(NN) = TPr;
            FDRt(NN) = FDr;
            FNRt(NN) = FNr;

            FF{MM,NN} = F;
        end
        TPR{MM} = TPRt;
        FDR{MM} = FDRt;
        FNR{MM} = FNRt;
        TPLoc{MM} = TPLoct;
        FPLoc{MM} = FPLoct;
        FNLoc{MM} = FNLoct;
    end
close all
end