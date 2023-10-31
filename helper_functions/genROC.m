function [TPR,FDR,FNR,TPLoc,FPLoc,FNLoc,signalTrace,smoothedNoiseTrace,trueLoc,FF] = genROC(smoothingList,thresholdList,amplitude,mpw,mpd,SL,noiseTrace,trainedNet,D,mode1,mode2,timeBeAf)
    
    noiseTrace = noiseTrace - movmean(noiseTrace, 10000); % smooth the noise trace                                 % smooth the noise
    disp(['mini of amp: ', num2str(amplitude)]);
    [peaks, ~, ~] = MakeMiniMatFS(D, SL, amplitude, mode1); % Pearsondist mean=1
    mininewForBench = -peaks;

    if mode2 == 1       % equal spacing
        MNtrace = reshape(mininewForBench', [1 numel(mininewForBench)]);
        MNpZtrace = MakeMNpZ(MNtrace,SL);% alternating minis and zeros
        trueLoc = 71:600:71+600*(D/2-1);
    else                % random spacing    
        spacing = normrnd(500,200,[1,D-1]);
        spacing = round(spacing);
        MNpZtrace(1:SL) = peaks(1,:);
        trueLoc = zeros(1,D);
        trueLoc(1) = 71;
        for i = 2:1:D
            MNpZtrace = MNpZtrace(1:trueLoc(i-1)+SL-71);
            MNpZtrace = cat(2,MNpZtrace,zeros(1,SL-71+spacing(i-1)));
            MNpZtrace(trueLoc(i-1)+spacing(i-1)-70:trueLoc(i-1)+spacing(i-1)+SL-71) = MNpZtrace(trueLoc(i-1)+spacing(i-1)-70:trueLoc(i-1)+spacing(i-1)+SL-71) + peaks(i,:);
            trueLoc(i) = trueLoc(i-1) + spacing(i-1);
        end
        MNpZtrace = -MNpZtrace;
    end
    
    ASmod = reshape(noiseTrace, [1 length(noiseTrace)]);
    minL = min(length(MNpZtrace), length(noiseTrace));
    % Bench_AS = ASmod(1:minL) + MNpZtrace(1:minL);   % mini trace + noise
    smoothedNoiseTrace = ASmod(1:minL);                         % noise trace
    signalTrace = MNpZtrace(1:minL);                 % mini trace
    TN = trainedNet;                                % rename the trained network
    
    smoothedTrace = signalTrace + smoothedNoiseTrace;
    originalTrace = smoothedTrace;
    
    numSmooth = numel(smoothingList);
    numThreshold = numel(thresholdList);
    TPR = cell(1,numSmooth);
    FDR = cell(1,numSmooth);
    FNR = cell(1,numSmooth);
    TPLoc = cell(1,numSmooth);
    FPLoc = cell(1,numSmooth);
    FNLoc = cell(1,numSmooth);
    FF = cell(numSmooth, numThreshold);    % store the figures from each iteration

    for MM = 1:1:numSmooth
        % initialize cell arrays
        TPRt = zeros(1,numel(thresholdList));
        FDRt = zeros(1,numel(thresholdList));
        FNRt = zeros(1,numel(thresholdList));
        TPLoct = cell(1,numel(thresholdList));
        FPLoct = cell(1,numel(thresholdList));
        FNLoct = cell(1,numel(thresholdList));

        for NN = 1:1:numThreshold
            [amplitude MM NN]
            pks = 1;                        % store the found peaks in each iteration                                    
            Pk_list = [];                   % store the found peaks across all iterations 
            POS = 0;                        % parameter for plotting figures
            SL = 300;                       % sweeping length. number of data points for one peak
            Nn = noiseTrace;                % value to set the found peaks to after each iteration
            i = 1;                          % keep track of the number of iterations
            F = cell(1,20);                % cell array to store the figures

            while ~isempty(pks)
             
                % GENERATE CV
                CV = TN(MakeCircMatData(SL, smoothedTrace));
                CV = CV(1,:);
                CV(CV < 0) = 0;
                CV = movmean(CV, smoothingList(MM));
             
                % FIND CV PEAKS
                [~,pks,~,~]=findpeaks...
                    (CV,'MinPeakDistance', mpd,'MinPeakProminence',thresholdList(NN),'MinPeakWidth',mpw);
                Pk_list = [Pk_list, pks]; % add peaks found to list
             
                % SUBTRACT DETECTED MINI FROM TRACE
                if ~isempty(pks)
                    for n=1:numel(pks)
                        inval = pks(n)-50+70:pks(n)+70+100;
                        smoothedTrace(inval) = 0;   % zero out peak
                        smoothedTrace(inval+100) = Nn(1:numel(inval)); % set trace after mini peak to noise values
                    end
                end
             
                % PLOT CV ORIGINAL AND REMADE TRACE
                f = figure; 
                f.OuterPosition = [10, 900 - POS, 2200, 500];
                POS = POS + 300;
                plot(originalTrace(71:end));
                hold on; 
                plot(smoothedTrace(71:end));
                plot(10 * CV(1:end));
                scatter(Pk_list, originalTrace(Pk_list+71), 'filled');
                F{i} = f;
                i = i + 1;
            end
            F = F(1:i);
            
            Pk_list1=sort(Pk_list+71); % offset detected locations to compare to true locations
 
            AA = ismembertol(Pk_list1,trueLoc,timeBeAf,'DataScale',1); % are detected minis within 40 locations of true locations?
            lastTL = numel(find(trueLoc<numel(originalTrace))); %last true location in trace
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
            disp(num2str(TPr));
            FDRt(NN) = FDr;
            disp(num2str(FDr));
            FNRt(NN) = FNr;
            disp(num2str(FNr));

            FF{MM,NN} = F;
        end
        TPR{MM} = TPRt;
        FDR{MM} = FDRt;
        FNR{MM} = FNRt;
        TPLoc{MM} = TPLoct;
        FPLoc{MM} = FPLoct;
        FNLoc{MM} = FNLoct;
    end

end