import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.io import savemat
from sklearn.externals import joblib

def genROC5b(MMList, NNList, A, minpeakwidth, mpd, SL, AS, trainedNet, D, mode1, mode2, timeBeAf):
    AS = AS - np.convolve(AS, np.ones(10000) / 10000, mode='valid')  # Smooth the noise
    print('mini of amp:', A)
    mininewForBench, _, _ = MakeMiniMatFS(D, SL, A, mode1)  # Pearsondist mean=1
    mininewForBench = -mininewForBench

    if mode2 == 1:  # equal spacing
        MNtrace = mininewForBench.reshape((1, -1))
        MNpZtrace = MakeMNpZ(MNtrace, SL)  # alternating minis and zeros
        trueLoc = np.arange(71, 71 + 600 * (D // 2), 600)
    else:  # random spacing
        spacing = np.round(np.random.normal(500, 200, D - 1)).astype(int)
        MNpZtrace = np.zeros(SL)
        MNpZtrace[:SL] = mininewForBench[0, :]
        trueLoc = [71]
        for i in range(1, D):
            MNpZtrace = MNpZtrace[:trueLoc[i - 1] + SL - 71]
            MNpZtrace = np.concatenate([MNpZtrace, np.zeros(SL - 71 + spacing[i - 1])])
            MNpZtrace[trueLoc[i - 1] + spacing[i - 1] - 70:trueLoc[i - 1] + spacing[i - 1] + SL - 71] += mininewForBench[i, :]
            trueLoc.append(trueLoc[i - 1] + spacing[i - 1])
        MNpZtrace = -MNpZtrace

    ASmod = AS.reshape((1, -1))
    minL = min(len(MNpZtrace), len(AS))
    ASmod2 = ASmod[:, :minL]
    MNpZtrace2 = MNpZtrace[:minL]
    TN = trainedNet  # rename the trained network

    tmp = MNpZtrace2 + ASmod2
    tmp1 = tmp.copy()

    TPR = np.zeros((len(MMList), len(NNList)))
    FDR = np.zeros((len(MMList), len(NNList)))
    FNR = np.zeros((len(MMList), len(NNList)))
    TPLoc = []
    FPLoc = []
    FNLoc = []
    MNpZtrace2_all = []
    ASmod2_all = []
    trueLoc_all = []
    FF = np.empty((len(MMList), len(NNList)), dtype=object)

    for MM_idx, MM in enumerate(MMList):
        TPLoc_MM = []
        FPLoc_MM = []
        FNLoc_MM = []
        for NN_idx, NN in enumerate(NNList):
            print(A, MM, NN)
            pks = [1]
            Pk_list = []
            POS = 0
            F = []

            iteration = 0
            tmp = tmp1.copy()
            while pks:
                iteration += 1
                # GENERATE CV
                CVtmpa = TN(MakeCircMatData(SL, tmp))
                CVtmp = CVtmpa[0, :]
                CVtmp[CVtmp < 0] = 0
                CVtmp = np.convolve(CVtmp, np.ones(MM) / MM, mode='valid')

                # FIND CV PEAKS
                pks, _ = find_peaks(CVtmp, distance=mpd, prominence=NN, width=minpeakwidth)
                Pk_list = np.concatenate([Pk_list, pks])  # add pk to list

                # SUBTRACT DETECTED MINI FROM TRACE
                if pks.size > 0:
                    for n in range(pks.size):
                        inval = slice(pks[n] - 50 + 70, pks[n] + 70 + 100)
                        tmp[inval] = 0  # zero peak
                        tmp[inval.start + 100] = ASmod2[0, :len(inval)]  # set trace after mini peak to noise values

                # PLOT CV ORIGINAL AND REMADE TRACE
                f, ax = plt.subplots()
                f.set_size_inches(22, 5)
                ax.plot(tmp1[70:])
                ax.hold(True)
                ax.plot(tmp[70:])
                ax.plot(10 * CVtmp)
                ax.scatter(Pk_list, tmp1[Pk_list + 71], c='r', marker='o', label='Peaks')
                lastTL = len([loc for loc in trueLoc if loc < len(tmp1)])  # last true location
                ax.scatter(trueLoc[:lastTL] - 70, tmp1[trueLoc[:lastTL] - 70], c='k', marker='o', label='True Peaks')
                F.append(f)
                plt.close(f)

            Pk_list1 = np.sort(Pk_list + 71)  # offset detected locations to compare to true locations

            AA = np.isclose(Pk_list1, trueLoc, rtol=timeBeAf, atol=1)  # are detected minis within 40 locations of true locations?
            lastTL = len([loc for loc in trueLoc if loc < len(tmp1)])  # last true location in trace
            BB = np.isclose(trueLoc[:lastTL], Pk_list1, rtol=timeBeAf, atol=1)  # true locations not detected

            TP = np.sum(AA == 1)
            TPr = TP / lastTL  # true pos rate = TP/total minis
            FP = np.sum(AA == 0)
            FDr = FP / len(Pk_list1)  # false detection rate= % detected minis not TP
            FN = np.sum(BB == 0)
            FNr = FN / lastTL  # false neg rate = FN/(total minis)

            TPLoc_MM.append(Pk_list1[AA])
            FPLoc_MM.append(Pk_list1[~AA])
            FNLoc_MM.append(trueLoc[~BB])

            TPR[MM_idx, NN_idx] = TPr
            FDR[MM_idx, NN_idx] = FDr
            FNR[MM_idx, NN_idx] = FNr

            FF[MM_idx, NN_idx] = F

        TPLoc.append(TPLoc_MM)
        FPLoc.append(FPLoc_MM)
        FNLoc.append(FNLoc_MM)

        MNpZtrace2_all.append(MNpZtrace2)
        ASmod2_all.append(ASmod2)
        trueLoc_all.append(trueLoc)

    return TPR, FDR, FNR, TPLoc, FPLoc, FNLoc, MNpZtrace2_all, ASmod2_all, trueLoc_all, FF
