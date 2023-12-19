import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from sklearn.neural_network import MLPClassifier
from tkinter import Tk, filedialog

# Load the neural net
wd = os.path.dirname(os.path.abspath(__file__))
load_path = os.path.join(wd, 'pre_trained_networks', 'trained_net.mat')

data = loadmat(load_path)
trainedNetN = data['trainedNetN']

# Set the parameters
D = 1000
mode1 = 1
mode2 = 2
SL = 300
mpw = 3
mpd = 4
traceLength = 700000
s = 4
numAmp = 1
timeBeAf = 10 * 2
thresholdList = [0.1, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999]
smoothingList = [1, 3, 10, 20, 30, 40, 50, 60, 70, 80]
noiseScale = 2.5

# Select the noise file to be used for testing
root = Tk()
root.withdraw()
fileName = filedialog.askopenfilename(title="Select Noise File", filetypes=[("MAT files", "*.mat")])
root.destroy()

data = loadmat(fileName)
ASnoise = data['ASnoise'] * noiseScale
noiseSD = np.std(ASnoise)

# Pad the noise trace to a longer trace
T = int(round((0.5 * traceLength / ASnoise.size), 0))
ASnoise2 = ASnoise.copy()
for n in range(T):
    ASnoise2 = np.concatenate((ASnoise2, ASnoise[100:]))

# Initialize arrays
TPR = []
FDR = []
FNR = []
TPLOC = []
FPLOC = []
FNLOC = []
SIGNALTRACE = []
NOISETRACE = []
TRUELOC = []
FF = []

for M in range(numAmp):
    Ms = s * (M + 1)
    print('running genROC loop')

    TPRt, FDRt, FNRt, TPLoc, FPLoc, FNLoc, signalTrace, noiseTrace, trueLoc, FFt = genROC5b(
        smoothingList, thresholdList, Ms, mpw, mpd, SL, ASnoise2, trainedNetN, D, mode1, mode2, timeBeAf
    )

    TPR.append(TPRt)
    FDR.append(FDRt)
    FNR.append(FNRt)
    TPLOC.append(TPLoc)
    FPLOC.append(FPLoc)
    FNLOC.append(FNLoc)
    SIGNALTRACE.append(signalTrace)
    NOISETRACE.append(noiseTrace)
    TRUELOC.append(trueLoc)
    FF.append(FFt)

# Plot ROC for different smoothing values
fig, axes = plt.subplots(numAmp, 1, figsize=(10, 10))

for M in range(numAmp):
    for MM in range(len(smoothingList)):
        axes[M].plot(FDR[M][MM], TPR[M][MM], label=f'Smoothing={smoothingList[MM]}', linewidth=2)

    axes[M].set_xlim([0, 1])
    axes[M].set_ylim([0, 1])
    axes[M].set_xlabel('False Discovery Rate')
    axes[M].set_ylabel('True Positive Rate')
    axes[M].set_title(f'Mini Amplitude = {M * s}')
    axes[M].legend()

plt.tight_layout()
plt.show()

# Save results
result_data = {
    'TPR': TPR,
    'FDR': FDR,
    'FNR': FNR,
    'TPLOC': TPLOC,
    'FPLOC': FPLOC,
    'FNLOC': FNLOC,
    'SIGNALTRACE': SIGNALTRACE,
    'NOISETRACE': NOISETRACE,
    'TRUELOC': TRUELOC,
    'FF': FF,
}

save_path = os.path.join(wd, 'result_data.mat')
savemat(save_path, result_data)