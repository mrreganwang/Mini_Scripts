import numpy as np
from scipy.signal import resample
from scipy.stats import pearson3, uniform

def make_mini_mat_fs(D, SL, K, mode, signal_shape):
    V = SL
    miniF = signal_shape
    pk, _ = np.unravel_index(np.argmax(miniF,axis=None), miniF.shape)
    strt = 53

    mini = miniF
    AMP = np.zeros(D)
    Fused = np.zeros(D)
    mininew = np.zeros((SL, D))

    for n in range(D):
        FF = 0.6 + 0.4 * np.random.rand()
        minits = miniF
        data_range1 = FF * np.arange(1, (1 / FF) * V + 1)
        data_range1[0] = 1
        minirs = resample(minits, len(data_range1), t=data_range1)
        data_range = np.arange(np.floor((1 / FF) * pk) - pk + 1, V + np.floor((1 / FF) * pk) - pk + 1)
        start = data_range.astype(int)[0]
        end = data_range.astype(int)[-1]
        minitmp = minirs[0][start:end+1]
        
        if mode == 1:
            AMP[n] = K
        elif mode == 2:
            AMP[n] = K * (pearson3.rvs(0.5)**2) / 100
        else:
            AMP[n] = uniform.rvs(loc=K-0.5, scale=1.0)

        Fused[n] = FF
        mininew[:,n] = np.reshape(AMP[n] * minitmp, (300,))

    mininewF = mininew.T
    AMPval = np.mean(AMP)

    return mininewF, AMPval, AMP

'''
used if there are two canonical shapes to detect
'''

def make_mini_mat_fs(D, SL, K, mode, signal_shape, signal_shape1):
    V = SL
    D = int(D/2)
    miniF = signal_shape
    miniS = signal_shape1
    pk, _ = np.unravel_index(np.argmax(miniF,axis=None), miniF.shape)
    strt = 53

    mini = miniF
    AMP = np.zeros(D)
    Fused = np.zeros(D)
    mininew = np.zeros((SL, D))

    for n in range(D):
        FF = 0.6 + 0.4 * np.random.rand()
        minits = miniF
        data_range1 = FF * np.arange(1, (1 / FF) * V + 1)
        data_range1[0] = 1
        minirs = resample(minits, len(data_range1), t=data_range1)
        data_range = np.arange(np.floor((1 / FF) * pk) - pk + 1, V + np.floor((1 / FF) * pk) - pk + 1)
        start = data_range.astype(int)[0]
        end = data_range.astype(int)[-1]
        minitmp = minirs[0][start:end+1]
        
        if mode == 1:
            AMP[n] = K
        elif mode == 2:
            AMP[n] = K * (pearson3.rvs(0.5)**2) / 100
        else:
            AMP[n] = uniform.rvs(loc=K-0.5, scale=1.0)

        Fused[n] = FF
        mininew[:,n] = np.reshape(AMP[n] * minitmp, (300,))

    mininewF = mininew.T
    AMPval = np.mean(AMP)

    pk, _ = np.unravel_index(np.argmax(miniS,axis=None), miniS.shape)
    strt = 53

    mini = miniS
    AMP = np.zeros(D)
    Fused = np.zeros(D)
    mininew = np.zeros((SL, D))

    for n in range(D):
        FF = 0.6 + 0.4 * np.random.rand()
        minits = miniS.T
        data_range1 = FF * np.arange(1, (1 / FF) * V + 1)
        data_range1[0] = 1
        minirs = resample(minits, len(data_range1), t=data_range1)
        data_range = np.arange(np.floor((1 / FF) * pk) - pk + 1, V + np.floor((1 / FF) * pk) - pk + 1)
        start = data_range.astype(int)[0]
        end = data_range.astype(int)[-1]
        minitmp = minirs[0][start:end+1]
        
        if mode == 1:
            AMP[n] = K
        elif mode == 2:
            AMP[n] = K * (pearson3.rvs(0.5)**2) / 100
        else:
            AMP[n] = uniform.rvs(loc=K-0.5, scale=1.0)

        Fused[n] = FF
        mininew[:,n] = np.reshape(AMP[n] * minitmp, (SL,))

    mininewF = np.concatenate((mininewF, mininew.T), axis = 0)
    print(f"mininewF shape: {mininewF.shape}")
    AMPval = np.mean(AMP)

    return mininewF, AMPval, AMP



