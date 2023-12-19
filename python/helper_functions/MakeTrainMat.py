import numpy as np

def make_train_mat(mininew, AS, D, SL):
    # D = D / 2  # Uncomment this line if D needs to be halved as in the MATLAB script
    miniTr = np.zeros((D, SL))
    noiseTr = np.zeros((D, SL))

    for n in range(D):
        tmp = np.random.randint(0, len(AS) - SL)
        miniTr[n, :] = AS[tmp:tmp + SL] + (mininew[n, :]) * 10

    for n in range(D):
        tmp2 = np.random.randint(0, len(AS) - SL)
        noiseTr[n, :] = AS[tmp2:tmp2 + SL]

    # Concatenate mini and noise sweeps into the training matrix
    mnTr = np.concatenate((miniTr, noiseTr), axis=0)

    # Target matrix

    mnTa = np.concatenate((np.ones((D,1), dtype=int), np.zeros((D,1), dtype=int)), axis = 0)
    mnTa = np.ravel(mnTa)

    return mnTr, mnTa

