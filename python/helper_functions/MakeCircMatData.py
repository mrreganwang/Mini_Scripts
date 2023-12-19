import numpy as np

def make_circ_mat_data(SL, AS):
    tmp1 = len(AS) % SL
    AS = AS[:-tmp1] if tmp1 > 0 else AS

    tmp_mat2 = np.zeros((SL, len(AS) - SL))
    for n in range(len(AS) - SL):
        tmp_mat2[:, n] = AS[n:n+SL]

    circ_dat_mat = tmp_mat2

    return circ_dat_mat.T