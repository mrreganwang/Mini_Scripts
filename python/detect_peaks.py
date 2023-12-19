import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks
import joblib
import pyabf
from parameters import *
from helper_functions.MakeCircMatData import make_circ_mat_data
import matplotlib.pyplot as plt
from mpl_interactions import ioff, panhandler, zoom_factory



# Load the trained neural network
trainedNetN = joblib.load(trained_net_path)

# Load the data trace to best tested
trace = pyabf.ABF(trace_file_path)
trace = trace.sweepY
trace = trace - np.convolve(trace, np.ones(10000) / 10000, mode='same')  # Remove slow data deviations

# Load the noise file
noise_data = np.load(noise_file_path)
ASnoise = noise_data * 2.5              # scale the noise
noiseSD = np.std(ASnoise)               

# Find peaks
data_length = len(trace)
smoothed_trace = trace.copy()
original_trace = smoothed_trace.copy()
TN = trainedNetN
pks = np.array([1])
pk_list = []
nn = ASnoise
iteration = 1
num_peak_table = []

while pks.any():
    # Generate CV
    cv = TN.predict_proba(make_circ_mat_data(SL, smoothed_trace))

    cv = cv[:,1]
    cv[cv < 0] = 0
    cv = np.convolve(cv, np.ones(10) / 10, mode='same')
    # Find CV peaks
    pks, _ = find_peaks(cv, distance=mpd, prominence=0.95, width=mpw)
    pk_list.extend(pks)

    # Subtract detected mini from trace
    if pks.size > 0:
        for n in range(pks.size):
            zero_inval = np.arange(pks[n] - 50 + 70, pks[n] + 70 + 100, 1, dtype=int)
            noise_inval = zero_inval + 100
            noise_inval = np.reshape(noise_inval, (150,1))
            smoothed_trace[zero_inval] = 0
            smoothed_trace[noise_inval] = nn[:len(zero_inval)]

    # Plot CV, original, and remade trace
    # fig = plt.figure()
    # plt.plot(original_trace[70:])
    # plt.plot(smoothed_trace[70:])
    # plt.plot(10 * cv)
    # plt.scatter(pk_list, original_trace[np.array(pk_list) + 71], color='red', marker='o', label='Peaks')
    # plt.legend()
    # plt.show()
    with plt.ioff():
        figure, axis = plt.subplots()
    plt.scatter(pk_list, original_trace[np.array(pk_list) + 71], color='red', marker='o', label='Peaks')
    plt.plot(original_trace[70:])
    
    plt.legend()
    disconnect_zoom = zoom_factory(axis)
    # Enable scrolling and panning with the help of MPL
    # Interactions library function like panhandler.
    pan_handler = panhandler(figure)
    # display(figure.canvas)
    # plt.show()

    print(f"{len(pks)} peaks are found in iteration {iteration}")
    print(pks)
    iteration += 1
    num_peak_table.append([len(pks), len(pk_list)])

# # Save the result
# result_data = {'numPeakTable': np.array(num_peak_table), 'trace': trace, 'Pk_list': np.array(pk_list)}
# save_path = os.path.join(wd, 'result_data.mat')
# savemat(save_path, result_data)