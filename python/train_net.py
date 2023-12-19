import os
import numpy as np
from scipy.io import loadmat
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
import joblib
from tkinter import Tk, filedialog
from helper_functions.MakeMiniMatFS import make_mini_mat_fs
from helper_functions.MakeTrainMat import make_train_mat
from parameters import *

signal = np.load(peak_file_path)
signal1 = np.load("helper_functions/SlowMini.npy")

noise_trace = np.load(noise_file_path)
noise_trace = noise_trace * 2.5

# Create training data

x = make_mini_mat_fs(N, SL, amp, mode, signal, signal1)[0]
X, y = make_train_mat(-x, noise_trace.flatten(), N, SL)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)

# Initialize neural net
mlp = MLPClassifier(hidden_layer_sizes=(200, 100, 100), activation="logistic", max_iter=10000)

# Train the neural net
trainedNetN = mlp.fit(X, y)

# Save the trained model
joblib.dump(trainedNetN, result_network_path)