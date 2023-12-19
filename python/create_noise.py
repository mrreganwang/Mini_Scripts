import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from tkinter import Tk, filedialog
import pyabf

# Set the correct working directory
wd = os.path.dirname(os.path.abspath(__file__))
os.chdir(wd)

# Prompt window to select the file
root = Tk()
root.withdraw()
fileName = filedialog.askopenfilename(title="Select Data File", filetypes=[("ABF files", "*.abf")])
root.destroy()

# Load the data
AS = np.load(fileName)
AS = AS - np.convolve(AS, np.ones(10000) / 10000, mode='valid')  # Remove slow data deviations

# Create the noise trace
fig, ax = plt.subplots()
ax.plot(AS)
ax.set(title="Select the region with only noise (left-click start, right-click end)")
noise_regions = plt.ginput(n=-1, timeout=-1)
plt.close()

ASnoise = np.zeros_like(AS)
for start, end in noise_regions:
    start, end = int(start), int(end)
    ASnoise[start:end] = AS[start:end]

ASnoise = ASnoise - np.convolve(ASnoise, np.ones(10000) / 10000, mode='valid')  # Remove slow data deviations

# Save the noise trace
noise_file_name = os.path.splitext(os.path.basename(fileName))[0] + '_noise.mat'
save_path = os.path.join(wd, 'data_to_test', noise_file_name)
savemat(save_path, {"ASnoise": ASnoise})