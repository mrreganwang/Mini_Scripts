# for train_net.py
peak_file_path = "helper_functions/FastMini.npy"    # stores the shape of the peak
noise_file_path = "noise_files/2022_02_28_0000_noise.npy"
training_name = "trainedNetN"
result_network_path = "pre_trained_networks/" + training_name + ".pkl"
N = 10000   # training set contains n positives and n negatives
SL = 300    # sweeping length. number of data points for one peak
amp = 2     # hyperparameter, mean peak amplitude of the signals in the training set
mode = 1    # same (1), pearson distributed (2), or uniformly random (3) peak amplitude


# for detect_peaks.py
trace_file_path = "data_to_test/2022_02_28_0000.abf"
trained_net_path = "pre_trained_networks/trainedNetN.pkl"
noise_file_path = "noise_files/2022_02_28_0000_noise.npy"
         
SL = 300    # sweeping length. number of data points for one peak
mpd = 4     # minimum peak distance 
mpw = 3     # minimum peak width