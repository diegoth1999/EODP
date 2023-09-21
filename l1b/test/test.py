from common.io.writeToa import readToa
import numpy as np
import matplotlib.pyplot as plt

## Test report script
# Compare outputs
# 1. Read LUSS
# 2. Read your outputs
# 3. Compare
directory_luss = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\output'
directory_myoutputs = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\myoutputs'
directory_myoutputsFalseEq = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\myoutputsEqFalse'

tol = 0.01e-2
three_sigma = 0.997

bands = ['VNIR-0','VNIR-1','VNIR-2','VNIR-3']
for band in bands:
    toa_luss = readToa(directory_luss,'l1b_toa_'+band+'.nc')
    toa_myoutputs = readToa(directory_myoutputs,'l1b_toa_'+band+'.nc')
    toa_myoutputsFalseEq = readToa(directory_myoutputsFalseEq,'l1b_toa_'+band+'.nc')

    # Calculate the absolute difference between the matrices
    absolute_difference = np.abs(toa_luss - toa_myoutputs)

    # Count the number of entries where the difference is less than 0.01 percent
    num_entries_below_threshold = np.sum(absolute_difference < 0.0001)

    # Calculate the total number of entries in the matrices
    total_entries = toa_myoutputs.size

    # Calculate the percentage of entries below the threshold
    percentage_below_threshold = (num_entries_below_threshold / total_entries) * 100

    # Check if at least 97.3 percent of the entries are below the threshold
    if percentage_below_threshold >= 97.3:
        print("At least 97.3% of the entries are less than 0.01% different.")
    else:
        print("Less than 97.3% of the entries are less than 0.01% different.")






