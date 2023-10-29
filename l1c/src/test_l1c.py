import self
from matplotlib import pyplot as plt, cm

from common.io.writeToa import readToa
import numpy as np

## Test report script
# Compare outputs
# 1. Read LUSS
# 2. Read your outputs
# 3. Compare

# Statement:
# ❑ Check for all bands that the differences with respect to the output TOA are <0.01% for 3-sigma of
# the points.

directory_luss = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1C\output'
directory_input = r"C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1C\input\gm_alt100_act_150,C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1C\input\l1b_output"
directory_myoutputs = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1C\myoutput'


plt.figure(figsize=(12, 8))  # Adjust the figure size as needed
bands = ['VNIR-0','VNIR-1','VNIR-2','VNIR-3']


for band in bands:

    toa_luss_prev = readToa(directory_luss,'l1c_toa_'+band+'.nc')
    toa_myoutputs_prev = readToa(directory_myoutputs,'l1c_toa_'+band+'.nc')
    toa_luss_sort = np.sort(toa_luss_prev)
    toa_myoutputs_sort = np.sort(toa_myoutputs_prev)

    # Calculate the absolute difference between the matrices
    absolute_difference = np.abs(toa_luss_sort - toa_myoutputs_sort)
    print("Absolute difference is:", np.mean(absolute_difference))

    # Count the number of entries where the difference is less than 0.01 percent
    num_entries_below_threshold = np.sum(absolute_difference < 0.0001)
    # Calculate the total number of entries in the matrices
    total_entries = toa_myoutputs_sort.size


    # Calculate the percentage of entries below the threshold
    percentage_below_threshold = (num_entries_below_threshold / total_entries) * 100

    # OPT - Check if at least 97.3 percent of the entries are below the threshold
    if percentage_below_threshold >= 97.3:
        print("At least 97.3% of the entries are less than 0.01% different.")
    else:
        print("Less than 97.3% of the entries are less than 0.01% different.")

#     # Plot the values of toa_luss_prev and toa_myoutputs_prev
#     plt.subplot(2, 2, bands.index(band) + 1)  # Create a 2x2 grid of subplots
#     x_values = np.arange(toa_luss_prev.size)  # Create x-axis values
#     plt.plot(x_values, toa_luss_prev, label='TOA Reference')
#     plt.plot(x_values, toa_myoutputs_prev, label='Tested TOA')
#
#     plt.xlabel('ACT [-]')
#     plt.ylabel('Radiances [mW/m2/sr]')
#     plt.title('Comparison TOA L1C (' + band + ')')
#     plt.legend()
#     plt.grid(True)
#
# plt.tight_layout()  # Adjust subplot layout for better spacing
# plt.show()
#
