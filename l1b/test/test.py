from common.io.writeToa import readToa
import numpy as np
import matplotlib.pyplot as plt

## Test report script
# Compare outputs
# 1. Read LUSS
# 2. Read your outputs
# 3. Compare
directory_luss = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\output'
directory_input = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\input'
directory_myoutputs = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\myoutputs'
directory_myoutputsFalseEq = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\myoutputsEqFalse'

plt.figure(figsize=(12, 8))  # Adjust the figure size as needed
bands = ['VNIR-0','VNIR-1','VNIR-2','VNIR-3']
for band in bands:
    toa_luss = readToa(directory_luss,'l1b_toa_'+band+'.nc')
    toa_input = readToa(directory_input,'ism_toa_isrf_'+band+'.nc')
    toa_myoutputs = readToa(directory_myoutputs,'l1b_toa_'+band+'.nc')
    toa_myoutputsFalseEq = readToa(directory_myoutputsFalseEq,'l1b_toa_'+band+'.nc')

    # Calculate the absolute difference between the matrices
    absolute_difference = np.abs(toa_luss - toa_myoutputs)
    print("Absolute difference is:", np.mean(absolute_difference))

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

    # Plot the midpoint values for each band
    plt.subplot(2, 2, bands.index(band) + 1)  # Create a 2x2 grid of subplots
    plt.plot(toa_input[int(toa_input.shape[0] / 2), :], label='TOA Reference')
    plt.plot(toa_myoutputs[int(toa_myoutputs.shape[0] / 2), :], label='Tested TOA')
    plt.plot(toa_myoutputsFalseEq[int(toa_myoutputsFalseEq.shape[0] / 2), :], label='No EQ tested TOA')

    plt.xlabel('ACT [-]')
    plt.ylabel('Radiances [mW/m2/sr]')
    plt.title('Comparison L1B TOA (' + band + ')')
    plt.legend()
    plt.grid(True)

plt.tight_layout()  # Adjust subplot layout for better spacing
plt.show()




