from matplotlib import pyplot as plt

from common.io.writeToa import readToa
import numpy as np

## Test report script
# Compare outputs
# 1. Read LUSS
# 2. Read your outputs
# 3. Compare

# Statement:
# ❑ Check for all bands that the differences with respect to the output TOA (ism_toa_isrf) are <0.01% for
# at least 3-sigma of the points.
# ❑ Check for all bands that the differences with respect to the output TOA (ism_toa_optical) are <0.01%
# for at least 3-sigma of the points.

directory_luss = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-ISM\output'
directory_input = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-ISM\input\gradient_alt100_act150'
directory_myoutputs = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-ISM\myoutput'

plt.figure(figsize=(12, 8))  # Adjust the figure size as needed
bands = ['VNIR-0','VNIR-1','VNIR-2','VNIR-3']

print('Results of the Test Optical Module')
print('-------------------------------------------------')

for band in bands:
    # OPT
    toa_luss_opt = readToa(directory_luss,'ism_toa_optical_'+band+'.nc')
    toa_myoutputs_opt = readToa(directory_myoutputs,'ism_toa_optical_'+band+'.nc')

    # OPT - Calculate the absolute difference between the matrices for the optical case
    absolute_difference_opt = np.abs(toa_luss_opt - toa_myoutputs_opt)
    print("Absolute difference for the 'ism_toa_optical'  is:", np.mean(absolute_difference_opt))

    # OPT - Count the number of entries where the difference is less than 0.01 percent
    num_entries_below_threshold_opt = np.sum(absolute_difference_opt < 0.0001)
    # OPT - Calculate the total number of entries in the matrices
    total_entries_opt = toa_myoutputs_opt.size


    # OPT - Calculate the percentage of entries below the threshold
    percentage_below_threshold_opt = (num_entries_below_threshold_opt / total_entries_opt) * 100

    # OPT - Check if at least 97.3 percent of the entries are below the threshold
    if percentage_below_threshold_opt >= 97.3:
        print("At least 97.3% of the entries are less than 0.01% different.")
    else:
        print("Less than 97.3% of the entries are less than 0.01% different.")

    # OPT - Plot the midpoint values for each band
    plt.figure(1)
    plt.subplot(2, 2, bands.index(band) + 1)  # Create a 2x2 grid of subplots
    plt.plot(toa_luss_opt[int(toa_luss_opt.shape[0] / 2), :], label='TOA OPT Reference')
    plt.plot(toa_myoutputs_opt[int(toa_myoutputs_opt.shape[0] / 2), :], label='Tested OPT TOA')

    plt.xlabel('ACT [-]')
    plt.ylabel('Radiances [mW/m2/sr]')
    plt.title('Comparison TOA OPT (' + band + ')')
    plt.legend()
    plt.grid(True)

    # ISRF
    toa_luss_isrf = readToa(directory_luss, 'ism_toa_isrf_' + band + '.nc')
    toa_myoutputs_isrf = readToa(directory_myoutputs, 'ism_toa_isrf_' + band + '.nc')

    # ISRF - Calculate the absolute difference between the matrices for the isrf case
    absolute_difference_isrf = np.abs(toa_luss_isrf - toa_myoutputs_isrf)
    print("Absolute difference for the 'ism_toa_isrf'  is:", np.mean(absolute_difference_isrf))

    # ISRF - Count the number of entries where the difference is less than 0.01 percent
    num_entries_below_threshold_isrf = np.sum(absolute_difference_isrf < 0.0001)
    # ISRF - Calculate the total number of entries in the matrices
    total_entries_isrf = toa_myoutputs_isrf.size

    # ISRF - Calculate the percentage of entries below the threshold
    percentage_below_threshold_isrf = (num_entries_below_threshold_isrf / total_entries_isrf) * 100

    # ISRF - Check if at least 97.3 percent of the entries are below the threshold
    if percentage_below_threshold_isrf >= 97.3:
        print("At least 97.3% of the entries are less than 0.01% different.")
    else:
        print("Less than 97.3% of the entries are less than 0.01% different.")

    # ISRF - Plot the midpoint values for each band
    plt.figure(2)
    plt.subplot(2, 2, bands.index(band) + 1)  # Create a 2x2 grid of subplots
    plt.plot(toa_luss_isrf[int(toa_luss_isrf.shape[0] / 2), :], label='TOA ISRF Reference')
    plt.plot(toa_myoutputs_isrf[int(toa_myoutputs_isrf.shape[0] / 2), :], label='Tested ISRF TOA')

    plt.xlabel('ACT [-]')
    plt.ylabel('Radiances [mW/m2/sr]')
    plt.title('Comparison TOA ISRF (' + band + ')')
    plt.legend()
    plt.grid(True)


    # TOA Test Detection Module

    toa_luss_DM = readToa(directory_luss, 'ism_toa_' + band + '.nc')
    toa_myoutputs_DM = readToa(directory_myoutputs, 'ism_toa_' + band + '.nc')

    # TOA Test Detection Module  - Calculate the absolute difference between the matrices for the optical case
    absolute_difference = np.abs(toa_luss_DM - toa_myoutputs_DM)
    print("Absolute difference for the 'ism_toa' for the Test Detection Module  is:", np.mean(absolute_difference))

    # Count the number of entries where the difference is less than 0.01 percent
    num_entries_below_threshold = np.sum(absolute_difference < 0.0001)
    # Calculate the total number of entries in the matrices
    total_entries = toa_myoutputs_DM.size

    # Calculate the percentage of entries below the threshold
    percentage_below_threshold = (num_entries_below_threshold / total_entries) * 100

    # Check if at least 97.3 percent of the entries are below the threshold
    if percentage_below_threshold >= 97.3:
        print("At least 97.3% of the entries are less than 0.01% different.")
    else:
        print("Less than 97.3% of the entries are less than 0.01% different.")

    # Plot the midpoint values for each band
    plt.figure(3)
    plt.subplot(2, 2, bands.index(band) + 1)  # Create a 2x2 grid of subplots
    plt.plot(toa_luss_DM[int(toa_luss_DM.shape[0] / 2), :], label='TOA Detection Module Reference')
    plt.plot(toa_myoutputs_DM[int(toa_myoutputs_DM.shape[0] / 2), :], label='Tested Detection Module TOA')

    plt.xlabel('ACT [-]')
    plt.ylabel('Radiances [mW/m2/sr]')
    plt.title('Comparison TOA Detection Module (' + band + ')')
    plt.legend()
    plt.grid(True)

plt.tight_layout()  # Adjust subplot layout for better spacing
plt.show()



