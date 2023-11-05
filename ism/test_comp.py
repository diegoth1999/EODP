from matplotlib import pyplot as plt

from common.io.writeToa import readToa

directory_isrf = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-ISM\myoutput\newOutput'
directory_l1b = r'C:\Users\diego\PycharmProjects\EODP_TER_2021\EODP-TS-L1B\myoutputs\newOutput'

plt.figure(figsize=(12, 8))  # Adjust the figure size as needed
bands = ['VNIR-0', 'VNIR-1', 'VNIR-2', 'VNIR-3']

for band in bands:
    toa_isrf = readToa(directory_isrf, 'ism_toa_isrf_' + band + '.nc')
    toa_l1b = readToa(directory_l1b, 'l1b_toa_' + band + '.nc')


    # Plot the midpoint values for each band
    plt.subplot(2, 2, bands.index(band) + 1)  # Create a 2x2 grid of subplots
    plt.plot(toa_isrf[int(toa_isrf.shape[0] / 2), :], label='ISRF TOA')
    plt.plot(toa_l1b[int(toa_l1b.shape[0] / 2), :], label='L1B TOA')

    plt.xlabel('ACT [-]')
    plt.ylabel('Radiances [mW/m2/sr]')
    plt.title('Comparison between TOA radiances for (' + band + ')')
    plt.legend()
    plt.grid(True)

plt.tight_layout()  # Adjust subplot layout for better spacing
plt.show()
