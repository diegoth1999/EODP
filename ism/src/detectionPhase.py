import os

import scipy.constants
from ism.src.initIsm import initIsm
import numpy as np
from common.io.writeToa import writeToa
from common.plot.plotMat2D import plotMat2D
from common.plot.plotF import plotF


class detectionPhase(initIsm):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

        # Initialise the random see for the PRNU and DSNU
        np.random.seed(self.ismConfig.seed)

    def compute(self, toa, band):

        self.logger.info("EODP-ALG-ISM-2000: Detection stage")

        # Irradiance to photons conversion
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-2010: Irradiances to Photons")
        area_pix = self.ismConfig.pix_size * self.ismConfig.pix_size  # [m2]
        toa = self.irrad2Phot(toa, area_pix, self.ismConfig.t_int, self.ismConfig.wv[int(band[-1])], band)

        self.logger.debug("TOA [0,0] " + str(toa[0, 0]) + " [ph]")

        # Photon to electrons conversion
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-2030: Photons to Electrons")
        toa = self.phot2Electr(toa, self.ismConfig.QE, self.ismConfig.FWC, band)

        self.logger.debug("TOA [0,0] " + str(toa[0, 0]) + " [e-]")

        if self.ismConfig.save_after_ph2e:
            saveas_str = self.globalConfig.ism_toa_e + band
            writeToa(self.outdir, saveas_str, toa)

        # PRNU
        # -------------------------------------------------------------------------------
        if self.ismConfig.apply_prnu:

            self.logger.info("EODP-ALG-ISM-2020: PRNU")
            toa = self.prnu(toa, self.ismConfig.kprnu)

            self.logger.debug("TOA [0,0] " + str(toa[0, 0]) + " [e-]")

            if self.ismConfig.save_after_prnu:
                saveas_str = self.globalConfig.ism_toa_prnu + band
                writeToa(self.outdir, saveas_str, toa)

        # Dark-signal
        # -------------------------------------------------------------------------------
        if self.ismConfig.apply_dark_signal:

            self.logger.info("EODP-ALG-ISM-2020: Dark signal")
            toa = self.darkSignal(toa, self.ismConfig.kdsnu, self.ismConfig.T, self.ismConfig.Tref,
                                  self.ismConfig.ds_A_coeff, self.ismConfig.ds_B_coeff)

            self.logger.debug("TOA [0,0] " + str(toa[0, 0]) + " [e-]")

            if self.ismConfig.save_after_ds:
                saveas_str = self.globalConfig.ism_toa_ds + band
                writeToa(self.outdir, saveas_str, toa)

        # Bad/dead pixels
        # -------------------------------------------------------------------------------
        if self.ismConfig.apply_bad_dead:
            self.logger.info("EODP-ALG-ISM-2050: Bad/dead pixels")
            toa = self.badDeadPixels(toa,
                                     self.ismConfig.bad_pix,
                                     self.ismConfig.dead_pix,
                                     self.ismConfig.bad_pix_red,
                                     self.ismConfig.dead_pix_red,
                                     band)

        # Write output TOA
        # -------------------------------------------------------------------------------
        if self.ismConfig.save_detection_stage:
            saveas_str = self.globalConfig.ism_toa_detection + band

            writeToa(self.outdir, saveas_str, toa)

            title_str = 'TOA after the detection phase [e-]'
            xlabel_str = 'ACT'
            ylabel_str = 'ALT'
            plotMat2D(toa, title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

            idalt = int(toa.shape[0] / 2)
            saveas_str = saveas_str + '_alt' + str(idalt)
            plotF([], toa[idalt, :], title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

        return toa

    def irrad2Phot(self, toa, area_pix, tint, wv, band):
        """
        Conversion of the input Irradiances to Photons
        :param toa: input TOA in irradiances [mW/m2]
        :param area_pix: Pixel area [m2]
        :param tint: Integration time [s]
        :param wv: Central wavelength of the band [m]
        :return: Toa in photons
        """


        h = scipy.constants.h
        c = scipy.constants.c

        # Convert toa to W/m2
        toa_watts = toa / 1000
        # Calculate toa in photons
        toa_ph = toa_watts * area_pix * tint / (h * c / wv)
        # Irradiation to photons
        I2P = toa_ph / toa * 1000

        output_directory = 'C:/Users/diego/PycharmProjects/EODP_TER_2021/EODP-TS-ISM/myoutput/TXT_OUT'
        file_path = os.path.join(output_directory, 'I2P.txt')

        if band == 'VNIR-0':
            file = open(file_path, 'w')
            file.truncate(0)
            file.close()

        with open(file_path, 'a') as the_file:
            the_file.write(band + '\n')
            the_file.write('Irradiation to Photon factor ' + '=' + str(I2P) + '\n')

        return toa_ph
        # TODO
    def phot2Electr(self, toa, QE, FWC, band):
        """
        Conversion of photons to electrons
        :param toa: input TOA in photons [ph]
        :param QE: Quantum efficiency [e-/ph]
        :return: toa in electrons
        """
        # TODO
        toae = toa * QE

        output_directory = 'C:/Users/diego/PycharmProjects/EODP_TER_2021/EODP-TS-ISM/myoutput/TXT_OUT'
        qe_file_path = os.path.join(output_directory, 'QE.txt')

        if band == 'VNIR-0':
            with open(qe_file_path, 'w') as file:
                file.truncate(0)

        with open(qe_file_path, 'a') as the_file:
            the_file.write(f'{band}\n')
            the_file.write(f'QE = {QE}\n')
            counter = 0
            for i in range(toae.shape[0]):
                for j in range(toae.shape[1]):
                    if toae[i, j] > FWC:
                        toae[i, j] = FWC
                        counter += 1

            per_saturated = (counter / (toae.shape[0] * toae.shape[1])) * 100
            the_file.write(f'\nPercentage of saturated pixels = {per_saturated}%\n')

        return toae

    def badDeadPixels(self, toa, bad_pix, dead_pix, bad_pix_red, dead_pix_red, band):
        """
        Bad and dead pixels simulation
        :param toa: input toa in [e-]
        :param bad_pix: Percentage of bad pixels in the CCD [%]
        :param dead_pix: Percentage of dead pixels in the CCD [%]
        :param bad_pix_red: Reduction in the quantum efficiency for the bad pixels [-, over 1]
        :param dead_pix_red: Reduction in the quantum efficiency for the dead pixels [-, over 1]
        :return: toa in e- including bad & dead pixels
        """
        # TODO
        # toa[:, 5] = toa[:, 5] * (1 - bad_pix_red)

        output_directory = 'C:/Users/diego/PycharmProjects/EODP_TER_2021/EODP-TS-ISM/myoutput/TXT_OUT'

        qfile_pathBP = os.path.join(output_directory, 'BadP.txt')
        qfile_pathDP = os.path.join(output_directory, 'DeadP.txt')

        if band == 'VNIR-0':
            with open(qfile_pathBP, 'w+') as file:
                file.truncate(0)

            with open(qfile_pathDP, 'w+') as file:
                file.truncate(0)

        with open(qfile_pathDP, 'a') as the_file:
            the_file.write(f'{band}\n')
            num_pixels_act = toa.shape[1]
            tot_DP = int(dead_pix / 100 * num_pixels_act)  # Total dead pixels
            if tot_DP == 0:
                the_file.write('No dead pixels index found \n')

            else:
                step_ded = int(num_pixels_act / tot_DP)
                i_DP = np.arange(0, num_pixels_act, step_ded)
                the_file.write('\n Index of bad pixels:\n')

                for ii in range(i_DP.shape[0]):
                    # Computation of the TOA
                    toa[:, i_DP[ii]] = toa[:, i_DP[ii]] * (1 - dead_pix_red)
                    the_file.write(f'index {ii} = {i_DP[ii]}\n')
                    the_file.write('\n')

            with open(qfile_pathBP, 'a') as file:
                file.write(f'{band}\n')
                num_pixels_act = toa.shape[1]
                tot_BP = int(bad_pix / 100 * num_pixels_act)  # Total bad pixels
                if tot_BP == 0:
                    file.write('No bad pixels index found \n')

                else:
                    step_bad = int(num_pixels_act / tot_BP)
                    i_BP = np.arange(5, num_pixels_act, step_bad)
                    file.write('\n Index of bad pixels:\n')

                    for ii in range(i_BP.shape[0]):
                        # Computation of the TOA
                        toa[:, i_BP[ii]] = toa[:, i_BP[ii]] * (1 - bad_pix_red)
                        file.write(f'index {ii} = {i_BP[ii]}\n')
                        file.write('\n')
        return toa

    def prnu(self, toa, kprnu):
        """
        Adding the PRNU effect
        :param toa: TOA pre-PRNU [e-]
        :param kprnu: multiplicative factor to the standard normal deviation for the PRNU
        :return: TOA after adding PRNU [e-]
        """
        # TODO
        normal = np.random.normal(0., 1., toa.shape[1])
        for act in range(toa.shape[1]):
            toa[:, act] = toa[:, act] * (1 + kprnu * normal[act])

        return toa

    def darkSignal(self, toa, kdsnu, T, Tref, ds_A_coeff, ds_B_coeff):
        """
        Dark signal simulation
        :param toa: TOA in [e-]
        :param kdsnu: multiplicative factor to the standard normal deviation for the DSNU
        :param T: Temperature of the system
        :param Tref: Reference temperature of the system
        :param ds_A_coeff: Empirical parameter of the model 7.87 e-
        :param ds_B_coeff: Empirical parameter of the model 6040 K
        :return: TOA in [e-] with dark signal
        """
        # TODO
        DSNU = np.abs(np.random.normal(0, 1, toa.shape[1]) * kdsnu)
        SD = ds_A_coeff * (T / Tref) ** 3 * np.exp(-ds_B_coeff * (1 / T - 1 / Tref))
        DS = SD * (1 + DSNU)

        for act in range(toa.shape[1]):
            toa[:, act] = toa[:, act] + DS[act]

        return toa
