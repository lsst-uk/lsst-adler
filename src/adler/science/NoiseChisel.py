import numpy as np
import pandas as pd
import sys
import os
import subprocess


class NoiseChisel:
    """
    Class for performing NoiseChisel source detection/segmentation analysis on a given fits image using gnuastro astnoisechisel.

    Attributes
    -----------
    fits_file : str
        Filename of fits file to be analysed.
    i_hdu: int
        Fits file HDU index containing the image to be analysed.
    """

    def __init__(
        self,
        fits_file,
        i_hdu,
        fits_file_suffix=".fits.fz",
        out_dir=".",
    ):

        self.fits_file = fits_file
        self.i_hdu = i_hdu

        self.fits_file_suffix = fits_file_suffix
        self.out_dir = out_dir
        self.file_suffix_out = "_detected.fits"  # output file suffix if noisechisel is successful
        self.file_suffix_check = "_detcheck.fits"  # output diagnostic file suffix

        self.file_root = fits_file.split("/")[-1].split(self.fits_file_suffix)[0]
        self.file_out = self.file_root + self.file_suffix_out
        self.file_check = self.file_root + self.file_suffix_check
        print("expected files:\n{}\n{}".format(self.file_out, self.file_check))

    def run_noise_chisel(self):

        ast_cmd = """astnoisechisel {} --hdu={} --checkdetection --continueaftercheck""".format(
            self.fits_file,
            self.i_hdu,
        )

        result = subprocess.run(ast_cmd, shell=True, capture_output=True, text=True)

        out = result.stdout
        err = result.stderr
        print(out)
        print(err)
