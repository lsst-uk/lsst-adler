import os
from astropy.table import Table
from astropy.io import fits

import adler.utilities.science_utilities as sci_utils


class NoiseChisel:
    """
    Class for performing NoiseChisel source detection/segmentation analysis on a given fits image using gnuastro astnoisechisel.

    Attributes
    -----------
    fits_file : str
        Filename of fits file to be analysed. Assumed to be like *.fit*
    i_hdu: int
        Fits file HDU index containing the image to be analysed.
    out_dir: str
        Path to directory for saving files. Directory will be created if necessary.
    """

    def __init__(
        self,
        fits_file,
        i_hdu,
        out_dir=".",
    ):

        self.fits_file = fits_file
        self.i_hdu = i_hdu
        self.out_dir = out_dir
        if not os.path.isdir(self.out_dir):
            os.path.mkdir(self.out_dir)

        # define the expected file names to be created during the noisechisel process
        # TODO: use os path join or similar to make sure paths work
        self.file_suffix_nc = "_detected.fits"  # output file suffix if noisechisel is successful
        self.file_suffix_check = "_detcheck.fits"  # output diagnostic file suffix
        self.file_root = self.out_dir + "/" + fits_file.split("/")[-1].split(".fit")[0]
        self.file_nc = self.file_root + self.file_suffix_nc
        self.file_check = self.file_root + self.file_suffix_check
        self.file_seg = self.file_nc.split(".fits")[0] + "_segmented.fits"
        self.file_cat = self.file_seg.split(".fits")[0] + "_cat.fits"

        # Define the noisechisel command
        self.astnoisechisel = "astnoisechisel"

    def noise_chisel(self, nc_flags="--checkdetection --continueaftercheck", pre_cmd=None):
        """
        Function to invoke the gnuastro noisechisel command. The results are stored in the file that is created (file_nc). This file contains a pixel map of detections, i.e. is a pixel signal or background?

        Returns
        ----------
        file_nc: str
            Name of the noisechisel results file
        """

        ast_cmd = """{} {} --hdu={} {}""".format(
            self.astnoisechisel,
            self.fits_file,
            self.i_hdu,
            nc_flags,
        )

        # add prerequisite commands if necessary
        if pre_cmd is not None:
            ast_cmd = pre_cmd + ast_cmd

        print(ast_cmd)
        out, err = sci_utils.execute_subprocess(ast_cmd)
        print(out)
        print(err)

        return self.file_nc

    def segment_image(self, pre_cmd=None):
        """
        Function to invoke the gnuastro image segmentation routine command. This step is required to separate the noisechisel chisel detections into individual objects/clumps. The results are stored in the file that is created (file_seg).

        Returns
        ----------
        file_seg: str
            Name of the image segmentation results file
        """
        ast_cmd = "astsegment {} --clumpsnthresh=5".format(self.file_nc)

        # add prerequisite commands if necessary
        if pre_cmd is not None:
            ast_cmd = pre_cmd + ast_cmd

        print(ast_cmd)
        out, err = sci_utils.execute_subprocess(ast_cmd)
        print(out)
        print(err)

        return self.file_seg

    def make_catalogue(self, i_cat=1, pre_cmd=None):
        """
        Function to invoke the gnuastro make catalogue command. The results are stored in the file that is created (file_cat)


        Parameters
        -----------
        i_cat : int
            Use either the object (i_cat=1) or the clump (i_cat=2) detections to make the catalogue

        Returns
        ----------
        df_cat: DataFrame
            Dataframe containing the measured properties of the clumps
        """
        # TODO: use a gnuastro conf file to determine which columns are calculated?

        ast_cmd = "astmkcatalog {} --clumpscat --ids -x -y --ra --dec --magnitude --sn --axis-ratio --geo-axis-ratio --geo-position-angle --geo-semi-major --geo-semi-minor --position-angle --semi-major --semi-minor".format(
            self.file_seg
        )

        # add prerequisite commands if necessary
        if pre_cmd is not None:
            ast_cmd = pre_cmd + ast_cmd

        print(ast_cmd)
        out, err = sci_utils.execute_subprocess(ast_cmd)
        print(out)
        print(err)

        hdu_cat = fits.open(self.file_cat)
        # print(hdu_cat.info())
        dat = hdu_cat[i_cat].data
        df_cat = Table(dat).to_pandas()

        return df_cat

    def clean_up(self):
        """
        Function to remove all .fits created as part of the noisechisel, segmentation and catalogue process.
        """

        ast_cmd = "rm {}_detected_*_.fits".format(self.file_root)

        print(ast_cmd)
        out, err = sci_utils.execute_subprocess(ast_cmd)
        print(out)
        print(err)

        return

    def run_noise_chisel(self, conda_start=None, conda_env=None, keep_files=False):
        """
        Wrapper function that calls each step to go from an input image to measurements of detections made by noisechisel.

        Parameters
        -----------
        keep_files : float
            Optional, flag to either remove all noisechisel associated files (by default) or keep them.

        Returns
        ----------
        df_cat: DataFrame
            Dataframe containing the measured properties of the clumps
        """

        # Explicitly start conda and run in the environment if required
        if (conda_start is not None) & (conda_env is not None):
            conda_run = "{}; ".format(conda_start)
        else:
            conda_run = None
        if conda_env is not None:
            conda_run += "conda run -n {} ".format(conda_env)

        self.noise_chisel(pre_cmd=conda_run)
        self.segment_image(pre_cmd=conda_run)
        df_cat = self.make_catalogue(pre_cmd=conda_run)

        if not keep_files:
            self.clean_up()

        return df_cat
