import numpy as np
import pandas as pd
from io import StringIO
from astropy.io import fits

import adler.utilities.science_utilities as sci_utils


class WedgePhot:
    """
    Class for performing wedge photometry on a given fits image using gnuastro astscript-radial-profile.
    An individual call to astscript-radial-profile calculates the radial profile in a given azimuthal angle bin, from the center to edge of the image.
    This class calculates the radial profile in a number of bins, N_wedge, from angles 0 to 360 degrees.

    Attributes
    -----------
    fits_file : str
        Filename of fits file to be analysed.
    i_hdu: int
        Fits file HDU index containing the image to be analysed.
    x,y: float,float
        Pixel coords of the target, i.e. the centre of the radial profile. If not passed the centre of the image is used
    N_wedge: int
        Number of azimuthal bins to calculate radial profiles for.
    measure: str
        Statistic(s) to be calculated on each radial profile, can be individual or a list.
        E.g.: sum,mean,median,sigclip-mean,sigclip-std
    ap_rad_out: float
        (Optional) Maximum radius (pixels) to calculate radial profile over. Will default to the cropped image size.
    """

    def __init__(
        self,
        fits_file,
        i_hdu,
        x=None,
        y=None,
        N_wedge=10,
        measure="sum",
        out_dir=".",
        ap_rad_out=None,
    ):

        self.fits_file = fits_file
        self.i_hdu = i_hdu
        self.measure = measure
        self.N_wedge = N_wedge

        # define the directory to save any output files
        self.out_dir = out_dir

        # set up azimuthal bin edges (degrees)
        self.az = np.linspace(0, 360, self.N_wedge + 1)

        # gnuastro setup, find the path to astscript-radial-profile executable
        # ENVBIN = sys.exec_prefix
        # self.ast_radial_profile = os.path.join(ENVBIN, "bin", "astscript-radial-profile")
        self.ast_radial_profile = "astscript-radial-profile"

        self.col_fmt = (
            "# Column"  # string to identify column names when astscript-radial-profile is printed out
        )

        # Explicitly set x and y and for the astscript-radial-profile --center argument
        if (x is None) | (y is None):
            # use the image centre as default # TODO: check pixel conventions!
            hdu = fits.open(self.fits_file)
            self.y, self.x = np.array(hdu[self.i_hdu].shape) / 2
        else:
            self.x = x
            self.y = y
        # Set ap_rad_out
        if ap_rad_out is None:
            self.ap_rad_out = np.min(
                [self.x, self.y]
            )  # use the minimum image size (from target position) by default
        else:
            self.ap_rad_out = ap_rad_out  # set the radius passed to WedgePhot class

    def astscript_radial_profile(
        self,
        az_min,
        az_max,
        out_file,
        conda_start=None,
        conda_env=None,
        keep_files=False,
        extra_options=None,
    ):
        """
        Get a radial profile in a given azimuthal range using gnuastro astscript-radial-profile.
        https://www.gnu.org/software/gnuastro/manual/html_node/Generate-radial-profile.html

        By default gnuastro assumes the centre of the image.
        This function use subprocess to launch and run a radial profile for a single azimuthal bin.

        Parameters
        -----------
        az_min : float
            Minimum azimuthal bin edge.
        az_max : float
            Maximum azimuthal bin edge.
        out_file: str
            File to store the results of astscript-radial-profile (deleted after results are extracted).
        conda_start: str
            Optional, command to launch conda in the subprocess virtual environment.
            This might be needed when running WedgePhot in a jupyter notebook (e.g. on the RSP conda_start = ". /opt/lsst/software/stack/conda/etc/profile.d/conda.sh")
        conda_env: str
            Optional, name of the conda environment to use in the subprocess virtual environment (TODO: this makes WedgePhot slow).
        keep_files: Boolean
            Optional, flag that can be set to true in order to keep the output file(s).
        extra_options: str
            Optional, pass any additional astscript-radial-profile command line options here (e.g. --keeptmp).
        Returns
        ----------
        df_results: DataFrame
            A Pandas DataFrame containing the radial profile statistics for the given azimuthal bin.
        """

        # Explicitly start conda and run in the environment if required
        if (conda_start is not None) & (conda_env is not None):
            conda_run = "{}; ".format(conda_start)
        else:
            conda_run = ""
        if conda_env is not None:
            conda_run += "conda run -n {} ".format(conda_env)

        # Construct the gnuastro command
        ast_cmd = "{}{} -q -h{} {} -a {},{} --measure={} --center={},{} -o {}".format(
            conda_run,
            self.ast_radial_profile,
            self.i_hdu,
            self.fits_file,
            az_min,
            az_max,
            self.measure,
            self.x,
            self.y,
            out_file,
        )

        # add any additional options
        if extra_options:
            ast_cmd += " " + extra_options

        # Set the maximum aperture radius if required
        if self.ap_rad_out is not None:
            ast_cmd += " -R {};".format(self.ap_rad_out)
        else:
            ast_cmd += ";"

        # To get the output we use "cat" to print the out_file contents
        ast_cmd += " cat {};".format(out_file)
        # We remove the temporary file(s) by default
        # TODO: option to also keep output fits files?
        if not keep_files:
            ast_cmd += " rm {};".format(out_file)

        # run the command
        print(ast_cmd)
        out, err = sci_utils.execute_subprocess(ast_cmd)

        # get the results as a dataframe
        if out != "":
            col_names = [x.split(": ")[-1].split(" ")[0] for x in out.split("\n") if self.col_fmt in x]
            df_results = pd.read_csv(StringIO(out), sep="\s+", names=col_names, comment="#")
        else:
            df_results = None
            print(err)

        # TODO: properly log the output and error from astscript-radial-profile!

        return df_results

    def run_wedge_phot(self, conda_start=None, conda_env=None, keep_files=False, extra_options=None):
        """
        Function to calculate radial profiles across all azimuthal bins and compile results into a dict.

        Returns
        ----------
        wp_results: dict
            Dictionary containing azimuthal bin edges and radial profile statistics for each wedge.
        """

        # initialise the results dict
        wp_results = {}

        # loop over each wedge
        for i in range(len(self.az) - 1):
            az_min = self.az[i]
            az_max = self.az[i + 1]
            outfile = "{}/wedge_out_{}.txt".format(self.out_dir, i)

            # run radial profile for a single bin
            df = self.astscript_radial_profile(
                az_min, az_max, outfile, conda_start, conda_env, keep_files, extra_options
            )

            # store the results in a dict
            wp_results[i] = {}
            wp_results[i]["az_min"] = az_min
            wp_results[i]["az_max"] = az_max
            wp_results[i]["data"] = df

        return wp_results
