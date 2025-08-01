import numpy as np
import pandas as pd
import sys
import os
import subprocess
from io import StringIO


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
    N_wedge: int
        Number of azimuthal bins to calculate radial profiles for.
    measure: str
        Statistic(s) to be calculated on each radial profile, can be individual or a list.
        E.g.: sum,mean,median,sigclip-mean,sigclip-std
    ap_rad_out: float
        (Optional) Maximum radius (pixels) to calculate radial profile over.
    """

    def __init__(
        self,
        fits_file,
        i_hdu,
        N_wedge=10,
        measure="sum",
        ap_rad_out=None,
    ):

        self.fits_file = fits_file
        self.i_hdu = i_hdu
        self.measure = measure
        self.N_wedge = N_wedge
        self.ap_rad_out = ap_rad_out

        # define the current working directory, i.e. where subprocess is running
        # this is needed for working with the output files
        self.cwd = (
            os.getcwd()
        )  # TODO: include option to set cwd? use cd <cwd> in subprocess to control where everything runs?

        # set up azimuthal bin edges (degrees)
        self.az = np.linspace(0, 360, self.N_wedge + 1)

        # gnuastro setup, find the path to astscript-radial-profile executable
        ENVBIN = sys.exec_prefix
        self.ast_radial_profile = os.path.join(ENVBIN, "bin", "astscript-radial-profile")
        self.col_fmt = (
            "# Column"  # string to identify column names when astscript-radial-profile is printed out
        )

    def astscript_radial_profile(
        self,
        az_min,
        az_max,
        out_file,
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
            Maximum azimuthal bin edge
        out_file: str
            File to store the results of astscript-radial-profile (deleted after results are extracted).

        Returns
        ----------
        df_results: DataFrame
            A Pandas DataFrame containing the radial profile statistics for the given azimuthal bin.
        """

        # Construct the gnuastro command
        ast_cmd = "{} -q -h{} {} -a {},{} --measure={} -o {}".format(
            self.ast_radial_profile, self.i_hdu, self.fits_file, az_min, az_max, self.measure, out_file
        )

        # Set the maximum aperture radius if required
        if self.ap_rad_out is not None:
            ast_cmd += " -R {};".format(self.ap_rad_out)
        else:
            ast_cmd += ";"

        # To get the output we use "cat" to print the out_file contents
        # We also remove the temporary file # TODO: option to keep output files?
        ast_cmd += " cat {};\
                rm {};".format(
            out_file, out_file
        )

        print(ast_cmd)

        result = subprocess.Popen("pwd", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.communicate()[0].decode("utf-8"))

        # run the command
        result = subprocess.Popen(ast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # print(result)

        # get the results as a dataframe
        data = result.communicate()[0].decode("utf-8")
        # print(data)
        col_names = [x.split(": ")[-1].split(" ")[0] for x in data.split("\n") if self.col_fmt in x]
        # print(col_names)
        df_results = pd.read_csv(StringIO(data), sep="\s+", names=col_names, comment="#")

        # TODO: properly log the output from astscript-radial-profile

        return df_results

    def run_wedge_phot(self):
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
            outfile = "{}/wedge_out_{}.txt".format(self.cwd, i)

            # run radial profile for a single bin
            df = self.astscript_radial_profile(az_min, az_max, outfile)

            # store the results in a dict
            wp_results[i] = {}
            wp_results[i]["az_min"] = az_min
            wp_results[i]["az_max"] = az_max
            wp_results[i]["data"] = df

        return wp_results
