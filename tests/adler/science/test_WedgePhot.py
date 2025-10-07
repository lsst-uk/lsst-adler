import numpy as np
import pandas as pd
import subprocess
import os
from adler.science.WedgePhot import WedgePhot
from numpy.testing import assert_array_equal

# test image of the DART mission impact on Didymos system, from ZTF
dir_path = os.path.dirname(os.path.realpath(__file__))  # dir path of this file
test_dir_path = "/".join(dir_path.split("/")[:-2])
infits = "{}/data/ztf_Didymos-system-barycenter20065803_20221201402616_000567_zr_c10_o_q1_scimrefdiffimg.fits.fz".format(
    test_dir_path
)
inhdu = 1
N_wedge = 10

wp = WedgePhot(
    fits_file=infits,
    i_hdu=inhdu,
    N_wedge=N_wedge,
    measure="sum,mean,median,sigclip-mean,sigclip-std",
)


def test_WedgePhot_init():
    """Test intialising the wedge photometry class."""

    # check azimuthal bins are correct
    assert_array_equal(
        wp.az, np.array([0.0, 36.0, 72.0, 108.0, 144.0, 180.0, 216.0, 252.0, 288.0, 324.0, 360.0])
    )


def test_astscript_radial_profile():
    """Test that the gnuastro astscript_radial_profile can be executed"""

    cmd = "which {}".format(wp.ast_radial_profile)
    print(cmd)
    result = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = result.communicate()
    stdout = stdout.decode("utf-8")
    stderr = stderr.decode("utf-8")
    assert "astscript-radial-profile" in stdout


def test_run_wedge_phot():
    """Test that the radial profile is run for all bins and results are compiled"""

    wp_results = wp.run_wedge_phot()  # TODO: the subprocess cat call doesn't seem to work here?
    df = wp_results[0]["data"]
    assert len(wp_results) == N_wedge
    assert isinstance(df, pd.DataFrame)
    for x in ["RADIUS", "SUM", "MEAN", "MEDIAN", "SIGCLIP_MEAN", "SIGCLIP_STD"]:
        assert x in list(df)


if __name__ == "__main__":
    import time

    start = time.time()
    test_run_wedge_phot()
    end = time.time()
    length = end - start
    print("time = ", length, "seconds")
