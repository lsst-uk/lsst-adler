import pandas as pd
import os
import subprocess
from adler.science.NoiseChisel import NoiseChisel

# test image of the DART mission impact on Didymos system, from ZTF
dir_path = os.path.dirname(os.path.realpath(__file__))  # dir path of this file
test_dir_path = "/".join(dir_path.split("/")[:-2])
infits = "{}/data/ztf_Didymos-system-barycenter20065803_size256_20221201402616_000567_zr_c10_o_q1_scimrefdiffimg.fits.fz".format(
    test_dir_path
)
inhdu = 1

nc = NoiseChisel(
    fits_file=infits,
    i_hdu=inhdu,
)


def test_astnoisechisel():
    """Test that the gnuastro astnoisechisel command can be executed"""

    cmd = "which {}".format(nc.astnoisechisel)
    print(cmd)
    result = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = result.communicate()
    stdout = stdout.decode("utf-8")
    stderr = stderr.decode("utf-8")
    assert "astnoisechisel" in stdout


def test_run_noise_chisel():
    """Test that the noisechisel routine runs and that a catalogue of detections is produced"""
    df_cat = nc.run_noise_chisel()

    assert isinstance(df_cat, pd.DataFrame)  # check that a dataframe is returned
    assert len(df_cat) > 0  # check that the dataframe has at least one row
