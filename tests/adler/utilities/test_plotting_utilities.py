import pytest
import matplotlib

from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid
from adler.utilities.tests_utilities import get_test_data_filepath
import adler.utilities.plotting_utilities as plot_utils

# set up test planetoid object
ssoid = 8268570668335894776
test_db_path = get_test_data_filepath("testing_database.db")
test_planetoid = AdlerPlanetoid.construct_from_SQL(
    ssoid, test_db_path, filter_list=["u", "g", "r", "i", "z", "y"]
)


def test_plot_errorbar_return(planetoid=test_planetoid, filt_list=["r"]):
    fig = plot_utils.plot_errorbar(planetoid, filt_list=filt_list)
    assert type(fig) == matplotlib.figure.Figure
