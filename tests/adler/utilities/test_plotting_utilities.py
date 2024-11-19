import pytest
import matplotlib

from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.utilities.tests_utilities import get_test_data_filepath
from adler.utilities.plotting_utilities import plot_errorbar

# set up test planetoid object
ssoid = 8268570668335894776
test_db_path = get_test_data_filepath("testing_database.db")
test_planetoid = AdlerPlanetoid.construct_from_SQL(
    ssoid, test_db_path, filter_list=["u", "g", "r", "i", "z", "y"]
)


def test_plot_errorbar_return(planetoid=test_planetoid, filt_list=["r"]):

    # make the fig object
    fig = plot_errorbar(planetoid, filt_list=filt_list)
    print(isinstance(fig, matplotlib.figure.Figure))

    # check that fig is of the correct type
    assert isinstance(fig, matplotlib.figure.Figure)


if __name__ == "__main__":
    print(test_planetoid.__dict__)

    test_plot_errorbar_return()
