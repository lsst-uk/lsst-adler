import pytest
import matplotlib
import os
from numpy.testing import assert_array_almost_equal

from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.utilities.tests_utilities import get_test_data_filepath
from adler.utilities.plotting_utilities import plot_errorbar

# set up test planetoid object
ssoid = 8268570668335894776
test_db_path = get_test_data_filepath("testing_database.db")
test_planetoid = AdlerPlanetoid.construct_from_SQL(
    ssoid, test_db_path, filter_list=["u", "g", "r", "i", "z", "y"]
)


def test_plot_errorbar_return():

    # make the fig object
    fig = plot_errorbar(test_planetoid, filt_list=["r"])
    print(isinstance(fig, matplotlib.figure.Figure))

    # check that fig is of the correct type
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_errorbar_file(tmp_path):

    # define a pytest tmp_path
    filename = "{}/test_file.png".format(tmp_path)

    # make the fig object
    fig = plot_errorbar(test_planetoid, filt_list=["r"], filename=filename)
    print(os.path.isfile(filename))

    # check that the file exists
    assert os.path.isfile(filename)


def test_plot_errorbar_xy_label():

    # make the fig object
    x_plot = "midPointMjdTai"
    y_plot = "mag"
    fig = plot_errorbar(test_planetoid, filt_list=["r"], y_plot=y_plot, x_plot=x_plot)

    # check the labels have been set
    assert fig.gca().get_xlabel() == x_plot
    assert fig.gca().get_ylabel() == y_plot

    # get the object data
    obs_filt = test_planetoid.observations_in_filter("r")
    x = getattr(obs_filt, x_plot)
    y = getattr(obs_filt, y_plot)

    # check the correct data has been plotted
    x_fig = fig.gca().lines[0].get_xdata()
    y_fig = fig.gca().lines[0].get_ydata()
    assert_array_almost_equal(x_fig, x)
    assert_array_almost_equal(y_fig, y)


def test_plot_errorbar_fig():

    # make the fig object
    fig = plot_errorbar(test_planetoid, filt_list=["r"])
    # add to the first fig object
    fig = plot_errorbar(test_planetoid, filt_list=["g"], fig=fig)
    print(fig.gca().lines)
    # check that there are two sets of data plotted
    assert len(fig.gca().lines) == 2


def test_plot_errorbar_legend():

    # make the fig object
    fig = plot_errorbar(test_planetoid, filt_list=["r", "g"], label_list=["r", "g"])

    # make the legend
    ax1 = fig.gca()
    ax1.legend()

    # check the legend labels
    assert ax1.get_legend_handles_labels()[1] == ["r", "g"]


# TODO: test no yerr_plot
# TODO: test col_list
# TODO: test plot_phasecurve

if __name__ == "__main__":
    print(test_planetoid.__dict__)

    # test_plot_errorbar_return()
    # test_plot_errorbar_file(".")
    # test_plot_errorbar_xy_label()
    # test_plot_errorbar_fig()
    test_plot_errorbar_legend()
