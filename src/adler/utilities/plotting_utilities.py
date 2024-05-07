import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def plot_phasecurve(
    planetoid,
    filt_list=["r"],
    x_plot="phaseAngle",
    y_plot="reduced_mag",
    xerr_plot="magErr",
    fig=None,
    label_list=None,
    col_list=None,
):
    """Make an errorbar scatter plot of reduced magnitude against phase angle to show the phase curve of an Adler object.

    planetoid: AdlerPlanetoid
        AdlerPlanetoid object containing the observational data to be plotted
    filt_list: list
        List of filters to be plotted
    x_plot: str
        Name of the AdlerPlanetoid attribute to be plotted on the x axis
    y_plot: str
        Name of the AdlerPlanetoid attribute to be plotted on the y axis
    xerr_plot: str
        Name of the AdlerPlanetoid attribute for the x axis uncertainties
    fig: matplotlib.figure.Figure
        Optional, pass an existing figure object to be added to
    label_list: list
        Optional, labels for errorbar plot elements
    col_list: list
        Optional, colors for errorbar scatter points

    Returns
    -----------
    fig: matplotlib.figure.Figure
        The figure object which can be manipulated further if required

    """

    if fig:
        # use the figure object that was passed
        ax1 = fig.axes[0]
    else:
        # set up a new figure object
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0, 0])
        ax1.invert_yaxis()
        ax1.set_xlabel(x_plot)
        ax1.set_ylabel(y_plot)

    for i, filt in enumerate(filt_list):
        # get the object data
        obs_filt = planetoid.observations_in_filter(filt)
        x = getattr(obs_filt, x_plot)
        y = getattr(obs_filt, y_plot)
        xerr = getattr(obs_filt, xerr_plot)

        # label the errorbars?
        if label_list is not None:
            l = label_list[i]
        else:
            l = None

        # select colours?
        if col_list is not None:
            c = col_list[i]
        else:
            c = None

        # plot the errorbars
        ax1.errorbar(x, y, xerr, color=c, fmt="o", label=l)

    return fig
