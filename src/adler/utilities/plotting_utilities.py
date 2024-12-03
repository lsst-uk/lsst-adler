import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def plot_errorbar(
    planetoid,
    filt_list=[],
    x_plot="phaseAngle",
    y_plot="reduced_mag",
    xerr_plot="magErr",
    c_plot=None,
    fig=None,
    label_list=None,
    col_list=None,
    filename=None,
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
    c_plot: str
        Name of the AdlerPlanetoid attribute used for the colour scale
    fig: matplotlib.figure.Figure
        Optional, pass an existing figure object to be added to
    label_list: list
        Optional, labels for errorbar plot elements
    col_list: list
        Optional, colors for errorbar scatter points
    filename: str
        Optional, if provided save the figure with this filename

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

        # label the errorbars?
        if label_list is not None:
            l = label_list[i]
        else:
            l = None

        # select colours?
        if (col_list is not None) and (c_plot is None):
            c = col_list[i]
        elif c_plot is not None:
            c = getattr(obs_filt, c_plot)
            # TODO: use kwargs to pass the cmap?
        else:
            c = None

        if xerr_plot is not None:
            # plot the errorbars
            xerr = getattr(obs_filt, xerr_plot)
            s1 = ax1.errorbar(x, y, xerr, color=c, fmt="o", label=l)  # TODO: get cmap working with errorbar
        else:
            # just plot scatter
            s1 = ax1.scatter(x, y, c=c, label=l)

        if c_plot is not None:
            cbar1 = plt.colorbar(s1)

    # save the figure?
    if filename:
        fig.savefig(filename, facecolor="w", transparent=True, bbox_inches="tight")

    return fig
