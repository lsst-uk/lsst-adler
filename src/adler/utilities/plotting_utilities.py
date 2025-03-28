import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


def plot_errorbar(
    planetoid,
    filt_list=[],
    x_plot="phaseAngle",
    y_plot="reduced_mag",
    yerr_plot="magErr",
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
    yerr_plot: str
        Name of the AdlerPlanetoid attribute for the y axis uncertainties
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

        if yerr_plot is not None:
            # plot the errorbars
            yerr = getattr(obs_filt, yerr_plot)
            s1 = ax1.errorbar(x, y, yerr, color=c, fmt="o", label=l)  # TODO: get cmap working with errorbar
        else:
            # just plot scatter
            s1 = ax1.scatter(x, y, c=c, label=l)

        if c_plot is not None:
            cbar1 = plt.colorbar(s1)

    # save the figure?
    if filename:
        fig.savefig(filename, facecolor="w", transparent=True, bbox_inches="tight")

    return fig


def plot_phasecurve(
    pc,
    x=np.radians(np.linspace(0, 30)),
    x_plot="phaseAngle",
    y_plot="reduced_mag",
    fig=None,
    label=None,
    col=None,
    alpha=None,
    filename=None,
    bounds_std=None,
):
    """Display the possible range of model values for a PhaseCurve object with uncertainty.

    pc: PhaseCurve
        PhaseCurve object containing the PhaseCurve model parameters
    filt_list: list
        List of filters to be plotted
    x : array
           Array of x values (phase angle) at which to evaluate the PhaseCurve model
    x_plot: str
        x axis label
    y_plot: str
        y axis label
    fig: matplotlib.figure.Figure
        Optional, pass an existing figure object to be added to
    label: str
        Optional, label for the plot element
    col: str
        Optional, color for the plot element
    alpha: float
        Optional, transparency for the plot element
    filename: str
        Optional, if provided save the figure with this filename
    bounds_std : float
           Optional, if included plot the bounds of the PhaseCurve function. Number of standard deviations we want to consider when determining the minimum/maximum parameter values, where we assume that the parameter uncertainty is representative of a 1 sigma std in a gaussian distribution.

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

    if bounds_std:
        pc_bounds = pc.ReducedMagBounds(x, std=bounds_std)
        ax1.fill_between(x, pc_bounds["mag_min"], pc_bounds["mag_max"], alpha=alpha, color=col, label=label)
    else:
        y = pc.ReducedMag(x)
        ax1.plot(x, y, alpha=alpha, color=col, label=label)

    # save the figure?
    if filename:
        fig.savefig(filename, facecolor="w", transparent=True, bbox_inches="tight")

    return fig
