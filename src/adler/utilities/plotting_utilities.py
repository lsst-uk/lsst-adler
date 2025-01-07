import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u

from adler.science.PhaseCurve import PhaseCurve


def plot_errorbar(
    planetoid,
    filt_list=[],
    x_plot="phaseAngle",
    y_plot="reduced_mag",
    yerr_plot="magErr",
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
        Name of the AdlerPlanetoid attribute for the x axis uncertainties
    fig: matplotlib.figure.Figure
        Optional, pass an existing figure object to be added to
    label_list: list
        Optional, labels for errorbar plot elements. The user can add the legend manually after the fact in case additional elements are added to the figure.
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
        yerr = getattr(obs_filt, yerr_plot)

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
        ax1.errorbar(x, y, yerr, color=c, fmt="o", label=l)

    # save the figure?
    if filename:
        fig.savefig(filename, facecolor="w", transparent=True, bbox_inches="tight")

    return fig


def plot_phasecurve(
    adler_data,
    filt_list=[],
    x_plot="phaseAngle",
    y_plot="reduced_mag",
    model_name="HG12_Pen16",
    fig=None,
    label_list=None,
    col_list=None,
    filename=None,
):
    """Make an errorbar scatter plot of reduced magnitude against phase angle to show the phase curve of an Adler object.

    adler_data: AdlerData
        AdlerData object containing the model phasecurve parameters to be plotted
    filt_list: list
        List of filters to be plotted
    x_plot: str
        Name of the AdlerPlanetoid attribute to be plotted on the x axis
    y_plot: str
        Name of the AdlerPlanetoid attribute to be plotted on the y axis
    model_name: str
        Name of the phasecurve model in AdlerData
    fig: matplotlib.figure.Figure
        Optional, pass an existing figure object to be added to
    label_list: list
        Optional, labels for errorbar plot elements. The user can add the legend manually after the fact in case additional elements are added to the figure.
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
        # Convert the AdlerData phasecurve parameters into a PhaseCurve object
        ad = adler_data.get_phase_parameters_in_filter(filt, model_name)
        pc = PhaseCurve().InitModelDict(ad.__dict__)
        print(pc.__dict__)

        # define the x values (degrees)
        # go from zero to the maximum phaseAngle
        x_max = (
            adler_data.get_phase_parameters_in_filter(filt, model_name).phaseAngle_min
            + adler_data.get_phase_parameters_in_filter(filt, model_name).phaseAngle_range
        )
        x = np.linspace(0, x_max)

        # calculate y from the PhaseCurve model
        # TODO: how to deal with units?
        if hasattr(pc.H, "unit"):
            if (pc.H.unit is None) or (
                pc.H.unit == ""
            ):  # if there are no units (or weird blank units?) use radians
                print("no H units")
                y = pc.ReducedMag(np.radians(x))
            elif pc.H.unit == u.mag:  # if the units of H are mag, then ensure phaseAngles are degrees
                print("H has mag units")
                y = pc.ReducedMag(x * u.deg)
            else:  # if H has weird units, make it unitless and pass phaseAngle in radians
                print("H has weird units")
                pc.H = pc.H.quantity
                y = pc.ReducedMag(np.radians(x))
        else:
            print("no H units")
            y = pc.ReducedMag(np.radians(x))

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
        ax1.plot(x, y, color=c, label=l)

    # save the figure?
    if filename:
        fig.savefig(filename, facecolor="w", transparent=True, bbox_inches="tight")

    return fig
