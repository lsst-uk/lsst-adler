import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from adler.science.PhaseCurve import PhaseCurve
from adler.utilities.science_utilities import get_df_obs_filt


def col_obs_ref(
    planetoid,
    adler_cols,
    filt_obs="g",
    filt_ref="r",
    N_ref=3,
    x_col="midPointMjdTai",
    y_col="AbsMag",
    yerr_col="magErr",
    x1=None,
    x2=None,
    plot_dir=None,
):
    """A function to calculate the colour of an Adler planetoid object.
    An observation in a given filter (filt_obs) is compared to previous observation in a reference filter (filt_ref).
    By setting N_ref one can control how many reference observations to include in the colour calculation.
    Note that the observations are considered to be in reduced magnitudes, hence why an adlerData object with phase curve models is passed.

    Parameters
    -----------
    planetoid : object
        Adler planetoid object

    adler_cols : object
        AdlerData object

    filt_obs: str
        Filter name of the new observation for which to calculate a filt_obs - filt_ref colour

    filt_ref:str
        Filter name of the reference observations when calculating the filt_obs - filt_ref colour

    N_ref: int
        Number of reference observations to use when calculated the reference absolute magnitude.
        Set to 1 to use only the most recent of ref observations.
        Set to None to use all past ref obs.

    x_col: str
        Time (or phase angle) column name

    y_col:str
        Magnitude column name (e.g. AbsMag or reduced_mag)

    yerr_col: str
        Magnitude uncertainty column name

    x1: float
        Lower limit on x_col values (>=)

    x2: float
        Upper limit on x_col values (<=)

    plot_dir: str
        Directory in which to save the colour plot

    Returns
    ----------

    col_dict : dict
        Dictionary containing the colour information for the most recent obs in filt_obs

    """
    # define fields to be recorded as a colour dictionary
    colour = "{}-{}".format(filt_obs, filt_ref)
    colErr = "{}-{}Err".format(filt_obs, filt_ref)
    delta_t_col = "delta_t_{}".format(colour)
    y_ref_col = "{}_{}".format(y_col, filt_ref)
    x1_ref_col = "{}1_{}".format(x_col, filt_ref)
    x2_ref_col = "{}2_{}".format(x_col, filt_ref)

    # get the stored AdlerData parameters for this filter # TODO: do this bit better
    ad_obs = adler_cols.get_phase_parameters_in_filter(filt_obs, "HG12_Pen16")
    # make the PhaseCurve object from AdlerData
    pc_obs = PhaseCurve().InitModelDict(ad_obs.__dict__)
    ad_ref = adler_cols.get_phase_parameters_in_filter(filt_ref, "HG12_Pen16")
    pc_ref = PhaseCurve().InitModelDict(ad_ref.__dict__)
    df_obs = get_df_obs_filt(planetoid, filt_obs, x1=x1, x2=x2, col_list=[y_col, yerr_col], pc_model=pc_obs)
    df_obs_ref = get_df_obs_filt(
        planetoid, filt_ref, x1=x1, x2=x2, col_list=[y_col, yerr_col], pc_model=pc_ref
    )

    # select the values of the new observation in the selected filter
    i = -1
    x_obs = df_obs.iloc[i][x_col]
    y_obs = df_obs.iloc[i][y_col]
    yerr_obs = df_obs.iloc[i][yerr_col]

    # select observations in the reference filter from before the new obs
    ref_mask = df_obs_ref[x_col] < x_obs

    # set the number of ref obs to use
    if N_ref is None:
        _N_ref = len(df_obs_ref[ref_mask])  # use all available ref obs
    else:
        _N_ref = N_ref

    # select only the N_ref ref obs for comparison
    _df_obs_ref = df_obs_ref[ref_mask].iloc[-_N_ref:]
    if len(_df_obs_ref) == 0:
        print("no reference observations")  # TODO: add proper error handling and logging here
        return df_obs

    # determine reference observation values
    y_ref = np.mean(_df_obs_ref[y_col])  # TODO: add option to choose statistic, e.g. mean or median?
    yerr_ref = np.std(_df_obs_ref[y_col])  # TODO: propagate ref uncertainty properly
    # determine the ref obs time range
    x1_ref = np.array(_df_obs_ref[x_col])[0]
    x2_ref = np.array(_df_obs_ref[x_col])[-1]

    # Create the colour dict
    col_dict = {}
    col_dict[colour] = y_obs - y_ref
    col_dict[delta_t_col] = x_obs - x2_ref
    col_dict[colErr] = np.sqrt((yerr_obs**2.0) + (yerr_ref**2.0))
    col_dict[y_ref_col] = y_ref
    col_dict[x1_ref_col] = x1_ref
    col_dict[x2_ref_col] = x2_ref

    # TODO:
    # could also record phase angle diff and phase curve residual?
    # need to test error case where there are no r filter obs yet

    # TODO: add a plotting option?
    if plot_dir:
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0, 0])

        ax1.errorbar(df_obs[x_col], df_obs[y_col], df_obs[yerr_col], fmt="o", label=filt_obs)
        ax1.errorbar(df_obs_ref[x_col], df_obs_ref[y_col], df_obs_ref[yerr_col], fmt="o", label=filt_ref)

        # plot some lines to show the colour and mean reference
        ax1.vlines(x_obs, y_obs, col_dict[y_ref_col], color="k", ls=":")
        ax1.hlines(col_dict[y_ref_col], col_dict[x1_ref_col], col_dict[x2_ref_col], color="k", ls="--")

        ax1.set_xlabel(x_col)
        ax1.set_ylabel(y_col)
        ax1.legend()
        ax1.invert_yaxis()

        fname = "{}/colour_plot_{}_{}-{}_{}.png".format(
            plot_dir, planetoid.ssObjectId, filt_obs, filt_ref, int(x_obs)
        )
        print("Save figure: {}".format(fname))
        plt.savefig(fname, facecolor="w", transparent=True, bbox_inches="tight")

        # plt.show() # TODO: add option to display figure, or to return the fig object?
        plt.close()

    return col_dict
