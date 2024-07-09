import numpy as np
from adler.science.PhaseCurve import PhaseCurve
from adler.utilities.science_utilities import get_df_obs_filt

# TODO: add colour estimate function


def col_obs_ref(
    planetoid,
    adler_cols,
    filt_obs="g",
    filt_ref="r",
    N_mean=3,
    x_col="midPointMjdTai",
    y_col="AbsMag",
    yerr_col="magErr",
    x1=None,
    x2=None,
):
    """

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

    N_mean: int
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
    new_obs_cols = [colour, colErr, delta_t_col, y_ref_col, x1_ref_col, x2_ref_col]

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
    if N_mean is None:
        _N_mean = len(df_obs_ref[ref_mask])  # use all available ref obs
    else:
        _N_mean = N_mean

    # select only the N_mean ref obs for comparison
    _df_obs_ref = df_obs_ref[ref_mask].iloc[-_N_mean:]
    if len(_df_obs_ref) == 0:
        print("no reference observations")  # TODO: add proper error handling and logging here
        return df_obs

    # determine reference observation values
    y_ref = np.mean(_df_obs_ref[y_col])
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

    return col_dict
