import numpy as np
import pandas as pd
from astropy.stats import sigma_clip as astropy_sigma_clip
import astropy.units as u


def outlier_diff(new_res, diff_cut=1.0):
    """Test whether new data point(s) is an outlier compared to the model by considering the residual difference.

    Parameters
    -----------

    new_res: array
        The residuals of the new data points compared to the model
    diff_cut: float
        The threshold difference value for outlier detection.

    Returns
    -----------
    outlier_flag : array
       Array of flag indicating if data point is an outlier (True)

    """

    if not isinstance(new_res, np.ndarray):
        new_res = np.array(new_res, ndmin=1)

    outlier_flag = np.array([False] * len(new_res))
    outlier_flag[np.abs(new_res) >= diff_cut] = True

    return outlier_flag


def outlier_std(new_res, data_res, std_cut=3.0):
    """Test whether new data point(s) is an outlier compared to the model by considering the standard deviation of the residuals.

    Parameters
    -----------

    new_res: array
        The residuals of the new data point(s) compared to the model
    data_res: array
        The residuals of the data compared to the model.
    std_cut: float
       The threshold standard deviation for outlier detection.

    Returns
    -----------
    outlier_flag : array
       Array of flag indicating if data point is an outlier (True)

    """

    # calculate the standard deviation of the data - model residuals
    data_std = np.std(data_res)

    if not isinstance(new_res, np.ndarray):
        new_res = np.array(new_res, ndmin=1)

    outlier_flag = np.array([False] * len(new_res))
    outlier_flag[np.abs(new_res) >= (data_std * std_cut)] = True

    return outlier_flag


def zero_func(x, axis=None):
    """Dummy function to return a zero.
    Can be used as the centre function in astropy.stats.sigma_clip to get std relative to zero rather than median/mean value.

    Parameters
    -----------

    x:
        Dummy variable
    axis:
        required to match the syntax of numpy functions such as np.median
    """
    return 0


# def sigma_clip(data_res, kwargs={"sigma":3, "maxiters": 1, "cenfunc": zero_func}):
def sigma_clip(data_res, **kwargs):
    """Wrapper function for astropy.stats.sigma_clip, here we define the default centre of the data (the data - model residuals) to be zero

    Parameters
    -----------

    data_res: array
        The residuals of the data compared to the model.
    kwargs: dict
        Dictionary of keyword arguments from astropy.stats.sigma_clip

    Returns
    -----------
    sig_clip_mask: array
        returns only the mask from astropy.stats.sigma_clip
    """

    sig_clip = astropy_sigma_clip(data_res, **kwargs)
    sig_clip_mask = sig_clip.mask

    return sig_clip_mask


def outlier_sigma_diff(data_res, data_sigma, std_sigma=1):
    """Function to identify outliers by comparing the uncertainty of measurements to their residuals

    Parameters
    -----------

    data_res: array
        The residuals of the data compared to the model.
    data_sigma: array
        The uncertainties of the data points
    std_sigma: float
        Number of standard deviations to identify outliers, assuming the data uncertainties represent one standard deviation

    Returns
    -----------
    outlier_flag : array
       Array of flag indicating if data point is an outlier (True)
    """

    outlier_flag = np.abs(data_res) > (std_sigma * data_sigma)

    return outlier_flag


def apparition_gap_finder(x, dx=100.0):
    """Function to find gaps in a data series. E.g. given an array of observation times, find the different apparitions of an asteroid from gaps between observations larger than a given value.

    Parameters
    -----------

    x: array
        The SORTED data array to search for gaps
    dx: float
        The size of gap to identify in data series

    Returns
    -----------
    x_gaps: array
        Values of x which define the groups in the data, where each group is x_gaps[i] <= x < x_gaps[i+1]
    """

    ## TODO: check that x is sorted?

    # find the difference between each data point
    x_diff = np.diff(x)

    # select only differences greater than value dx
    x_diff_mask = x_diff > dx

    # find the values in x that correspond to the start of each group
    x_gaps = x[1:][x_diff_mask]

    # add the very first observation to the array
    x_gaps = np.insert(x_gaps, 0, x[0])

    # add the very last observation to the array, including a shift of dx to ensure the boundary inequalities work
    x_gaps = np.insert(x_gaps, len(x_gaps), x[-1] + dx)

    return x_gaps


def get_df_obs_filt(planetoid, filt, x_col="midPointMjdTai", x1=None, x2=None, col_list=None, pc_model=None):
    """Retrieve a dataframe of observations in a given filter. Has the option to limit the observations to a range of values, e.g. times/phase angles, if required.

    Parameters
    -----------

    planetoid: object
        Adler planetoid object containging the observations
    filt: str
        The filter to query
    x_col: str
        Column name to use for ordering values and limiting observations to a range: x1 <= df_obs[x_col] <= x2
    x1: float
        Lower limit value for x_col
    x2: float
        Upper limit value for x_col
    col_list: list
        List of column names to retrieve in addition to x_col, otherwise all columns are retrieved.
        N.B. if AbsMag is included in col_list then the absolute magnitude is calculated for the phase curve model pc_model
    pc_model: object
        Adler PhaseCurve model used to calculate AbsMag if required

    Returns
    -----------
    df_obs: DataFrame
        DataFrame of observations in the requested filter, ordered by x_col and with any x_col limits applied
    """

    # get observations in filter as a dataframe
    obs = planetoid.observations_in_filter(filt)
    # TODO: add option to get all filters? calculating AbsMag gets complicated...
    # TODO: split the dataframe functions into a separate function?
    df_obs = pd.DataFrame(obs.__dict__)

    # if no list of columns is provided, retrieve all columns
    if col_list is None:
        col_list = list(df_obs)
    else:
        col_list = [x_col] + col_list

    # calculate the phase curve absolute magnitudes (i.e. remove phase curve variation from reduced_mag)
    if "AbsMag" in col_list:
        # calculate the model absolute magnitude
        # TODO: add robustness to the units, phaseAngle and reduced_mag must match pc_model
        # For now we must assume that there are no units, and that degrees have been passed...
        # df_obs["AbsMag"] = pc_model.AbsMag(obs.phaseAngle * u.deg, obs.reduced_mag * u.mag).value
        df_obs["AbsMag"] = pc_model.AbsMag(np.radians(obs.phaseAngle), obs.reduced_mag)

    # select only the required columns
    df_obs = df_obs[col_list]

    # sort into x_col order
    df_obs = df_obs.sort_values(x_col)
    # limit to range in x_col if limits x1, x2 are supplied
    if (x1 is not None) & (x2 is not None):  # TODO: add functionality to set only upper/lower limit
        x_mask = (df_obs[x_col] >= x1) & (df_obs[x_col] <= x2)
        df_obs = df_obs[x_mask]
    # reset the dataframe index
    df_obs = df_obs.reset_index(drop=True)

    return df_obs
