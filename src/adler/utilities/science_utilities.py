import numpy as np
from astropy.stats import sigma_clip as astropy_sigma_clip


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


def sigma_clip(data_res, kwargs={"maxiters": 1, "cenfunc": zero_func}):
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
