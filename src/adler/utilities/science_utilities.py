import numpy as np
from astropy.stats import sigma_clip as astropy_sigma_clip


def outlier_diff(new_res, diff_cut=1.0):
    """Test whether new data point(s) is an outlier compared to the model by considering the residual difference.

    new_res: array
        The residuals of the new data points compared to the model
    diff_cut:
        The threshold difference value for outlier detection.

    Returns
    -----------
    outlier_flag : bool
        Flag indicating if data point is an outlier (True)

    """

    if not isinstance(new_res, np.ndarray):
        new_res = np.array(new_res)

    outlier_flag = np.array([False] * len(new_res))
    outlier_flag[np.abs(new_res) >= diff_cut] = True

    return outlier_flag


def outlier_std(new_res, data_res, std_cut=3.0):
    """Test whether new data point(s) is an outlier compared to the model by considering the standard deviation of the residuals.

    new_res: float or array
        The residuals of the new data point(s) compared to the model
    data_res: array
        The residuals of the data compared to the model.
    std_cut: float
       The threshold standard deviation for outlier detection.

    Returns
    -----------
    outlier_flag : bool
       Flag indicating if data point is an outlier (True)

    """

    # calculate the standard deviation of the data - model residuals
    data_std = np.std(data_res)

    if not isinstance(new_res, np.ndarray):
        new_res = np.array(new_res)

    outlier_flag = np.array([False] * len(new_res))
    outlier_flag[np.abs(new_res) >= (data_std * std_cut)] = True

    return outlier_flag


def zero_func(x, axis=None):
    """Dummy function to return a zero.
    Can be used as the centre function in astropy.stats.sigma_clip to get std relative to zero rather than median/mean value.

    axis:
        required to match the syntax of numpy functions such as np.median
    """
    return 0


def sigma_clip(data_res, kwargs={"maxiters": 1, "cenfunc": zero_func}):
    """Wrapper function for astropy.stats.sigma_clip, here we define the default centre of the data (the data - model residuals) to be zero

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
