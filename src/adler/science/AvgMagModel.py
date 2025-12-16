import numpy as np


class AvgMagModel:
    # TODO docstring

    def __init__(self, avg_mag=np.nan, std_mag=np.nan, model_name="median"):
        """# TODO docstring"""
        self.avg_mag = avg_mag
        self.std_mag = std_mag
        self.model_name = model_name

    # TODO what functions?
    # TODO how does np.median (nanmedian?) and np.mean (nanmean) etc come into this
    # TODO what do we want from this class

    def InitModelObs(self, mag, magErr=None, model_name="median"):
        """# TODO docstring"""
        if model_name == "median":
            model_dict = {"avg_mag":np.nanmedian(mag), "std_mag":np.nanstd(mag), "model_name":model_name}
        elif model_name == "mean":
            model_dict = {"avg_mag":np.nanmean(mag), "std_mag":np.nanstd(mag), "model_name":model_name}
        else:
            print("Invalid model selected. Options are median or mean")
        
        model = AvgMagModel(**model_dict)

        return model
