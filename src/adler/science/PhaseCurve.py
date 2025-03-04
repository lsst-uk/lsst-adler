from sbpy.photometry import HG, HG1G2, HG12, HG12_Pen16, LinearPhaseFunc
import astropy.units as u
import numpy as np
from astropy.modeling.fitting import LevMarLSQFitter
import itertools

# translation between sbpy and adler field names
SBPY_ADLER_DICT = {
    "H": "H",
    "G": "phase_parameter_1",
    "G12": "phase_parameter_1",
    "G1": "phase_parameter_1",
    "G2": "phase_parameter_2",
    "S": "phase_parameter_1",
}
# translation between adler and sbpy fields, depends on model
ADLER_SBPY_DICT = {
    "HG": {"phase_parameter_1": "G"},
    "HG1G2": {"phase_parameter_1": "G1", "phase_parameter_2": "G2"},
    "HG12": {"phase_parameter_1": "G12"},
    "HG12_Pen16": {"phase_parameter_1": "G12"},
    "LinearPhaseFunc": {"phase_parameter_1": "S"},
}


class PhaseCurve:
    """A class to define the phasecurve model and associated functions.
    Units - by default no units are set but astropy units can be passed.
    It is up to the user to ensure that units are correct for the relevant phasecurve model.

    Attributes
    -----------
    H : float
       Absolute magnitude, H, of the phasecurve model (often units of mag).

    phase_parameter_1: float
       The first phase parameter of the phasecurve model, e.g. G from HG
       (often dimensionless units unless S from LinearPhaseFunc, which has units mag/deg or mag/rad).

    phase_parameter_2: float
       The second phase parameter, only used for the 3 parameter HG1G2 phasecurve model.

    model_name: str
       Label for the phasecurve model to be used.
       Choice of: "HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc"

    bounds: bool
       Flag for using the default bounds on the sbpy model (True) or ignoring bounds (False).

    H_err : float
       Uncertainty in absolute magnitude.

    phase_parameter_1_err : float
       Uncertainty in the (first) phase parameter.

    phase_parameter_2_err : float
       Uncertainty in the (second, if required) phase parameter.

    """

    def __init__(
        self,
        H=18,
        phase_parameter_1=0.2,
        phase_parameter_2=0.2,
        model_name="HG",
        H_err=None,
        phase_parameter_1_err=None,
        phase_parameter_2_err=None,
    ):
        self.H = H
        self.phase_parameter_1 = phase_parameter_1
        self.phase_parameter_2 = phase_parameter_2
        self.model_name = model_name
        self.H_err = H_err
        self.phase_parameter_1_err = phase_parameter_1_err
        self.phase_parameter_2_err = phase_parameter_2_err

        if model_name == "HG":
            self.model_function = HG(H=H, G=self.phase_parameter_1)
        elif model_name == "HG1G2":
            self.model_function = HG1G2(H=H, G1=self.phase_parameter_1, G2=self.phase_parameter_2)
        elif model_name == "HG12":
            self.model_function = HG12(H=H, G12=self.phase_parameter_1)
        elif model_name == "HG12_Pen16":
            self.model_function = HG12_Pen16(H=H, G12=self.phase_parameter_1)
        elif model_name == "LinearPhaseFunc":
            self.model_function = LinearPhaseFunc(H=H, S=self.phase_parameter_1)
        else:
            print("no model selected")

    def SetModelBounds(self, param, bound_vals=(None, None)):
        """Set the "bounds" attribute of an sbpy model parameter, i.e. the lower and upper constraints for the fitter.

        Parameters
        -----------
        param : str
           Parameter name bounds to be set fix, if not an sbpy parameter name (e.g. G) the corresponding adler name is looked up (e.g. phase_parameter_1)

        bound_vals : tuple
           Set the fitter constraints to (upper, lower) for the param.

        """
        model_sbpy = self.model_function
        if param not in model_sbpy.__dict__:
            param = ADLER_SBPY_DICT[self.model_name][param]
        x = getattr(model_sbpy, param)
        setattr(x, "bounds", bound_vals)
        return

    def FixParam(self, param, fix_flag=True):
        """Set the "fixed" attribute of an sbpy model parameter.
        E.g. use this to fit for absolute magnitude whilst keeping phase parameter fixed.

        Parameters
        -----------
        param : str
           Parameter name to fix, if not an sbpy parameter name (e.g. G) the corresponding adler name is looked up (e.g. phase_parameter_1)

        fix_flag : bool
           Set True to keep the param fixed when fitting

        """

        model_sbpy = self.model_function
        if param not in model_sbpy.__dict__:
            param = ADLER_SBPY_DICT[self.model_name][param]
        x = getattr(model_sbpy, param)
        setattr(x, "fixed", fix_flag)
        return

    def ReturnModelDict(self):
        """Return the values for the PhaseCurve class as a dict

        Returns
        ----------

        self.__dict__ : dict
           The dict of PhaseCurve object parameters.

        """

        return self.__dict__

    def InitModelDict(self, model_dict):
        """Set up a new PhaseCurve model object from a dictionary.
        This could be written by the user or generated from another PhaseCurve object using ReturnModelDict

        Parameters
        -----------
        model_dict : dict
           Dictionary containing the PhaseCurve parameters you wish to set, e.g. H, phase_parameter_1

        Returns
        ----------

        model : object
           The new PhaseCurve class object

        """

        # clean the input dictionary
        del_keys = []
        if "model_function" in model_dict.keys():
            del_keys.append("model_function")
        for key, value in model_dict.items():
            if not hasattr(self, key):
                del_keys.append(key)
        model_dict = model_dict.copy()  # make a copy to avoid changing the original dict
        for key in del_keys:
            model_dict.pop(key, None)

        # initialise a new model
        model = PhaseCurve(**model_dict)

        return model

    def InitModelSbpy(self, model_sbpy):
        """Set up a new PhaseCurve model object from an existing sbpy model
        ### TODO or create dict from sbpy model and then use InitModelDict?

        Parameters
        -----------
        model_sbpy : object
           The sbpy model object, e.g. HG()

        Returns
        ----------

        model : object
           The new PhaseCurve class object

        """
        # get model name from the sbpy model object
        model_name = model_sbpy.__class__.name

        # get the sbpy model parameters
        param_names = list(model_sbpy.param_names)
        # we will also check for any uncertainties we have stored in the sbpy object
        param_names_err = ["{}_err".format(x) for x in param_names]
        param_names = param_names + param_names_err

        # create a dictionary of phase curve parameters from sbpy in a format accepted by PhaseCurve
        parameters = {}
        for p in param_names:
            if p in model_sbpy.__dict__:  # check that the parameter is available in the sbpy object
                x = getattr(model_sbpy, p)
                # try get the quantity (value with units)
                if hasattr(x, "unit"):
                    if (x.unit is None) or (
                        x.unit == ""
                    ):  # if there are no units (or weird blank units?) get just the value
                        x = x.value
                    else:
                        x = x.quantity
                # look up the correct adler parameter name (accounting for additional uncertainty, "_err", parameters)
                if p.endswith("_err"):  # assumes the uncertainty parameter always ends in "_err"
                    _p = SBPY_ADLER_DICT[p.split("_err")[0]] + "_err"
                    parameters[_p] = x
                else:
                    parameters[SBPY_ADLER_DICT[p]] = x

        # create a PhaseCurve object with the extracted parameters
        model = PhaseCurve(**parameters, model_name=model_name)

        return model

    def ReducedMag(self, phase_angle):
        """Return the reduced magnitude of the phasecurve model for a given phase angle(s)

        phase_angle - value or array, must have astropy units of degrees

        Parameters
        -----------
        phase_angle : float or array
           value or array of phase angles at which to evaluate the phasecurve model.
           Must have astropy units of degrees if sbpy model uses units, otherwise pass radian values.

        Returns
        ----------

        return_value : float or array
           The phasecurve model reduced magnitude at the given phase angle(s)

        """

        return self.model_function(phase_angle)

    def ReducedMagBounds(self, phase_angle, std=1):
        """Accounting for the parameter uncertainties, return the minimum and maximum reduced magnitudes of the phasecurve model for a given phase angle(s). This is done by determining the minimum and maximum values of each parameter given the uncertainty (and some std, assuming a gaussian distribution). A PhaseCurve model is made for all unique combinations of the min, max of each parameter. For each PhaseCurve model the ReducedMag is calculated; the overall min/max ReducedMag is determined at each phase angle.

        Parameters
        -----------
        phase_angle : float or array
           value or array of phase angles at which to evaluate the phasecurve model.
           Must have astropy units of degrees if sbpy model uses units, otherwise pass radian values.
        std : float
           Number of standard deviations we want to consider when determining the minimum/maximum parameter values, where we assume that the parameter uncertainty is representative of a 1 sigma std in a gaussian distribution.
        Returns
        ----------

        pc_bounds : dict
           mag_min - The minimum reduced magnitude of all PhaseCurve models at every phase angle
           mag_max - The maximum reduced magnitude of all PhaseCurve models at every phase angle
           mag_mean - The mean reduced magnitude of all PhaseCurve models at every phase angle
           PhaseCurves - A list of all unique PhaseCurve models given the possible combinations of parameters and their uncertainties

        """

        # Get all parameters in the model that have an uncertainty
        err_keys = [x for x in self.__dict__.keys() if "_err" in x and getattr(self, x) is not None]

        # calculate the min, max for each parameter with uncertainty
        params = {}
        for x_err in err_keys:
            x = x_err.replace("_err", "")
            _x = getattr(self, x)
            _x_err = getattr(self, x_err)
            plus = _x + (std * _x_err)
            minus = _x - (std * _x_err)
            params[x] = [plus, minus]

        # Get all possible combinations of the min, max for each parameter
        p_vals = [params[x] for x in params.keys()]
        all_p_vals = list(itertools.product(*p_vals))

        # make a new PhaseCurve model for each set of parameters
        # Determine the reduced magnitudes for the given phase angle range
        mags = []
        pcs = []
        for p in all_p_vals:
            _dict = self.__dict__.copy()
            for i, x in enumerate(params.keys()):
                _dict[x] = p[i]
            _pc = PhaseCurve().InitModelDict(_dict)
            pcs.append(_pc)
            _mag = _pc.ReducedMag(phase_angle)
            mags.append(_mag)
        mags = np.array(mags)

        # Calculate the magnitude statistics
        # Get the absolute maximum and minimum value at each phase angle for all possible models
        mag_min = np.amin(mags, axis=0)
        mag_max = np.amax(mags, axis=0)
        mag_mean = np.mean(mags, axis=0)

        # Store results in a dictionary
        pc_bounds = {"mag_min": mag_min, "mag_max": mag_max, "mag_mean": mag_mean, "PhaseCurves": pcs}

        return pc_bounds

    def ModelResiduals(self, phase_angle, reduced_mag):
        """For a set of phase curve observations, return the residuals to the PhaseCurve model.
        NB that units must match the sbpy model. E.g. phase_angle should be passed with units of degrees, or be in values of radians

        Parameters
        -----------
        phase_angle : float or array
           value or array of phase angles at which to evaluate the phasecurve model.

        reduced_mag : float or array
           value or array of reduced magnitudes at which to evaluate the phasecurve model.

        Returns
        ----------

        residuals : float or array
           The residuals of the observations minus PhaseCurve model values

        """

        residuals = reduced_mag - self.ReducedMag(phase_angle)

        return residuals

    def AbsMag(self, phase_angle, reduced_mag):
        """For a set of phase curve observations, return the absolute magnitude from the fitted phase curve model.
        I.e. this is the model residuals, shifted by the fitted absolute magnitude
        NB that units for phase_angle and reduced_mag must match the sbpy model. E.g. phase_angle should be passed with units of degrees, or be in values of radians

        Parameters
        -----------
        phase_angle : float or array
           value or array of phase angles at which to evaluate the phasecurve model.

        reduced_mag : float or array
           value or array of reduced magnitudes at which to evaluate the phasecurve model.

        Returns
        ----------

        abs_mag : float or array
           The residuals of the observations minus PhaseCurve model values shifted by the model absolute magnitude

        """

        # TODO: add option to pass model & filt instead of reduced_mag & phase angle - will calculate all absmags automatically
        # probably not possible as observations are contained in a separate object
        # if (phase_angle is None) and (reduced_mag is None) and (filt is not None)...

        abs_mag = reduced_mag - self.ReducedMag(phase_angle) + self.H

        return abs_mag

    def FitModel(self, phase_angle, reduced_mag, mag_err=None, fitter=None):
        """Fit the phasecurve model parameters to observations.
        starts with a phase curve model as an initial guess for parameters.
        fits model to phase angle and reduced magnitude.

        phase_angle - phase angle of each observations
        reduced_mag - distance corrected reduced magnitudes
        mag_err - photometric uncertainties to weight the measurements
        fitter - can pass a fitting function from astropy.modeling.fitting, defaults to astropy.modeling.fitting.LevMarLSQFitter

        Parameters
        -----------
        phase_angle : float or array
           The Sun-object-observer phase angles of the observations.

        reduced_mag : float or array
           The observed reduced magnitudes at the corresponding phase angles.

        mag_err : float or array
           Uncertainty on the reduced magnitude, used to weight the fit.

        fitter : object
           Select a fitting function from astropy.modeling.fitting, defaults to astropy.modeling.fitting.LevMarLSQFitter.
           N.B. that LevMarLSQFitter cannot handle inequality constraints for the HG1G2 model, use something like SLSQPLSQFitter from astropy.modeling.fitting (does not return covariance matrix!).

        Returns
        ----------

        model_fit : object
           The sbpy phasecurve model object

        """

        # use the LevMarLSQFitter by default
        if fitter is None:
            fitter = LevMarLSQFitter()

        if mag_err is not None:  # fit weighted by photometric uncertainty
            model_fit = fitter(self.model_function, phase_angle, reduced_mag, weights=1.0 / mag_err)
        else:  # unweighted fit
            model_fit = fitter(self.model_function, phase_angle, reduced_mag)

        # Add fitted uncertainties as an additional attribute within the sbpy object
        if "param_cov" in fitter.fit_info:
            # get the covariance matrix from the fit
            covariance = fitter.fit_info["param_cov"]
            if covariance is not None:
                # get fit uncertainties as square of the diagonal of the covariance
                fit_errs = np.sqrt(np.diag(covariance))
                # update only the uncertainties for parameters used in the fit
                param_names = np.array(model_fit.param_names)
                fit_mask = ~np.array([getattr(model_fit, x).fixed for x in param_names])
                for i, x in enumerate(param_names[fit_mask]):
                    setattr(model_fit, "{}_err".format(x), fit_errs[i])
            # else:
            ### TODO log covariance is None error here

        ### TODO
        # else:
        #     log lack of uncertainties for fitter
        #     run an MCMC estimate of uncertainty?

        ### if overwrite_model: # add an overwrite option?
        # redo __init__ with the new fitted parameters
        # this would then return an adler PhaseCurve object rather than an sbpy object

        return model_fit
