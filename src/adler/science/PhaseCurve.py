from sbpy.photometry import HG, HG1G2, HG12, HG12_Pen16, LinearPhaseFunc
import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter


class PhaseCurve:
    """A class to define the phasecurve model and associated functions.
    Units - by default no units are set but astropy units can be passed.
    It is up to the user to ensure that units are correct for the relevant phasecurve model.

    Attributes
    -----------
    abs_mag : float
       Absolute magnitude, H, of the phasecurve model (often units of mag).

    phase_param: float
       The first phase parameter of the phasecurve model, e.g. G from HG
       (often dimensionless units unless S from LinearPhaseFunc, which has units mag/deg or mag/rad).

    phase_param2: float
       The second phase parameter, only used for the 3 parameter HG1G2 phasecurve model.

    model_name: str
       Label for the phasecurve model to be used.
       Choice of: "HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc"

    bounds: bool
       Flag for using the default bounds on the sbpy model (True) or ignoring bounds (False).

    """

    def __init__(self, abs_mag=18, phase_param=0.2, phase_param2=None, model_name="HG"):
        self.abs_mag = abs_mag
        self.phase_param = phase_param
        self.phase_param2 = phase_param2
        self.model_name = model_name

        if model_name == "HG":
            self.model_function = HG(H=abs_mag, G=self.phase_param)
        elif model_name == "HG1G2":
            self.model_function = HG1G2(H=abs_mag, G1=self.phase_param, G2=self.phase_param)
        elif model_name == "HG12":
            self.model_function = HG12(H=abs_mag, G12=self.phase_param)
        elif model_name == "HG12_Pen16":
            self.model_function = HG12_Pen16(H=abs_mag, G12=self.phase_param)
        elif model_name == "LinearPhaseFunc":
            self.model_function = LinearPhaseFunc(H=abs_mag, S=self.phase_param)
        else:
            print("no model selected")

    def SetModelBounds(self, param, bound_vals=(None, None)):
        model_sbpy = self.model_function
        param_names = model_sbpy.param_names
        x = getattr(model_sbpy, param)
        setattr(x, "bounds", bound_vals)

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
           Dictionary containing the PhaseCurve parameters you wish to set, e.g. abs_mag, phase_param

        Returns
        ----------

        model : object
           The new PhaseCurve class object

        """

        model = PhaseCurve()
        for key, value in model_dict.items():
            setattr(model, key, value)
        return model

    def InitModelSbpy(self, model_sbpy):
        """Set up a new PhaseCurve model object from an existing sbpy model
        ### or create dict from sbpy model and then use InitModelDict?

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
        param_names = model_sbpy.param_names
        parameters = []
        for p in param_names:
            # try get the quantity (value with units)
            x = getattr(model_sbpy, p).quantity
            # if there are no units get just the value
            if x is None:
                x = getattr(model_sbpy, p).value
            parameters.append(x)
        # print(param_names, parameters)

        # create a PhaseCurve object with the extracted parameters
        model = PhaseCurve(*parameters, model_name=model_name)

        return model

    def ReducedMag(self, phase_angle):
        """Return the reduced magnitude of the phasecurve model for a given phase angle(s)

        phase_angle - value or array, must have astropy units of degrees

        Parameters
        -----------
        phase_angle : float or array
           value or array of phase angles at which to evaluate the phasecurve model, must have astropy units of degrees.

        Returns
        ----------

        return_value : float or array
           The phasecurve model reduced magnitude at the given phase angle(s)

        """

        return self.model_function(phase_angle)

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
           Select a fitting function from astropy.modeling.fitting, defaults to astropy.modeling.fitting.LevMarLSQFitter

        Returns
        ----------

        model_fit : object
           The sbpy phasecurve model object

        """

        # use the LevMarLSQFitter by default
        if fitter is None:
            fitter = LevMarLSQFitter()
        # print(fitter)

        if mag_err is not None:  # fit weighted by photometric uncertainty
            model_fit = fitter(self.model_function, phase_angle, reduced_mag, weights=1.0 / mag_err)
        else:  # unweighted fit
            model_fit = fitter(self.model_function, phase_angle, reduced_mag)

        ### if overwrite_model: # add an overwrite option?
        # redo __init__ with the new fitted parameters

        return model_fit
