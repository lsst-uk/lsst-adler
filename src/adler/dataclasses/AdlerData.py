from dataclasses import dataclass, field
import numpy as np


@dataclass
class AdlerData:
    """
    Class for storing Adler-calculated values.

    Attributes:
    -----------
    filter_list : list
        List of filters under investigation.

    model_lists : list
        List of lists of models per-filter. Length = len(filter_list) such that [[filter1_model1, filter1_model2], [filter2_model1]]. Used to index values
        in the H_adler, phase_parameter_1, phase_parameter_1_err, phase_parameter_2 and phase_parameter_2_err list structures.

    phaseAngle_min_adler : array_like
        Minimum phase angle of observations used in fitting model (degrees). Size = len(filter_list).

    phaseAngle_range_adler : array_like
        Max minus min phase angle range of observations used in fitting model (degrees). Size = len(filter_list).

    nobs_adler : array_like
        Number of observations used in fitting model. Size = len(filter_list).

    arc_adler : array_like
        Observational arc used to fit model (days). Size = len(filter_list).

    H_adler : list
        Absolute magnitude. List of lists arranged identically to model_lists, with per-filter and per-model values.

    H_err_adler : list
        Error in absolute magnitude. List of lists arranged identically to model_lists, with per-filter and per-model values.

    phase_parameters : list
        Phase parameters of the model in question. List of lists arranged as model_lists, with per-filter and per-model values: however, phase_parameters[filter_index][model_index] is itself a
        list of either one or two values, depending on the number of phase parameters.

    phase_parameter_err : list
        Error in the first phase parameter. List of lists arranged identically to phase_parameters.

    """

    filter_list: list
    model_lists: list = field(default_factory=list)

    phaseAngle_min_adler: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=float))
    phaseAngle_range_adler: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=float))
    nobs_adler: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    arc_adler: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=float))

    # Note that these are lists and not arrays because appending to Numpy arrays (i.e. to add parameters for a different model) is not recommended.
    # They can be cast to Numpy arrays later if needed.
    H_adler: list = field(default_factory=list, init=False)
    H_err_adler: list = field(default_factory=list, init=False)
    phase_parameters: list = field(default_factory=list, init=False)
    phase_parameters_err: list = field(default_factory=list, init=False)

    def __post_init__(self):
        """This runs post-initialisation and creates all of the class attributes where one dimension is "filters" to ensure the arrays
        and lists have the correct size. This makes population a little easier.
        """
        filter_length = len(self.filter_list)

        self.model_lists = [[] for a in range(0, filter_length)]

        self.phaseAngle_min_adler = np.zeros(filter_length, dtype=float)
        self.phaseAngle_range_adler = np.zeros(filter_length, dtype=float)
        self.nobs_adler = np.zeros(filter_length, dtype=int)
        self.arc_adler = np.zeros(filter_length, dtype=float)

        self.H_adler = [[] for a in range(0, filter_length)]
        self.H_err_adler = [[] for a in range(0, filter_length)]
        self.phase_parameters = [[] for a in range(0, filter_length)]
        self.phase_parameters_err = [[] for a in range(0, filter_length)]

    def populate_phase_parameters(
        self,
        filter_name,
        model_name,
        phaseAngle_min,
        phaseAngle_range,
        nobs,
        arc,
        H,
        H_err,
        parameters,
        parameters_err,
    ):
        """Convenience method to correctly populate phase curve arrays/lists.

        Parameters
        -----------

        filter_name : str
            The one-letter name of the filter in which the phase curve was calculated.

        model_name : str
            The name of the model used to calculate the phase curve.

        phaseAngle_min : float
            Minimum phase angle of observations used in fitting model (degrees)

        phaseAngle_range : float
            Max minus min phase angle range of observations used in fitting model (degrees).

        nobs : int
            Number of observations used in fitting model.

        arc : float
            Observational arc used to fit model (days).

        H : float
            Absolute magnitude in model.

        H_err : float
            Error on the absolute magnitude.

        parameters : list
            Phase parameters of the model.

        parameters_err : list
            Phase parameter errors.

        """

        # Make sure the supplied filter is in the filter list.
        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            raise Exception("Filter {} is not in supplied filter list.".format(filter_name))

        # If parameters and/or parameters_err are not lists, error out.
        if not isinstance(parameters, list) or not isinstance(parameters_err, list):
            raise TypeError("Both parameters and parameters_err arguments must be lists.")

        self.phaseAngle_min_adler[filter_index] = phaseAngle_min
        self.phaseAngle_range_adler[filter_index] = phaseAngle_range
        self.nobs_adler[filter_index] = nobs
        self.arc_adler[filter_index] = arc

        if model_name not in self.model_lists[filter_index]:
            self.model_lists[filter_index].append(model_name)

            self.H_adler[filter_index].append(H)
            self.H_err_adler[filter_index].append(H_err)
            self.phase_parameters[filter_index].append(parameters)
            self.phase_parameters_err[filter_index].append(parameters_err)

        else:
            model_index = self.model_lists[filter_index].index(model_name)

            self.H_adler[filter_index][model_index] = H
            self.H_err_adler[filter_index][model_index] = H_err
            self.phase_parameters[filter_index][model_index] = parameters
            self.phase_parameters_err[filter_index][model_index] = parameters_err
