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

    phaseAngle_min_adler: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    phaseAngle_range_adler: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    nobs_adler: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    arc_adler: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))

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

        self.phaseAngle_min_adler = np.empty(filter_length, dtype=float)
        self.phaseAngle_range_adler = np.empty(filter_length, dtype=float)
        self.nobs_adler = np.empty(filter_length, dtype=int)
        self.arc_adler = np.empty(filter_length, dtype=float)

        self.H_adler = [[] for a in range(0, filter_length)]
        self.H_err_adler = [[] for a in range(0, filter_length)]
        self.phase_parameters = [[] for a in range(0, filter_length)]
        self.phase_parameters_err = [[] for a in range(0, filter_length)]

    def populate_phase_parameters(
        self,
        filter_name,
        model_name=None,
        phaseAngle_min=None,
        phaseAngle_range=None,
        nobs=None,
        arc=None,
        H=None,
        H_err=None,
        parameters=[],
        parameters_err=[],
    ):
        """Convenience method to correctly populate phase curve arrays/lists. Only the supplied arguments to the method will
        be updated, allowing for only some values to be populated if desired.

        Parameters
        -----------

        filter_name : str
            The one-letter name of the filter in which the phase curve was calculated.

        model_name : str, optional
            The name of the model used to calculate the phase curve. If not supplied, no model-dependent parameters will be updated.

        phaseAngle_min : float, optional
            Minimum phase angle of observations used in fitting model (degrees)

        phaseAngle_range : float, optional
            Max minus min phase angle range of observations used in fitting model (degrees).

        nobs : int, optional
            Number of observations used in fitting model.

        arc : float, optional
            Observational arc used to fit model (days).

        H : float, optional
            Absolute magnitude in model.

        H_err : float, optional
            Error on the absolute magnitude.

        parameters : list, optional
            Phase parameters of the model.

        parameters_err : list, optional
            Phase parameter errors.

        """

        # Make sure the supplied filter is in the filter list.
        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        # If parameters and/or parameters_err are not lists, error out.
        if not isinstance(parameters, list) or not isinstance(parameters_err, list):
            raise TypeError("Both parameters and parameters_err arguments must be lists.")

        self.phaseAngle_min_adler[filter_index] = self.update_value(
            self.phaseAngle_min_adler[filter_index], phaseAngle_min
        )
        self.phaseAngle_range_adler[filter_index] = self.update_value(
            self.phaseAngle_range_adler[filter_index], phaseAngle_range
        )
        self.nobs_adler[filter_index] = self.update_value(self.nobs_adler[filter_index], nobs)
        self.arc_adler[filter_index] = self.update_value(self.arc_adler[filter_index], arc)

        if model_name is None and (
            any(v is not None for v in [H, H_err]) or len(parameters) != 0 or len(parameters_err) != 0
        ):
            raise Exception("No model name given. Cannot update model-specific parameters.")

        elif model_name is None:
            pass

        elif model_name not in self.model_lists[filter_index]:
            self.model_lists[filter_index].append(model_name)

            self.H_adler[filter_index].append(H)
            self.H_err_adler[filter_index].append(H_err)
            self.phase_parameters[filter_index].append(parameters)
            self.phase_parameters_err[filter_index].append(parameters_err)

        else:
            model_index = self.model_lists[filter_index].index(model_name)

            self.H_adler[filter_index][model_index] = self.update_value(
                self.H_adler[filter_index][model_index], H
            )
            self.H_err_adler[filter_index][model_index] = self.update_value(
                self.H_err_adler[filter_index][model_index], H_err
            )
            self.phase_parameters[filter_index][model_index] = self.update_value(
                self.phase_parameters[filter_index][model_index], parameters
            )
            self.phase_parameters_err[filter_index][model_index] = self.update_value(
                self.phase_parameters_err[filter_index][model_index], parameters_err
            )

    def update_value(self, original_value, new_value):
        """Returns one of two values depending on whether one of them is None. Used to update class attributes.

        Parameters
        -----------
        original_value : any
            Original value of the attribute.

        new_value : any or None
            The value to replace it with, if this is not None.


        Returns
        -----------
        any
            Either original_value (if new_value is None) or new_value.

        """
        if not new_value:
            return original_value
        else:
            return new_value

    def print_data(self):
        """Convenience method to clearly print the stored values."""

        print("Phase parameters (per filter):\n")
        for f, filter_name in enumerate(self.filter_list):
            print("Filter: {}".format(filter_name))
            print("Phase angle minimum: {}".format(self.phaseAngle_min_adler[f]))
            print("Phase angle range: {}".format(self.phaseAngle_range_adler[f]))
            print("Number of observations: {}".format(self.nobs_adler[f]))
            print("Arc: {}".format(self.arc_adler[f]))

            for m, model_name in enumerate(self.model_lists[f]):
                print("Model: {}.".format(model_name))
                print("\tH: {}".format(self.H_adler[f][m]))
                print("\tH error: {}".format(self.H_err_adler[f][m]))
                print("\tPhase parameter(s): {}".format(self.phase_parameters[f][m]))
                print("\tPhase parameter(s) error: {}".format(self.phase_parameters_err[f][m]))

            print("\n")

    def get_phase_parameters_in_filter(self, filter_name, model_name=None):
        """Convenience method to return the phase parameters in a specific filter and model.

        Parameters
        -----------
        filter_name : str
            The filter of interest.

        model_name : str, optional
            The model name of the model of interest. If this is not supplied, the code automatically returns the
            parameters for the first model in the list.


        Returns
        -----------
        parameters_dict : dict
            Dictionary of the phase curve parameters for the specified filter and model.

        """

        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        if model_name:
            try:
                model_index = self.model_lists[filter_index].index(model_name)
            except ValueError:
                raise ValueError(
                    "Model {} does not exist for filter {} in AdlerData.model_lists.".format(
                        model_name, filter_name
                    )
                )
        else:
            model_index = 0
            print(
                "No model name specified. Returning phase parameters for first model in list: {}.".format(
                    self.model_lists[filter_index][model_index]
                )
            )

        parameters_dict = {
            "phaseAngle_min": self.phaseAngle_min_adler[filter_index],
            "phaseAngle_range": self.phaseAngle_range_adler[filter_index],
            "nobs": self.nobs_adler[filter_index],
            "arc": self.arc_adler[filter_index],
            "H": self.H_adler[filter_index][model_index],
            "H_err": self.H_err_adler[filter_index][model_index],
            "phase_parameters": self.phase_parameters[filter_index][model_index],
            "phase_parameters_err": self.phase_parameters_err[filter_index][model_index],
        }

        return parameters_dict
