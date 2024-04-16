from dataclasses import dataclass, field
import numpy as np


FILTER_DEPENDENT_KEYS = ["phaseAngle_min", "phaseAngle_range", "nobs", "arc"]
MODEL_DEPENDENT_KEYS = [
    "H",
    "H_err",
    "phase_parameter_1",
    "phase_parameter_1_err",
    "phase_parameter_2",
    "phase_parameter_2_err",
]


@dataclass
class AdlerData:
    """
    Class for storing Adler-calculated values.

    Attributes:
    -----------
    ssObjectId : str
        ssObjectId of the object of interest.

    filter_list : list of str
        List of filters under investigation.

    filter_dependent_values : list of FilterDependentAdler objects, optional
        List of FilterDependentAdler objects containing filter-dependent data in order of filter_list. Default empty list.

    """

    ssObjectId: str
    filter_list: list

    filter_dependent_values: list = field(default_factory=list)

    def __post_init__(self):
        """This runs post-initialisation and creates the class attribute where one dimension is "filters" to ensure the array
        has the correct size. This makes population a little easier.
        """

        # note that we don't do the same for model-dependent values as we don't know a priori how many models the user wishes
        # to calculate, but we do know how many filters the AdlerPlanetoid object was generated with
        self.filter_dependent_values = [FilterDependentAdler(filter_name) for filter_name in self.filter_list]

    def populate_phase_parameters(self, filter_name, **kwargs):
        """Convenience method to correctly populate phase curve parameters for a given filter and (if desired) model.
        Only the supplied arguments to the method will be updated, allowing for only some values to be populated if desired.

        Parameters
        -----------

        filter_name : str
            The one-letter name of the filter in which the phase curve was calculated.
        **kwargs : FilterDependentAdler and PhaseModelDependentAdler attributes
            The attribute names of the parameters you wish to update. See docs for FilterDependentAdler and PhaseModelDependentAdler
            classes for definitions of each attribute.
            Valid keyword arguments are: model_name, phaseAngle_min, phaseAngle_range, nobs, arc, H, H_err, phase_parameter_1,
            phase_parameter_1_err, phase_parameter_2, phase_parameter_2_err.
            Note that to update any of the model-dependent parameters (H, H_err, etc.), you WILL need to supply a model_name.

        """

        # make sure the supplied filter is in the filter list
        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        # if model-dependent parameters exist without a model name, return an error
        if not kwargs.get("model_name") and any(name in kwargs for name in MODEL_DEPENDENT_KEYS):
            raise NameError("No model name given. Cannot update model-specific phase parameters.")

        # update the value if it's in **kwargs
        for filter_key in FILTER_DEPENDENT_KEYS:
            if kwargs.get(filter_key):
                setattr(self.filter_dependent_values[filter_index], filter_key, kwargs.get(filter_key))

        # if no model_name is supplied, just end here
        # else, if the model does not exist for this filter, create it
        if not kwargs.get("model_name"):
            return
        elif kwargs.get("model_name") not in self.filter_dependent_values[filter_index].model_list:
            self.filter_dependent_values[filter_index].model_list.append(kwargs.get("model_name"))
            self.filter_dependent_values[filter_index].model_dependent_values.append(
                PhaseModelDependentAdler(filter_name, kwargs.get("model_name"))
            )

        # then get the model index
        model_index = self.filter_dependent_values[filter_index].model_list.index(kwargs.get("model_name"))

        # update the value if it's in **kwargs
        for model_key in MODEL_DEPENDENT_KEYS:
            if kwargs.get(model_key):
                setattr(
                    self.filter_dependent_values[filter_index].model_dependent_values[model_index],
                    model_key,
                    kwargs.get(model_key),
                )

    def print_data(self):
        """Convenience method to clearly print the stored values."""

        for f, filter_name in enumerate(self.filter_list):
            print("Filter: {}".format(filter_name))
            print("Phase angle minimum: {}".format(self.filter_dependent_values[f].phaseAngle_min))
            print("Phase angle range: {}".format(self.filter_dependent_values[f].phaseAngle_range))
            print("Number of observations: {}".format(self.filter_dependent_values[f].nobs))
            print("Arc: {}".format(self.filter_dependent_values[f].arc))

            for m, model_name in enumerate(self.filter_dependent_values[f].model_list):
                print("Model: {}.".format(model_name))
                print("\tH: {}".format(self.filter_dependent_values[f].model_dependent_values[m].H))
                print("\tH error: {}".format(self.filter_dependent_values[f].model_dependent_values[m].H_err))
                print(
                    "\tPhase parameter 1: {}".format(
                        self.filter_dependent_values[f].model_dependent_values[m].phase_parameter_1
                    )
                )
                print(
                    "\tPhase parameter 1 error: {}".format(
                        self.filter_dependent_values[f].model_dependent_values[m].phase_parameter_1_err
                    )
                )
                print(
                    "\tPhase parameter 2: {}".format(
                        self.filter_dependent_values[f].model_dependent_values[m].phase_parameter_2
                    )
                )
                print(
                    "\tPhase parameter 2 error: {}".format(
                        self.filter_dependent_values[f].model_dependent_values[m].phase_parameter_2_err
                    )
                )

            print("\n")

    def get_phase_parameters_in_filter(self, filter_name, model_name=None):
        """Convenience method to return the phase parameters in a specific filter and model.

        Parameters
        -----------
        filter_name : str
            The filter of interest.

        model_name : str, optional
            The model name of the model of interest. If this is not supplied, the code will not return any model-dependent
            parameters. Default None.


        Returns
        -----------
        output_obj : PhaseParameterOutput object
            Object containing phase curve parameters for the specified filter and model.

        """

        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        output_obj = PhaseParameterOutput()
        output_obj.filter_name = filter_name
        output_obj.phaseAngle_min = self.filter_dependent_values[filter_index].phaseAngle_min
        output_obj.phaseAngle_range = self.filter_dependent_values[filter_index].phaseAngle_range
        output_obj.nobs = self.filter_dependent_values[filter_index].nobs
        output_obj.arc = self.filter_dependent_values[filter_index].arc

        if not model_name:
            print("No model name specified. Returning non-model-dependent phase parameters.")
        else:
            try:
                model_index = self.filter_dependent_values[filter_index].model_list.index(model_name)
            except ValueError:
                raise ValueError(
                    "Model {} does not exist for filter {} in AdlerData.model_lists.".format(
                        model_name, filter_name
                    )
                )

            output_obj.model_name = model_name
            output_obj.H = self.filter_dependent_values[filter_index].model_dependent_values[model_index].H
            output_obj.H_err = (
                self.filter_dependent_values[filter_index].model_dependent_values[model_index].H_err
            )
            output_obj.phase_parameter_1 = (
                self.filter_dependent_values[filter_index]
                .model_dependent_values[model_index]
                .phase_parameter_1
            )
            output_obj.phase_parameter_1_err = (
                self.filter_dependent_values[filter_index]
                .model_dependent_values[model_index]
                .phase_parameter_1_err
            )
            output_obj.phase_parameter_2 = (
                self.filter_dependent_values[filter_index]
                .model_dependent_values[model_index]
                .phase_parameter_2
            )
            output_obj.phase_parameter_2_err = (
                self.filter_dependent_values[filter_index]
                .model_dependent_values[model_index]
                .phase_parameter_2_err
            )

        return output_obj


@dataclass
class FilterDependentAdler:
    """Dataclass containing filter-dependent values generated by Adler. Note that NaN indicates a value that has not yet been populated.

    Attributes:
    -----------
    filter_name : str
        The filter for which these values are calculated.

    phaseAngle_min : float, optional
        Minimum phase angle of observations used in fitting model (degrees).

    phaseAngle_range : float, optional
        Max minus min phase angle range of observations used in fitting model (degrees).

    nobs : int, optional
        Number of observations used in fitting model.

    arc: float, optional
        Observational arc used to fit model (days).

    model_list: list of str, optional
        List of the models for which phase curve parameters have been calculated. Default: empty list

    model_dependent_values: list of PhaseModelDependentAdler objects, optional
        List of PhaseModelDependentAdler objects storing phase-model parameters for each model, given in order of model_list. Default: empty list.

    """

    filter_name: str
    phaseAngle_min: float = np.nan
    phaseAngle_range: float = np.nan
    nobs: int = 0
    arc: float = np.nan
    model_list: list = field(default_factory=list)
    model_dependent_values: list = field(default_factory=list)


@dataclass
class PhaseModelDependentAdler:
    """Dataclass containing phase-model-dependent values generated by Adler. Note that NaN indicates a value that has not yet been populated.

    Attributes:
    -----------
    filter_name : str
        The filter for which these values are calculated.

    model_name : str
        The phase model for which these values were calculated. Example: "HG", "HG1G2", "linear".

    H : float, optional
        The absolute magnitude. Default NaN.

    H_err : float, optional
        Error in absolute magnitude. Default NaN.

    phase_parameter_1 : float, optional
        The first parameter of the phase model. May be the only parameter. For example, G in the HG model. Default NaN.

    phase_parameter_1_err : float, optional
        The error on the first parameter of the phase model. Default NaN.

    phase_parameter_2 : float, optional
        The second parameter of the phase model. May not exist for this model. Default NaN.

    phase_parameter_2_err : float, optional
        The error on the second parameter of the phase model. Default NaN.

    """

    filter_name: str
    model_name: str
    H: float = np.nan
    H_err: float = np.nan
    phase_parameter_1: float = np.nan
    phase_parameter_1_err: float = np.nan
    phase_parameter_2: float = np.nan
    phase_parameter_2_err: float = np.nan


class PhaseParameterOutput:
    """Empty convenience class so that the output of AdlerData.get_phase_parameters_in_filter is an object."""

    pass
