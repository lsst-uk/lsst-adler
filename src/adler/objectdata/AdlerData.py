import os
import sqlite3
import logging
import re
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Optional
from astropy.time import Time

FILTER_DEPENDENT_KEYS = [
    "phaseAngle_min",
    "phaseAngle_range",
    "observationTime_max",
    "nobs",
    "arc",
    "n_outliers",
    "n_std_outliers",
    "sustained_outliers",
]
PHASE_MODEL_DEPENDENT_KEYS = [
    "H",
    "H_err",
    "phase_parameter_1",
    "phase_parameter_1_err",
    "phase_parameter_2",
    "phase_parameter_2_err",
]
VALID_PHASE_MODELS = ["HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc"]
AVG_MAG_MODEL_DEPENDENT_KEYS = [
    "avg_mag",
    "std_mag",
]
VALID_AVG_MAG_MODELS = ["median", "mean"]
ALL_FILTER_LIST = ["u", "g", "r", "i", "z", "y"]

VALID_MODELS = sorted(
    VALID_PHASE_MODELS + VALID_AVG_MAG_MODELS, key=len, reverse=True
)  # sorted to avoid partial matches when using _get_model_name

logger = logging.getLogger(__name__)

# Ensure that numpy dtypes correctly map to SQL types
sqlite3.register_adapter(np.float64, float)
sqlite3.register_adapter(np.float32, float)
sqlite3.register_adapter(np.int64, int)
sqlite3.register_adapter(np.int32, int)


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

    modelId : str, optional
        modelId for the model of interest that has been or is to be computed. Default: Empty str

    updatedMJD : float, optional
        Timestamp (in MJD) of the time that this AdlerData object was initialized.

    filter_dependent_values : list of FilterDependentAdler objects, optional
        List of FilterDependentAdler objects containing filter-dependent data in order of filter_list. Default empty list.

    """

    ssObjectId: str
    filter_list: list
    modelId: str = ""
    updatedMJD: float = np.nan
    filter_dependent_values: list = field(default_factory=list)

    def __post_init__(self):
        """This runs post-initialisation and creates the class attribute where one dimension is "filters" to ensure the array
        has the correct size. This makes population a little easier. We also generate the current MJD timestamp to record the time this AdlerData object was initialized.
        """

        # note that we don't do the same for model-dependent values as we don't know a priori how many models the user wishes
        # to calculate, but we do know how many filters the AdlerPlanetoid object was generated with
        self.filter_dependent_values = [FilterDependentAdler(filter_name) for filter_name in self.filter_list]

        self.updatedMJD = Time.now().mjd

    def _MJD_update(self):
        """
        Function for updating the updatedMJD value stored in AdlerData. This should be called whenever an update is made to AdlerData so that it is clear what is the most up-to-date version of AdlerData.
        """

        self.updatedMJD = Time.now().mjd

    def set_modelId(self, model_name, end_mjd, data_timespan, n_new_nights):
        """
        Function for setting the modelId parameter.

        Parameters
        -----------
        model_name : str
            The model name for the given model calculated. One of "HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc", "median", "mean".

        end_mjd : float
            The MJD set as the maximum MJD to consider in the model.

        data_timespan : float
            The number of nights of data that is considered for the given model.

        n_new_nights : float
            The number of nights of data that are considered as "new observations" in calculating outliers.
        """
        # N.B. double underscore is intentional to provide a point to break this string if searching for the model_name after generation (as HG12_Pen16 is a valid model name with a single underscore)
        self.modelId = f"{model_name}_{end_mjd:.1f}_{data_timespan}n_{n_new_nights}n"

        self._MJD_update()

    def populate_filter_dependent_parameters(self, filter_name, **kwargs):
        """Convenience method to correctly populate the filter-dependent parameters for a given filter.
        Only the supplied arguments to the method will be updated, allowing for only some values to be populated if desired.

        filter_name : str
            The one-letter name of the filter of interest.
        **kwargs : FilterDependentAdler attributes
            The attribute names of the parameters you wish to update. See docs for FilterDependentAdler class for definitions of each attribute.
            Valid keyword arguments are: model_name, phaseAngle_min, phaseAngle_range, observationTime_max, arc, nobs, n_outliers, n_std_outliers, sustained_outliers.

        """

        # make sure the supplied filter is in the filter list
        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        # update the value if it's in **kwargs
        for filter_key in FILTER_DEPENDENT_KEYS:
            if kwargs.get(filter_key):
                setattr(self.filter_dependent_values[filter_index], filter_key, kwargs.get(filter_key))

        self._MJD_update()

    def populate_phase_parameters(self, filter_name, **kwargs):
        """Convenience method to correctly populate phase curve parameters for a given filter and model.
        Only the supplied arguments to the method will be updated, allowing for only some values to be populated if desired.
        This method will automatically populate filter-dependent parameters also.

        Parameters
        -----------

        filter_name : str
            The one-letter name of the filter in which the phase curve was calculated.
        **kwargs : FilterDependentAdler and PhaseModelDependentAdler attributes
            The attribute names of the parameters you wish to update. See docs for FilterDependentAdler and PhaseModelDependentAdler
            classes for definitions of each attribute.
            Valid keyword arguments are: model_name, phaseAngle_min, phaseAngle_range, observationTime_max, arc, nobs, n_outliers, n_std_outliers, sustained_outliers, H, H_err, phase_parameter_1, phase_parameter_1_err, phase_parameter_2, phase_parameter_2_err.
            Note that to update any of the model-dependent parameters (H, H_err, etc.), you WILL need to supply a model_name.

        """

        # make sure the supplied filter is in the filter list
        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        # populate the filter dependent parameters
        self.populate_filter_dependent_parameters(filter_name, **kwargs)

        # if model-dependent parameters exist without a model name, return an error
        if not kwargs.get("model_name") and any(name in kwargs for name in PHASE_MODEL_DEPENDENT_KEYS):
            logger.error("NameError: No model name given. Cannot update model-specific phase parameters.")
            raise NameError("No model name given. Cannot update model-specific phase parameters.")

        # if no model_name is supplied, just end here
        # else, if the model does not exist for this filter, create it
        if not kwargs.get("model_name"):
            return
        elif kwargs.get("model_name") != self.filter_dependent_values[filter_index].model_name:
            logger.warning(
                f"Input model name {kwargs.get("model_name")} does not match model name {self.filter_dependent_values[filter_index].model_name} in AdlerData. Parameters will be overwritten."
            )
            self.filter_dependent_values[filter_index].model_name = kwargs.get("model_name")
            self.filter_dependent_values[filter_index].model_dependent_values = PhaseModelDependentAdler(
                filter_name, kwargs.get("model_name")
            )

        # update the value if it's in **kwargs
        for model_key in PHASE_MODEL_DEPENDENT_KEYS:
            if model_key in kwargs:
                setattr(
                    self.filter_dependent_values[filter_index].model_dependent_values,
                    model_key,
                    kwargs.get(model_key),
                )

        self._MJD_update()

    def populate_avg_mag_parameters(self, filter_name, **kwargs):
        """Convenience method to correctly populate average magnitude model parameters for a given filter and model.
        Only the supplied arguments to the method will be updated, allowing for only some values to be populated if desired.
        This method will automatically populate filter-dependent parameters also.

        Parameters
        -----------

        filter_name : str
            The one-letter name of the filter in which the model was calculated.
        **kwargs : FilterDependentAdler and AvgMagModelDependentAdler attributes
            The attribute names of the parameters you wish to update. See docs for FilterDependentAdler and AvgMagModelDependentAdler
            classes for definitions of each attribute.
            Valid keyword arguments are: model_name, phaseAngle_min, phaseAngle_range, observationTime_max, arc, nobs, n_outliers, n_std_outliers, sustained_outliers, avg_mag, std_mag.
            Note that to update any of the model-dependent parameters (avg_mag, std_mag), you WILL need to supply a model_name.

        """

        # make sure the supplied filter is in the filter list
        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        # populate the filter dependent parameters
        self.populate_filter_dependent_parameters(filter_name, **kwargs)

        # if model-dependent parameters exist without a model name, return an error
        if not kwargs.get("model_name") and any(name in kwargs for name in AVG_MAG_MODEL_DEPENDENT_KEYS):
            logger.error(
                "NameError: No model name given. Cannot update model-specific average magnitude parameters."
            )
            raise NameError("No model name given. Cannot update model-specific average magnitude parameters.")

        # if no model_name is supplied, just end here
        # else, if the model does not exist for this filter, create it
        if not kwargs.get("model_name"):
            return
        elif kwargs.get("model_name") != self.filter_dependent_values[filter_index].model_name:
            logger.warning(
                f"Input model name {kwargs.get("model_name")} does not match model name {self.filter_dependent_values[filter_index].model_name} in AdlerData. Parameters will be overwritten."
            )
            self.filter_dependent_values[filter_index].model_name = kwargs.get("model_name")
            self.filter_dependent_values[filter_index].model_dependent_values = AvgMagModelDependentAdler(
                filter_name, kwargs.get("model_name")
            )

        # update the value if it's in **kwargs
        for model_key in AVG_MAG_MODEL_DEPENDENT_KEYS:
            if model_key in kwargs:
                setattr(
                    self.filter_dependent_values[filter_index].model_dependent_values,
                    model_key,
                    kwargs.get(model_key),
                )

        self._MJD_update()

    def populate_source_flags(self, filter_name, modelId, df, **kwargs):
        """Convenience method to correctly populate the source outlier flags for a given filter and modelId.
        Only the supplied arguments to the method will be updated, allowing for only some values to be populated if desired.
        Observations detected as outliers (source flags) must be supplied in a pandas.DataFrame.
        This method will automatically populate filter-dependent parameters also.

        Parameters
        -----------

        filter_name : str
            The one-letter name of the filter in which the model was calculated.

        modelId : str
            modelId that the given outliers correspond to. This is used to check the supplied modelId from the user matches that stored in AdlerData.

        df : pandas.DataFrame
            DataFrame of the observations that are identified as outliers. Must contain columns diaSourceId, midPointMjdTai, mag_diff, std_diff.

        **kwargs : FilterDependentAdler and AvgMagModelDependentAdler attributes
            The attribute names of the parameters you wish to update. See docs for FilterDependentAdler
            classes for definitions of each attribute.
            Valid keyword arguments are: model_name, phaseAngle_min, phaseAngle_range, observationTime_max, arc, nobs, n_outliers, n_std_outliers, sustained_outliers.

        """

        if modelId != self.modelId:
            logger.error(
                f"ValueError: modelId {modelId} does not match the modelId in AdlerData.modelId: {self.modelId}"
            )
            raise ValueError(
                f"modelId {modelId} does not match the modelId in AdlerData.modelId: {self.modelId}"
            )

        # make sure the supplied filter is in the filter list
        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        # populate the filter dependent parameters
        # Add n_outliers, n_std_outliers to kwargs so it populates FilterDependentAdler
        kwargs.update({"n_outliers": len(df.loc[df.mag_diff != 0])})
        kwargs.update({"n_std_outliers": len(df.loc[df.std_diff != 0])})

        self.populate_filter_dependent_parameters(filter_name, **kwargs)

        self.filter_dependent_values[filter_index].source_flags = AdlerSourceFlags.construct_from_data_table(
            self.ssObjectId, filter_name, modelId, df
        )

        self._MJD_update()

    def populate_from_database(self, filepath, modelId=None):
        """Populates the AdlerData object with information from the most recent timestamped entry for the ssObjectId in a given database.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database. Note that for now, we assume only one table with all the data.

        modelId : str, optional
            modelId for the model of interest that should be recovered. Default: None.
        """

        if modelId:
            self.modelId = modelId

        con = self._get_database_connection(filepath)
        cursor = con.cursor()

        tbl_list = self._get_tables(con)

        # Ensure AdlerData is always the first table queried (thus populating modelId if not already specified)
        tbl_list.insert(0, tbl_list.pop(tbl_list.index("AdlerData")))

        for tbl_name in tbl_list:
            logger.info(f"Populating information from {tbl_name}.")
            # Specific query required for AdlerSourceFlags whereas the other tables can follow the same style (as they are unique on ssObjectId)
            if tbl_name == "AdlerSourceFlags":
                source_flags_query = f"SELECT * FROM AdlerSourceFlags WHERE ssObjectId='{self.ssObjectId}' and modelId='{self.modelId}'"
                cursor.execute(source_flags_query)
                rows = cursor.fetchall()
                columns = [desc[0] for desc in cursor.description]
                df = pd.DataFrame(rows, columns=columns)

                for filter_name in self.filter_list:
                    _df = df.loc[df.filter_name == filter_name]
                    self.populate_source_flags(filter_name, self.modelId, _df)

            else:
                if not modelId and tbl_name == "AdlerData":
                    # If modelId isn't specified and we're querying AdlerData (i.e. the default first table to be queried), we take the most recent entry
                    sql_query = f"""SELECT * from {tbl_name} WHERE ssObjectId='{self.ssObjectId}' ORDER BY updatedMJD DESC LIMIT 1"""
                else:
                    sql_query = f"""SELECT * from {tbl_name} WHERE ssObjectId='{self.ssObjectId}' AND modelId='{self.modelId}' ORDER BY updatedMJD DESC LIMIT 1"""
                query_result = cursor.execute(sql_query)

                try:
                    fetched_data_raw = query_result.fetchall()[0]
                except IndexError:
                    logger.error(
                        f"ValueError: No data found in table {tbl_name} in this database for the supplied ssObjectId/modelId."
                    )
                    raise ValueError(
                        f"ValueError: No data found in table {tbl_name} in this database for the supplied ssObjectId/modelId."
                    )

                fetched_data = [
                    np.nan if v is None else v for v in fetched_data_raw
                ]  # replaces Nones with nans
                column_list = self._get_database_columns(con, tbl_name)

                row_data = dict(zip(column_list, fetched_data))

                # Set modelId
                if not modelId:
                    self.modelId = row_data["modelId"]

                # TODO decide whether to update updatedMJD (and make it updatedMJD)
                self.updatedMJD = row_data["updatedMJD"]

                filter_bools = [
                    any((column_heading.startswith(filter + "_") for column_heading in column_list))
                    for filter in ALL_FILTER_LIST
                ]
                database_filter_list = [b for a, b in zip(filter_bools, ALL_FILTER_LIST) if a]

                if not all(
                    [requested_filter in database_filter_list for requested_filter in self.filter_list]
                ):
                    logger.error(
                        "ValueError: Data does not exist for some of the requested filters in this database. Filters in database for this object: {}".format(
                            database_filter_list
                        )
                    )
                    raise ValueError(
                        "Data does not exist for some of the requested filters in this database. Filters in database for this object: {}".format(
                            database_filter_list
                        )
                    )

                for filter_name in self.filter_list:
                    expected_filter_columns = [
                        filter_name + "_" + filter_key for filter_key in FILTER_DEPENDENT_KEYS
                    ]
                    filter_columns = [col for col in column_list if col in expected_filter_columns]
                    filter_values = [row_data[col] for col in filter_columns]
                    present_filter_columns = [col.strip(filter_name + "_") for col in filter_columns]
                    filter_dependent_info = dict(zip(present_filter_columns, filter_values))

                    self.populate_filter_dependent_parameters(filter_name, **filter_dependent_info)

                    model_name = self._get_model_name()

                    if model_name in VALID_PHASE_MODELS:
                        expected_model_columns = [
                            filter_name + "_" + model_name + "_" + model_key
                            for model_key in PHASE_MODEL_DEPENDENT_KEYS
                        ]
                        model_columns = [col for col in column_list if col in expected_model_columns]
                        model_values = [row_data[col] for col in model_columns]
                        model_dependent_info = dict(zip(PHASE_MODEL_DEPENDENT_KEYS, model_values))
                        model_dependent_info["model_name"] = model_name

                        self.populate_phase_parameters(filter_name, **model_dependent_info)
                    elif model_name in VALID_AVG_MAG_MODELS:
                        expected_model_columns = [
                            filter_name + "_" + model_name + "_" + model_key
                            for model_key in AVG_MAG_MODEL_DEPENDENT_KEYS
                        ]
                        model_columns = [col for col in column_list if col in expected_model_columns]
                        model_values = [row_data[col] for col in model_columns]
                        model_dependent_info = dict(zip(AVG_MAG_MODEL_DEPENDENT_KEYS, model_values))
                        model_dependent_info["model_name"] = model_name

                        self.populate_avg_mag_parameters(filter_name, **model_dependent_info)
                    else:
                        logger.error(
                            f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                        )
                        raise ValueError(
                            f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                        )

        con.close()

    def print_data(self):
        """Convenience method to clearly print the stored values."""

        for f, filter_name in enumerate(self.filter_list):
            print("Filter: {}".format(filter_name))
            print("Phase angle minimum: {}".format(self.filter_dependent_values[f].phaseAngle_min))
            print("Phase angle range: {}".format(self.filter_dependent_values[f].phaseAngle_range))
            print("Maximum observation time: {}".format(self.filter_dependent_values[f].observationTime_max))
            print("Number of observations: {}".format(self.filter_dependent_values[f].nobs))
            print("Arc: {}".format(self.filter_dependent_values[f].arc))
            print("Number of outliers detected: {}".format(self.filter_dependent_values[f].n_outliers))
            print(
                "Number of outliers in sigma-space detected: {}".format(
                    self.filter_dependent_values[f].n_std_outliers
                )
            )
            print(
                "Magnitude change of sustained outliers: {}".format(
                    self.filter_dependent_values[f].sustained_outliers
                )
            )

            model_name = self.filter_dependent_values[f].model_name

            print("Model: {}.".format(model_name))
            if model_name in VALID_PHASE_MODELS:
                print("H: {}".format(self.filter_dependent_values[f].model_dependent_values.H))
                print("H error: {}".format(self.filter_dependent_values[f].model_dependent_values.H_err))
                print(
                    "Phase parameter 1: {}".format(
                        self.filter_dependent_values[f].model_dependent_values.phase_parameter_1
                    )
                )
                print(
                    "Phase parameter 1 error: {}".format(
                        self.filter_dependent_values[f].model_dependent_values.phase_parameter_1_err
                    )
                )
                print(
                    "Phase parameter 2: {}".format(
                        self.filter_dependent_values[f].model_dependent_values.phase_parameter_2
                    )
                )
                print(
                    "Phase parameter 2 error: {}".format(
                        self.filter_dependent_values[f].model_dependent_values.phase_parameter_2_err
                    )
                )
            elif model_name in VALID_AVG_MAG_MODELS:
                print(
                    "Average magnitude {}".format(
                        self.filter_dependent_values[f].model_dependent_values.avg_mag
                    )
                )
                print(
                    "Standard deviation of magnitudes {}".format(
                        self.filter_dependent_values[f].model_dependent_values.std_mag
                    )
                )
            else:
                logger.error(
                    f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                )
                raise ValueError(
                    f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
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
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        output_obj = PhaseParameterOutput()
        output_obj.filter_name = filter_name
        output_obj.phaseAngle_min = self.filter_dependent_values[filter_index].phaseAngle_min
        output_obj.phaseAngle_range = self.filter_dependent_values[filter_index].phaseAngle_range
        output_obj.observationTime_max = self.filter_dependent_values[filter_index].observationTime_max
        output_obj.nobs = self.filter_dependent_values[filter_index].nobs
        output_obj.arc = self.filter_dependent_values[filter_index].arc
        output_obj.n_outliers = self.filter_dependent_values[filter_index].n_outliers
        output_obj.n_std_outliers = self.filter_dependent_values[filter_index].n_std_outliers
        output_obj.sustained_outliers = self.filter_dependent_values[filter_index].sustained_outliers

        if not model_name:
            logger.warning("No model name was specified. Returning non-model-dependent phase parameters.")
            print("No model name specified. Returning non-model-dependent phase parameters.")
        elif model_name == self.filter_dependent_values[filter_index].model_name:
            output_obj.model_name = model_name
            output_obj.H = self.filter_dependent_values[filter_index].model_dependent_values.H
            output_obj.H_err = self.filter_dependent_values[filter_index].model_dependent_values.H_err
            output_obj.phase_parameter_1 = self.filter_dependent_values[
                filter_index
            ].model_dependent_values.phase_parameter_1
            output_obj.phase_parameter_1_err = self.filter_dependent_values[
                filter_index
            ].model_dependent_values.phase_parameter_1_err
            output_obj.phase_parameter_2 = self.filter_dependent_values[
                filter_index
            ].model_dependent_values.phase_parameter_2
            output_obj.phase_parameter_2_err = self.filter_dependent_values[
                filter_index
            ].model_dependent_values.phase_parameter_2_err
        else:
            logger.error(
                f"Input model name {model_name} does not match model name {self.filter_dependent_values[filter_index].model_name} in AdlerData. Cannot get phase parameters for {model_name}."
            )
            raise ValueError(
                f"Input model name {model_name} does not match model name {self.filter_dependent_values[filter_index].model_name} in AdlerData. Cannot get phase parameters for {model_name}."
            )

        return output_obj

    def get_avg_mag_parameters_in_filter(self, filter_name, model_name=None):
        """Convenience method to return the average magnitude parameters in a specific filter and model.

        Parameters
        -----------
        filter_name : str
            The filter of interest.

        model_name : str, optional
            The model name of the model of interest. If this is not supplied, the code will not return any model-dependent
            parameters. Default None.


        Returns
        -----------
        output_obj : AvgMagParameterOutput object
            Object containing average magnitude model parameters for the specified filter and model.

        """

        try:
            filter_index = self.filter_list.index(filter_name)
        except ValueError:
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        output_obj = AvgMagParameterOutput()
        output_obj.filter_name = filter_name
        output_obj.phaseAngle_min = self.filter_dependent_values[filter_index].phaseAngle_min
        output_obj.phaseAngle_range = self.filter_dependent_values[filter_index].phaseAngle_range
        output_obj.observationTime_max = self.filter_dependent_values[filter_index].observationTime_max
        output_obj.nobs = self.filter_dependent_values[filter_index].nobs
        output_obj.arc = self.filter_dependent_values[filter_index].arc
        output_obj.n_outliers = self.filter_dependent_values[filter_index].n_outliers
        output_obj.n_std_outliers = self.filter_dependent_values[filter_index].n_std_outliers
        output_obj.sustained_outliers = self.filter_dependent_values[filter_index].sustained_outliers

        if not model_name:
            logger.warning("No model name was specified. Returning non-model-dependent phase parameters.")
            print("No model name specified. Returning non-model-dependent phase parameters.")
        elif model_name == self.filter_dependent_values[filter_index].model_name:
            output_obj.model_name = model_name
            output_obj.avg_mag = self.filter_dependent_values[filter_index].model_dependent_values.avg_mag
            output_obj.std_mag = self.filter_dependent_values[filter_index].model_dependent_values.std_mag
        else:
            logger.error(
                f"Input model name {model_name} does not match model name {self.filter_dependent_values[filter_index].model_name} in AdlerData. Cannot get average magnitude parameters for {model_name}."
            )
            raise ValueError(
                f"Input model name {model_name} does not match model name {self.filter_dependent_values[filter_index].model_name} in AdlerData. Cannot get average magnitude parameters for {model_name}."
            )

        return output_obj

    def _get_database_connection(self, filepath, table_name=None):
        """Returns the connection to the output SQL database, creating it and the given table if it does not exist.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.

        table_name : str, optional
            Name of the table to create if it doesn't exist. (This replaces the create_new argument)

        Returns
        ----------
        con : sqlite3 Connection object
            The connection to the output database.

        """

        database_exists = os.path.isfile(
            filepath
        )  # check this FIRST as the next statement creates the db if it doesn't exist

        if not database_exists and table_name:  # we need to make the table and a couple of starter columns
            con = sqlite3.connect(filepath)
            cur = con.cursor()
            cur.execute(f"CREATE TABLE {table_name}(ssObjectId PRIMARY KEY, modelId, updatedMJD REAL)")
        elif not database_exists and not table_name:
            logger.error("ValueError: Database cannot be found at given filepath.")
            raise ValueError("Database cannot be found at given filepath.")
        elif database_exists and not table_name:
            # If no table_name specified connect and return the connection
            con = sqlite3.connect(filepath)
        else:
            # If the database exists, connect to it
            con = sqlite3.connect(filepath)
            cur = con.cursor()
            # If the table doesn't exist (i.e. we're creating the table for writing, this will create it)
            cur.execute(
                f"CREATE TABLE IF NOT EXISTS {table_name}(ssObjectId PRIMARY KEY, modelId, updatedMJD REAL)"
            )

        return con

    def _get_row_data_and_columns(self, table_name):
        """Collects all of the data present in the AdlerData object as a list with a corresponding list of column names,
        in preparation for a row to be written to a SQL database in the given table_name.

        Returns
        -----------
        table_name : str
            Name of the table that we want to collect data for in preparation for writing to the database.

        """
        required_columns = ["ssObjectId", "modelId", "updatedMJD"]
        row_data = [self.ssObjectId, self.modelId, self.updatedMJD]

        for f, filter_name in enumerate(self.filter_list):
            if table_name == "AdlerData":
                columns_by_filter = [
                    "_".join([filter_name, filter_key])
                    for filter_key in FILTER_DEPENDENT_KEYS
                    if "outliers" in filter_key
                ]
                data_by_filter = [
                    getattr(self.filter_dependent_values[f], filter_key)
                    for filter_key in FILTER_DEPENDENT_KEYS
                    if "outliers" in filter_key
                ]

                required_columns.extend(columns_by_filter)
                row_data.extend(data_by_filter)

            elif table_name == "FilterDependentAdler":
                columns_by_filter = [
                    "_".join([filter_name, filter_key]) for filter_key in FILTER_DEPENDENT_KEYS
                ]
                data_by_filter = [
                    getattr(self.filter_dependent_values[f], filter_key)
                    for filter_key in FILTER_DEPENDENT_KEYS
                ]

                required_columns.extend(columns_by_filter)
                row_data.extend(data_by_filter)

            elif table_name == "PhaseModelDependentAdler":
                model_name = self.filter_dependent_values[f].model_name
                if model_name == "":
                    logger.warning(
                        f"No models calculated for filter {filter_name}, continuing to next filter"
                    )
                    continue
                else:
                    columns_by_model = [
                        "_".join([filter_name, model_name, model_key])
                        for model_key in PHASE_MODEL_DEPENDENT_KEYS
                    ]
                    data_by_model = [
                        getattr(self.filter_dependent_values[f].model_dependent_values, model_key)
                        for model_key in PHASE_MODEL_DEPENDENT_KEYS
                    ]

                    required_columns.extend(columns_by_model)
                    row_data.extend(data_by_model)

            elif table_name == "AvgMagModelDependentAdler":
                model_name = self.filter_dependent_values[f].model_name
                if model_name == "":
                    logger.warning(
                        f"No models calculated for filter {filter_name}, continuing to next filter"
                    )
                    continue
                else:
                    columns_by_model = [
                        "_".join([filter_name, model_name, model_key])
                        for model_key in AVG_MAG_MODEL_DEPENDENT_KEYS
                    ]
                    data_by_model = [
                        getattr(self.filter_dependent_values[f].model_dependent_values, model_key)
                        for model_key in AVG_MAG_MODEL_DEPENDENT_KEYS
                    ]

                    required_columns.extend(columns_by_model)
                    row_data.extend(data_by_model)

            else:
                logger.error(
                    f"{table_name} must be one of [AdlerData, FilterDependentAlder, PhaseModelDependentAdler, AvgMagModelDependentAdler]"
                )
                raise ValueError(
                    (
                        f"{table_name} must be one of [AdlerData, FilterDependentAlder, PhaseModelDependentAdler, AvgMagModelDependentAdler]"
                    )
                )

        return row_data, required_columns

    def _get_database_columns(self, con, table_name):
        """Gets a list of the current columns in a given table in a SQL database.

        Parameters
        -----------
        con : sqlite3 Connection object
            The connection to the output SQL database.

        table_name : str
            The name of the relevant table in the database.


        Returns
        ----------
        list of str
            List of current columns existing in the table.

        """

        cur = con.cursor()
        cur.execute(f"""SELECT * from {table_name} where 1=0""")
        return [d[0] for d in cur.description]

    def _get_tables(self, con):
        """Gets a list of the current tables in a SQL database.

        Parameters
        -----------
        con : sqlite3 Connection object
            The connection to the output SQL database.


        Returns
        ----------
        list of str
            List of current tables existing in the database.

        """

        cur = con.cursor()
        cur.execute(f"SELECT tbl_name FROM sqlite_schema WHERE type='table'")
        res = cur.fetchall()
        return [r[0] for r in res]

    def _ensure_columns(self, con, table_name, current_columns, required_columns):
        """Creates new columns in a given table of a SQL database as needed by checking the list of current columns against a list
        of required columns.


        Parameters
        -----------
        con : sqlite3 Connection object
            The connection to the output SQL database.

        table_name : str
            The name of the relevant table in the database.

        current_columns : list of str
            A list of the columns already existing in the database table.

        required_columns : list of str
            A list of the columns needed in the database table.

        """

        cur = con.cursor()
        for column_name in required_columns:
            if column_name not in current_columns:
                cur.execute(f"""ALTER TABLE {table_name} ADD COLUMN {column_name}""")

    def _write_table(self, filepath, table_name):
        """
        Function for writing information to a given table in the given database.
        Connects to the database and creates table if necessary, gathers data from AdlerData, ensures required columns are present and writes to the database.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.

        table_name : str
            Name of table to write to. Must be one of AdlerData, FilterDependentAdler, PhaseModelDependentAdler, AvgMagModelDependentAdler

        """

        con = self._get_database_connection(filepath, table_name=table_name)

        row_data, required_columns = self._get_row_data_and_columns(table_name=table_name)
        current_columns = self._get_database_columns(con, table_name)
        self._ensure_columns(con, table_name, current_columns, required_columns)

        column_names = ",".join(required_columns)
        column_spaces = ",".join(["?"] * len(required_columns))
        update_clause = ", ".join([f"{col} = excluded.{col}" for col in required_columns[1:]])
        sql_command = f"""
                        INSERT INTO {table_name} ({column_names})
                        VALUES ({column_spaces})
                        ON CONFLICT(ssObjectId) DO UPDATE SET {update_clause};
                        """
        cur = con.cursor()
        cur.execute(sql_command, row_data)
        con.commit()
        con.close()

    def _get_model_name(self):
        """Returns the model_name by parsing the modelId

        Returns
        -----------
            Name of the model specified in self.modelId
        """
        matches = [m for m in VALID_MODELS if self.modelId.startswith(m + "_")]
        if len(matches) == 1:
            return matches[0]
        if not matches:
            logger.error(f"No valid model found in: {self.modelId}")
            raise ValueError(f"No valid model found in: {self.modelId}")
        logger.error(f"Ambiguous model match in: {self.modelId}")
        raise ValueError(f"Ambiguous model match in: {self.modelId}")

    def write_to_database(self, filepath, write_model_data=False):
        """Writes all of the relevant data contained within the AdlerData object to a SQLite database.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.

        write_model_data : Boolean, optional
            A flag to set whether to write out specific model data to AdlerData. Default: False.

        """

        # Write default AdlerData information
        self._write_table(filepath=filepath, table_name="AdlerData")
        logger.info(f"Top-level information written to AdlerData table")

        if write_model_data:
            # Write FilterDependentAdler data
            self._write_table(filepath=filepath, table_name="FilterDependentAdler")
            logger.info(f"Filter-specific information written to FilterDependentAdler table")

            model_name = self._get_model_name()
            if model_name in VALID_PHASE_MODELS:
                # Write PhaseModelDependentAdler data
                self._write_table(filepath=filepath, table_name="PhaseModelDependentAdler")
                logger.info(
                    f"Phase Model-specific information for model {model_name} written to PhaseModelDependentAdler table"
                )
            elif model_name in VALID_AVG_MAG_MODELS:
                # Write AvgMagModelDependentAdler data
                self._write_table(filepath=filepath, table_name="AvgMagModelDependentAdler")
                logger.info(
                    f"Average Magnitude model-specific information for model {model_name} written to AvgMagModelDependentAdler table"
                )
            else:
                logger.error(
                    f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                )
                raise ValueError(
                    f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                )

            # Write AdlerSourceFlags data
            for f, filter_name in enumerate(self.filter_list):
                filter_source_flags = self.filter_dependent_values[f].source_flags
                if filter_source_flags:
                    filter_source_flags.write_flags_to_database(filepath=filepath)
                    logger.info(
                        f"Source flags information written to AdlerSourceFlags for filter '{filter_name}'"
                    )
                else:
                    logger.info(f"No source flags for filter '{filter_name}', continuing to next filter")


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


@dataclass
class AvgMagModelDependentAdler:
    """Dataclass containing model-dependent values for the simple average magnitude model pgenerated by Adler. Note that NaN indicates a value that has not yet been populated.

    Attributes:
    -----------
    filter_name : str
        The filter for which these values are calculated.

    model_name : str
        The model for which these values were calculated. Example: "median", "mean".

    avg_mag : float, optional
        Average magnitude of the measurements used to calculate the model. Default NaN.

    std_mag : float, optional
        Standard deviation of the measurements used to calculate the model. Default NaN.

    """

    filter_name: str
    model_name: str
    avg_mag: float = np.nan
    std_mag: float = np.nan


class PhaseParameterOutput:
    """Empty convenience class so that the output of AdlerData.get_phase_parameters_in_filter is an object."""

    pass


class AvgMagParameterOutput:
    """Empty convenience class so that the output of AdlerData.get_avg_mag_parameters_in_filter is an object."""

    pass


@dataclass
class AdlerSourceFlags:
    """
    Class for storing Adler-determined outlier information.

    Attributes:
    -----------
    ssObjectId : str
        ssObjectId of the object of interest.

    filter_name : str
        Filter the observation was taken in.

    modelId : str
        modelId for the model that the outliers are compared to.

    n_outliers : int
        Number of observations identified as outliers.

    n_std_outliers : int, optional
        Number of outliers detected for the given model in sigma space.

    diaSourceId : array_like of ints or strs
        Unique identifier of the observation.

    midPointMjdTai : array_like of floats
        Observation timestamps.

    mag_diff : array_like of floats
        Differences in (reduced) magnitude between the observations and the model.

    std_diff : array_like of floats
        Deviation (in terms of the observations uncertainties) between the observations and the model.

    """

    ssObjectId: str
    filter_name: str
    modelId: str
    n_outliers: int
    n_std_outliers: int
    diaSourceId: np.ndarray = field(default_factory=lambda: np.zeros(0))
    midPointMjdTai: np.ndarray = field(default_factory=lambda: np.zeros(0))
    mag_diff: np.ndarray = field(default_factory=lambda: np.zeros(0))
    std_diff: np.ndarray = field(default_factory=lambda: np.zeros(0))

    @classmethod
    def construct_from_data_table(cls, ssObjectId, filter_name, modelId, df):
        """Method for constructing the AdlerSourceFlags object from a dataframe.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_name : str
        Filter the observation was taken in.

        modelId : str
            modelId for the model that the outliers are compared to.

        df : pandas.DataFrame
            DataFrame of the observations that are identified as outliers. Must contain columns diaSourceId, midPointMjdTai, mag_diff, std_diff

        Returns
        -----------
        AdlerSourceFlags object
            Object containing the source flags information on outliers identified.

        """
        obs_dict = {"ssObjectId": ssObjectId, "filter_name": filter_name, "modelId": modelId}

        obs_dict.update(
            df.loc[:, ["diaSourceId", "midPointMjdTai", "mag_diff", "std_diff"]].to_dict(orient="list")
        )

        obs_dict.update({"n_outliers": len(df.loc[df.mag_diff != 0])})
        obs_dict.update({"n_std_outliers": len(df.loc[df.std_diff != 0])})

        return cls(**obs_dict)

    def _get_database_connection(self, filepath, create_new=False):
        """Returns the connection to the output SQL database, creating it and the AdlerSource Flags table if it does not exist.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.

        create_new : Boolean
            Whether to create the database if it doesn't already exist. Default is False.

        Returns
        ----------
        con : sqlite3 Connection object
            The connection to the output database.

        """

        database_exists = os.path.isfile(
            filepath
        )  # check this FIRST as the next statement creates the db if it doesn't exist

        if not database_exists and create_new:  # we need to make the table and a couple of starter columns
            con = sqlite3.connect(filepath)
            cur = con.cursor()
            cur.execute(
                "CREATE TABLE AdlerSourceFlags(ssObjectId, filter_name, modelId, diaSourceId, midPointMjdTai, mag_diff, std_diff)"
            )
        elif not database_exists and not create_new:
            logger.error("ValueError: Database cannot be found at given filepath.")
            raise ValueError("Database cannot be found at given filepath.")
        else:
            con = sqlite3.connect(filepath)
            cur = con.cursor()
            # Create the table if it doesn't exist (in case database was created previously without this table)
            cur.execute(
                "CREATE TABLE IF NOT EXISTS AdlerSourceFlags(ssObjectId, filter_name, modelId, diaSourceId, midPointMjdTai, mag_diff, std_diff)"
            )

        return con

    def write_flags_to_database(self, filepath, table_name="AdlerSourceFlags"):
        """
        Writes the information from AdlerSourceFlags to the given database.

        Parameters
        -----------
        filepath : path-like object
            Path to the output database.

        table_name : str, optional
            Name of the table to write the flags to. Default: AdlerSourceFlags.
        """

        con = self._get_database_connection(filepath, create_new=True)

        required_columns = [
            "ssObjectId",
            "filter_name",
            "modelId",
            "diaSourceId",
            "midPointMjdTai",
            "mag_diff",
            "std_diff",
        ]

        column_names = ",".join(required_columns)
        column_spaces = ",".join(["?"] * len(required_columns))
        sql_command = f"""
                        INSERT INTO {table_name} ({column_names})
                        VALUES ({column_spaces});
                        """

        row_data = list(
            zip(
                [self.ssObjectId] * len(self.diaSourceId),
                [self.filter_name] * len(self.diaSourceId),
                [self.modelId] * len(self.diaSourceId),
                self.diaSourceId,
                self.midPointMjdTai,
                self.mag_diff,
                self.std_diff,
            )
        )

        cur = con.cursor()
        cur.executemany(sql_command, row_data)
        con.commit()
        con.close()


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

    observationTime_max : float, optional
        Maximum time of observation used in fitting the model (in modified Julian day).

    arc: float, optional
        Observational arc used to fit model (days).

    nobs : int, optional
        Number of observations used in fitting model.

    n_outliers : int, optional
        Number of outliers detected for the given model.

    n_std_outliers : int, optional
        Number of outliers detected for the given model in sigma space.

    sustained_outliers : float, optional
        Magnitude difference between old and new observations.

    model_name : str, optional
        Name of the model computed. Must be one of "HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc", "median", "mean".

    model_dependent_values : PhaseModelDependentAdler or AvgMagModelDependentAdler object, optional
        PhaseModelDependentAdler or AvgMagModelDependentAdler object storing phase-model or average-magnitude-model parameters for the given model. Default: None.

    source_flags : AdlerSourceFlags, optional
        AdlerSourceFlags object storing the information on specific observations identified as outliers compared to the given model.

    """

    filter_name: str
    phaseAngle_min: float = np.nan
    phaseAngle_range: float = np.nan
    observationTime_max: float = np.nan
    arc: float = np.nan
    nobs: int = 0
    n_outliers: int = 0
    n_std_outliers: int = 0
    sustained_outliers: float = np.nan
    model_name: str = ""
    model_dependent_values: Optional[PhaseModelDependentAdler | AvgMagModelDependentAdler] = None
    source_flags: Optional[AdlerSourceFlags] = None
