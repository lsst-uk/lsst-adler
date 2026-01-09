import os
import sqlite3
import logging
import re
import numpy as np
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
    "sustained_outliers",
]
PHASE_MODEL_DEPENDENT_KEYS = [
    "H",
    "H_err",
    "phase_parameter_1",
    "phase_parameter_1_err",
    "phase_parameter_2",
    "phase_parameter_2_err",
    "modelFitMjd",
]
VALID_PHASE_MODELS = ["HG", "HG1G2", "HG12", "HG12_Pen16", "LinearPhaseFunc"]
AVG_MAG_MODEL_DEPENDENT_KEYS = [
    "avg_mag",
    "std_mag",
    "modelFitMjd",
]
VALID_AVG_MAG_MODELS = ["median", "mean"]
ALL_FILTER_LIST = ["u", "g", "r", "i", "z", "y"]

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

    filter_dependent_values : list of FilterDependentAdler objects, optional
        List of FilterDependentAdler objects containing filter-dependent data in order of filter_list. Default empty list.

    """

    ssObjectId: str
    filter_list: list

    modelId: str = (
        ""  # TODO need a decision on whether this should be a list or what, how will this interact, are we still going to have lists of modelDepedentAdlers
    )
    filter_dependent_values: list = field(default_factory=list)

    def __post_init__(self):
        """This runs post-initialisation and creates the class attribute where one dimension is "filters" to ensure the array
        has the correct size. This makes population a little easier.
        """

        # note that we don't do the same for model-dependent values as we don't know a priori how many models the user wishes
        # to calculate, but we do know how many filters the AdlerPlanetoid object was generated with
        self.filter_dependent_values = [FilterDependentAdler(filter_name) for filter_name in self.filter_list]

    def populate_filter_dependent_parameters(self, filter_name, **kwargs):
        """#TODO docstring"""

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

    def populate_phase_parameters(self, filter_name, **kwargs):
        # TODO fix docstring
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
        # TODO need error handling or clear message of overwriting as this stands
        if not kwargs.get("model_name"):
            return
        elif kwargs.get("model_name") != self.filter_dependent_values[filter_index].model_name:
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

    def populate_avg_mag_parameters(self, filter_name, **kwargs):
        """#TODO docstring"""

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

    # TODO develop this function further
    def populate_source_flags(self, filter_name, modelId, df, **kwargs):
        """#TODO docstring"""

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
        kwargs.update(
            {"n_outliers": len(df)}
        )  # Add n_outliers to kwargs so it populates FilterDependentAdler
        self.populate_filter_dependent_parameters(filter_name, **kwargs)

        self.filter_dependent_values[filter_index].source_flags = AdlerSourceFlags.construct_from_data_table(
            self.ssObjectId, filter_name, modelId, df
        )

    # TODO figure out how to edit this with the new modelId and no lists
    def populate_from_database(self, filepath):
        """Populates the AdlerData object with information from the most recent timestamped entry for the ssObjectId in a given database.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database. Note that for now, we assume only one table with all the data.
        """

        con = self._get_database_connection(filepath)
        cursor = con.cursor()
        # TODO potentially populate on self.modelId? might have to be a try/except or option of whether it's by ssObjectId or modelId
        # modelId would need to be supplied
        sql_query = f"""SELECT * from AdlerData where ssObjectId='{self.ssObjectId}' ORDER BY timestamp DESC LIMIT 1"""
        query_result = cursor.execute(sql_query)

        try:
            fetched_data_raw = query_result.fetchall()[0]
        except IndexError:
            logger.error("ValueError: No data found in this database for the supplied ssObjectId.")
            raise ValueError("No data found in this database for the supplied ssObjectId.")

        fetched_data = [np.nan if v is None else v for v in fetched_data_raw]  # replaces Nones with nans
        column_list = self._get_database_columns(con, "AdlerData")
        con.close()

        filter_bools = [
            any((column_heading.startswith(filter + "_") for column_heading in column_list))
            for filter in ALL_FILTER_LIST
        ]
        database_filter_list = [b for a, b in zip(filter_bools, ALL_FILTER_LIST) if a]

        if not all([requested_filter in database_filter_list for requested_filter in self.filter_list]):
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
            expected_filter_columns = [filter_name + "_" + filter_key for filter_key in FILTER_DEPENDENT_KEYS]
            filter_indices_list = [column_list.index(column_name) for column_name in expected_filter_columns]
            filter_values = [fetched_data[a] for a in filter_indices_list]
            filter_dependent_info = dict(zip(FILTER_DEPENDENT_KEYS, filter_values))

            r = re.compile("^(" + filter_name + "_).*_modelFitMjd$")
            model_column_list = list(filter(r.match, column_list))
            models_in_filter = [model[2:-12] for model in model_column_list]

            self.populate_filter_dependent_parameters(filter_name, **filter_dependent_info)

            # TODO remove loop
            for model_name in models_in_filter:
                if model_name in VALID_PHASE_MODELS:
                    expected_model_columns = [
                        filter_name + "_" + model_name + "_" + model_key
                        for model_key in PHASE_MODEL_DEPENDENT_KEYS
                    ]
                    model_indices_list = [
                        column_list.index(column_name) for column_name in expected_model_columns
                    ]
                    model_values = [fetched_data[a] for a in model_indices_list]
                    model_dependent_info = dict(zip(PHASE_MODEL_DEPENDENT_KEYS, model_values))
                    model_dependent_info["model_name"] = model_name

                    self.populate_phase_parameters(filter_name, **model_dependent_info)
                elif model_name in VALID_AVG_MAG_MODELS:
                    expected_model_columns = [
                        filter_name + "_" + model_name + "_" + model_key
                        for model_key in AVG_MAG_MODEL_DEPENDENT_KEYS
                    ]
                    model_indices_list = [
                        column_list.index(column_name) for column_name in expected_model_columns
                    ]
                    model_values = [fetched_data[a] for a in model_indices_list]
                    model_dependent_info = dict(zip(AVG_MAG_MODEL_DEPENDENT_KEYS, model_values))
                    model_dependent_info["model_name"] = model_name

                    self.populate_avg_mag_parameters(filter_name, **model_dependent_info)
                else:
                    # TODO improve error message
                    logger.error(
                        f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                    )
                    raise ValueError(
                        f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                    )

    # TODO edits here to print all new info
    def print_data(self):
        """Convenience method to clearly print the stored values."""

        for f, filter_name in enumerate(self.filter_list):
            print("Filter: {}".format(filter_name))
            print("Phase angle minimum: {}".format(self.filter_dependent_values[f].phaseAngle_min))
            print("Phase angle range: {}".format(self.filter_dependent_values[f].phaseAngle_range))
            print("Maximum observation time: {}".format(self.filter_dependent_values[f].observationTime_max))
            print("Number of observations: {}".format(self.filter_dependent_values[f].nobs))
            print("Arc: {}".format(self.filter_dependent_values[f].arc))

            for m, model_name in enumerate(self.filter_dependent_values[f].model_list):
                print("Model: {}.".format(model_name))
                if model_name in VALID_PHASE_MODELS:
                    print("\tH: {}".format(self.filter_dependent_values[f].model_dependent_values[m].H))
                    print(
                        "\tH error: {}".format(
                            self.filter_dependent_values[f].model_dependent_values[m].H_err
                        )
                    )
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
                elif model_name in VALID_AVG_MAG_MODELS:
                    print(
                        "\tAverage magnitude {}".format(
                            self.filter_dependent_values[f].model_dependent_values[m].avg_mag
                        )
                    )
                    print(
                        "\tStandard deviation of magnitudes {}".format(
                            self.filter_dependent_values[f].model_dependent_values[m].std_mag
                        )
                    )
                else:
                    # TODO improve error message
                    logger.error(
                        f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                    )
                    raise ValueError(
                        f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                    )
            print("\n")

    # TODO edit for modelid and no lists
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
        output_obj.sustained_outliers = self.filter_dependent_values[filter_index].sustained_outliers

        if not model_name:
            logger.warning("No model name was specified. Returning non-model-dependent phase parameters.")
            print("No model name specified. Returning non-model-dependent phase parameters.")
        else:
            try:
                model_index = self.filter_dependent_values[filter_index].model_list.index(model_name)
            except ValueError:
                logger.error(
                    "ValueError: Model {} does not exist for filter {} in AdlerData.model_lists.".format(
                        model_name, filter_name
                    )
                )
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

    # TODO edit for modelid and no lists
    def get_avg_mag_parameters_in_filter(self, filter_name, model_name=None):
        # TODO docstring edits/check
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
        output_obj.sustained_outliers = self.filter_dependent_values[filter_index].sustained_outliers

        if not model_name:
            logger.warning(
                "No model name was specified. Returning non-model-dependent average magnitude parameters."
            )
            print("No model name specified. Returning non-model-dependent average magnitude parameters.")
        else:
            try:
                model_index = self.filter_dependent_values[filter_index].model_list.index(model_name)
            except ValueError:
                logger.error(
                    "ValueError: Model {} does not exist for filter {} in AdlerData.model_lists.".format(
                        model_name, filter_name
                    )
                )
                raise ValueError(
                    "Model {} does not exist for filter {} in AdlerData.model_lists.".format(
                        model_name, filter_name
                    )
                )

            output_obj.model_name = model_name
            output_obj.avg_mag = (
                self.filter_dependent_values[filter_index].model_dependent_values[model_index].avg_mag
            )
            output_obj.std_mag = (
                self.filter_dependent_values[filter_index].model_dependent_values[model_index].std_mag
            )

        return output_obj

    def _get_database_connection(self, filepath, create_new=False):
        """Returns the connection to the output SQL database, creating it if it does not exist.

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
            cur.execute("CREATE TABLE AdlerData(ssObjectId, modelId PRIMARY KEY, timestamp REAL)")
        elif not database_exists and not create_new:
            logger.error("ValueError: Database cannot be found at given filepath.")
            raise ValueError("Database cannot be found at given filepath.")
        else:
            con = sqlite3.connect(filepath)
            cur = con.cursor()
            # Create the table if it doesn't exist (in case database was created through AdlerSourceFlags)
            cur.execute(
                "CREATE TABLE IF NOT EXISTS AdlerData(ssObjectId, modelId PRIMARY KEY, timestamp REAL)"
            )

        return con

    def _get_database_columns(self, con, tablename="AdlerData"):
        """Gets a list of the current columns in a given table in a SQL database.

        Parameters
        -----------
        con : sqlite3 Connection object
            The connection to the output SQL database.

        tablename : str
            The name of the relevant table in the database. Default is "AdlerData".


        Returns
        ----------
        list of str
            List of current columns existing in the table.

        """

        cur = con.cursor()
        cur.execute(f"""SELECT * from {tablename} where 1=0""")
        return [d[0] for d in cur.description]

    def _get_row_data_and_columns(self, write_model_data=False):
        """Collects all of the data present in the AdlerData object as a list with a corresponding list of column names,
        in preparation for a row to be written to a SQL database table.

        Parameters
        -----------
        write_model_data : Boolean, optional
            A flag to set whether to write out specific model data to AdlerData. Default: False.

        Returns
        -----------
        row_data : list
            A list containing all of the relevant data present in the AdlerData object.

        required_columns : list of str
            A list of the corresponding column names in the same order.

        """
        required_columns = ["ssObjectId", "modelId", "timestamp"]
        row_data = [self.ssObjectId, self.modelId, Time.now().mjd]

        for f, filter_name in enumerate(self.filter_list):
            columns_by_filter = ["_".join([filter_name, filter_key]) for filter_key in FILTER_DEPENDENT_KEYS]
            data_by_filter = [
                getattr(self.filter_dependent_values[f], filter_key) for filter_key in FILTER_DEPENDENT_KEYS
            ]

            required_columns.extend(columns_by_filter)
            row_data.extend(data_by_filter)
            if write_model_data:
                logger.info(
                    f"write_model_data={write_model_data}, calculated model-specific parameters will be written to AdlerData"
                )
                # TODO remove loop
                for m, model_name in enumerate(self.filter_dependent_values[f].model_list):
                    if model_name in VALID_PHASE_MODELS:
                        columns_by_model = [
                            "_".join([filter_name, model_name, model_key])
                            for model_key in PHASE_MODEL_DEPENDENT_KEYS
                        ]
                        data_by_model = [
                            getattr(self.filter_dependent_values[f].model_dependent_values[m], model_key)
                            for model_key in PHASE_MODEL_DEPENDENT_KEYS
                        ]

                        required_columns.extend(columns_by_model)
                        row_data.extend(data_by_model)
                    elif model_name in VALID_AVG_MAG_MODELS:
                        columns_by_model = [
                            "_".join([filter_name, model_name, model_key])
                            for model_key in AVG_MAG_MODEL_DEPENDENT_KEYS
                        ]
                        data_by_model = [
                            getattr(self.filter_dependent_values[f].model_dependent_values[m], model_key)
                            for model_key in AVG_MAG_MODEL_DEPENDENT_KEYS
                        ]

                        required_columns.extend(columns_by_model)
                        row_data.extend(data_by_model)
                    else:
                        # TODO improve error message
                        logger.error(
                            f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                        )
                        raise ValueError(
                            f"Invalid model name '{model_name}' provided. Model must be one of {VALID_PHASE_MODELS} or {VALID_AVG_MAG_MODELS}"
                        )
            else:
                logger.info(
                    f"write_model_data={write_model_data}, only filter-dependent/model metadata will be written to AdlerData"
                )

        return row_data, required_columns

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

    def write_row_to_database(self, filepath, table_name="AdlerData", write_model_data=False):
        """Writes all of the relevant data contained within the AdlerData object to a timestamped row in a SQLite database.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.

        table_name : str, optiona
            String containing the table name to write the data to. Default is "AdlerData".

        write_model_data : Boolean, optional
            A flag to set whether to write out specific model data to AdlerData. Default: False.

        """

        con = self._get_database_connection(filepath, create_new=True)

        row_data, required_columns = self._get_row_data_and_columns(write_model_data=write_model_data)
        current_columns = self._get_database_columns(con, table_name)
        self._ensure_columns(con, table_name, current_columns, required_columns)

        # TODO edit this (or create new one for metadata, that doesn't overwrite)
        column_names = ",".join(required_columns)
        column_spaces = ",".join(["?"] * len(required_columns))
        update_clause = ", ".join([f"{col} = excluded.{col}" for col in required_columns[1:]])
        sql_command = f"""
                        INSERT INTO {table_name} ({column_names})
                        VALUES ({column_spaces})
                        ON CONFLICT(modelId) DO UPDATE SET {update_clause};
                        """
        # Old command, keeping for now during changes
        # sql_command = f"""
        #                 INSERT INTO {table_name} ({column_names})
        #                 VALUES ({column_spaces})
        #                 ON CONFLICT(ssObjectId) DO UPDATE SET {update_clause};
        #                 """
        cur = con.cursor()
        cur.execute(sql_command, row_data)
        con.commit()
        con.close()


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

    modelFitMjd : float, optional
        The MJD of when the model fit was calculated. Default NaN.

    """

    filter_name: str
    model_name: str
    H: float = np.nan
    H_err: float = np.nan
    phase_parameter_1: float = np.nan
    phase_parameter_1_err: float = np.nan
    phase_parameter_2: float = np.nan
    phase_parameter_2_err: float = np.nan
    modelFitMjd: float = np.nan


@dataclass
class AvgMagModelDependentAdler:
    """Dataclass containing model-dependent values for the simple average magnitude model pgenerated by Adler. Note that NaN indicates a value that has not yet been populated.

    Attributes:
    -----------
    filter_name : str
        The filter for which these values are calculated.

    model_name : str
        The phase model for which these values were calculated. Example: "median", "mean".

    avg_mag : float, optional
        Average magnitude of the measurements used to calculate the model. Default NaN.

    std_mag : float, optional
        Standard deviation of the measurements used to calculate the model. Default NaN.

    modelFitMjd : float, optional
        The MJD of when the model fit was calculated. Default NaN.

    """

    filter_name: str
    model_name: str
    avg_mag: float = np.nan
    std_mag: float = np.nan
    modelFitMjd: float = np.nan


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

    modelId : str, optional
        modelId for the model that the outliers are compared to. Default: Empty str

    n_outliers : int
        Number of observations identified as outliers

    diaSourceId : array_like of ints or strs
        Unique identifier of the observation.

    mag_diff : array_like of floats
        Differences in (reduced) magnitude between the observations and the model. Default: np.nan

    std_diff : array_like of floats
        Differences in sigma-space between the observations and the model. Default: np.nan

    """

    ssObjectId: str
    filter_name: str
    modelId: str
    n_outliers: int  # TODO possibly a cleaner interaction between this and the value in AdlerData/FilterDependentAdler
    diaSourceId: np.ndarray = field(default_factory=lambda: np.zeros(0))
    mag_diff: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # std_diff: np.ndarray = field(default_factory=lambda: np.zeros(0))

    # TODO probably make this more like the one in Observations.py for consistency
    @classmethod
    def construct_from_data_table(cls, ssObjectId, filter_name, modelId, df):
        """#TODO docstring"""
        obs_dict = {"ssObjectId": ssObjectId, "filter_name": filter_name, "modelId": modelId}
        obs_dict.update({"n_outliers": len(df)})
        obs_dict.update(df.loc[:, ["diaSourceId", "mag_diff"]].to_dict(orient="list"))

        # TODO check this will still work if there is difference between what is an outlier for mag_diff vs std_diff
        # i.e. a different number of rows (need to make sure to NULL or 0 fill)

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

        # TODO does this need an explicit PRIMARY KEY? Currently the default columns are non-unique
        if not database_exists and create_new:  # we need to make the table and a couple of starter columns
            con = sqlite3.connect(filepath)
            cur = con.cursor()
            cur.execute(
                "CREATE TABLE AdlerSourceFlags(ssObjectId, filter_name, modelId, diaSourceId, mag_diff, std_diff)"
            )
        elif not database_exists and not create_new:
            logger.error("ValueError: Database cannot be found at given filepath.")
            raise ValueError("Database cannot be found at given filepath.")
        else:
            con = sqlite3.connect(filepath)
            cur = con.cursor()
            # Create the table if it doesn't exist (in case database was created through AdlerData)
            cur.execute(
                "CREATE TABLE IF NOT EXISTS AdlerSourceFlags(ssObjectId, filter_name, modelId, diaSourceId, mag_diff, std_diff)"
            )

        return con

    def write_flags_to_database(self, filepath, table_name="AdlerSourceFlags"):
        """#TODO docstring"""

        con = self._get_database_connection(filepath, create_new=True)

        # Don't need to write out n_outliers in this table
        # TODO implement std_diff
        required_columns = [
            "ssObjectId",
            "filter_name",
            "modelId",
            "diaSourceId",
            "mag_diff",
        ]  # , "std_diff"]

        column_names = ",".join(required_columns)
        column_spaces = ",".join(["?"] * len(required_columns))
        sql_command = f"""
                        INSERT INTO {table_name} ({column_names})
                        VALUES ({column_spaces});
                        """

        # TODO implement check for data in mag_diff, std_diff data

        row_data = list(
            zip(
                [self.ssObjectId] * self.n_outliers,
                [self.filter_name] * self.n_outliers,
                [self.modelId] * self.n_outliers,
                self.diaSourceId,
                self.mag_diff,
                # self.std_diff,
            )
        )

        cur = con.cursor()
        cur.executemany(sql_command, row_data)
        con.commit()
        con.close()

        # TODO sqlite ON CONFLICT on modelId and diaSourceId combined?
        #  (modelId can appear multiple times with different diaSourceIds)
        # (diaSourceId can appear multilpe times for different modelIds)
        # but we shouldn't have two rows with the same modelId and diaSourceId


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
        Maximum time of observation used in fitting model (modifided Julian day).

    arc: float, optional
        Observational arc used to fit model (days).

    nobs : int, optional
        Number of observations used in fitting model.

    n_outliers : int, optional
        Number of outliers detected for the given model.

    sustained_outliers : float, optional
        Magnitude difference between old and new observations

    model_list: list of str, optional
        List of the models for which phase curve or average magnitude parameters have been calculated. Default: empty list

    model_dependent_values: list of PhaseModelDependentAdler or AvgMagModelDependentAdler objects, optional
        List of PhaseModelDependentAdler or AvgMagModelDependentAdler objects storing phase-model or average-magnitude-model parameters for each model, given in order of model_list. Default: empty list.

    """

    filter_name: str
    phaseAngle_min: float = np.nan
    phaseAngle_range: float = np.nan
    observationTime_max: float = np.nan
    arc: float = np.nan
    nobs: int = 0
    n_outliers: int = 0
    sustained_outliers: float = np.nan
    model_name: str = ""
    model_dependent_values: Optional[PhaseModelDependentAdler | AvgMagModelDependentAdler] = None
    source_flags: Optional[AdlerSourceFlags] = None
