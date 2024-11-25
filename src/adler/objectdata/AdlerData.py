import os
import sqlite3
import logging
import re
import numpy as np
from dataclasses import dataclass, field
from astropy.time import Time


FILTER_DEPENDENT_KEYS = ["phaseAngle_min", "phaseAngle_range", "nobs", "arc"]
MODEL_DEPENDENT_KEYS = [
    "H",
    "H_err",
    "phase_parameter_1",
    "phase_parameter_1_err",
    "phase_parameter_2",
    "phase_parameter_2_err",
    "modelFitMjd",
]
ALL_FILTER_LIST = ["u", "g", "r", "i", "z", "y"]

logger = logging.getLogger(__name__)


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
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        # if model-dependent parameters exist without a model name, return an error
        if not kwargs.get("model_name") and any(name in kwargs for name in MODEL_DEPENDENT_KEYS):
            logger.error("NameError: No model name given. Cannot update model-specific phase parameters.")
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
            if model_key in kwargs:
                setattr(
                    self.filter_dependent_values[filter_index].model_dependent_values[model_index],
                    model_key,
                    kwargs.get(model_key),
                )

    def populate_from_database(self, filepath):
        """Populates the AdlerData object with information from the most recent timestamped entry for the ssObjectId in a given database.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database. Note that for now, we assume only one table with all the data.
        """

        con = self._get_database_connection(filepath)
        cursor = con.cursor()
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

            self.populate_phase_parameters(filter_name, **filter_dependent_info)

            r = re.compile("^(" + filter_name + "_).*_H$")
            model_column_list = list(filter(r.match, column_list))
            models_in_filter = [model[2:-2] for model in model_column_list]

            for model_name in models_in_filter:
                expected_model_columns = [
                    filter_name + "_" + model_name + "_" + model_key for model_key in MODEL_DEPENDENT_KEYS
                ]
                model_indices_list = [
                    column_list.index(column_name) for column_name in expected_model_columns
                ]
                model_values = [fetched_data[a] for a in model_indices_list]
                model_dependent_info = dict(zip(MODEL_DEPENDENT_KEYS, model_values))
                model_dependent_info["model_name"] = model_name

                self.populate_phase_parameters(filter_name, **model_dependent_info)

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
            logger.error("ValueError: Filter {} does not exist in AdlerData.filter_list.".format(filter_name))
            raise ValueError("Filter {} does not exist in AdlerData.filter_list.".format(filter_name))

        output_obj = PhaseParameterOutput()
        output_obj.filter_name = filter_name
        output_obj.phaseAngle_min = self.filter_dependent_values[filter_index].phaseAngle_min
        output_obj.phaseAngle_range = self.filter_dependent_values[filter_index].phaseAngle_range
        output_obj.nobs = self.filter_dependent_values[filter_index].nobs
        output_obj.arc = self.filter_dependent_values[filter_index].arc

        if not model_name:
            logger.warn("No model name was specified. Returning non-model-dependent phase parameters.")
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
            cur.execute("CREATE TABLE AdlerData(ssObjectId INTEGER PRIMARY KEY, timestamp TEXT)")
        elif not database_exists and not create_new:
            logger.error("ValueError: Database cannot be found at given filepath.")
            raise ValueError("Database cannot be found at given filepath.")
        else:
            con = sqlite3.connect(filepath)

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

    def _get_row_data_and_columns(self):
        """Collects all of the data present in the AdlerData object as a list with a corresponding list of column names,
        in preparation for a row to be written to a SQL database table.

        Returns
        -----------
        row_data : list
            A list containing all of the relevant data present in the AdlerData object.

        required_columns : list of str
            A list of the corresponding column names in the same order.

        """
        required_columns = ["ssObjectId", "timestamp"]
        row_data = [int(self.ssObjectId), Time.now().mjd]

        for f, filter_name in enumerate(self.filter_list):
            columns_by_filter = ["_".join([filter_name, filter_key]) for filter_key in FILTER_DEPENDENT_KEYS]
            data_by_filter = [
                getattr(self.filter_dependent_values[f], filter_key) for filter_key in FILTER_DEPENDENT_KEYS
            ]

            required_columns.extend(columns_by_filter)
            row_data.extend(data_by_filter)

            for m, model_name in enumerate(self.filter_dependent_values[f].model_list):
                columns_by_model = [
                    "_".join([filter_name, model_name, model_key]) for model_key in MODEL_DEPENDENT_KEYS
                ]
                data_by_model = [
                    getattr(self.filter_dependent_values[f].model_dependent_values[m], model_key)
                    for model_key in MODEL_DEPENDENT_KEYS
                ]

                required_columns.extend(columns_by_model)
                row_data.extend(data_by_model)

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

    def write_row_to_database(self, filepath, table_name="AdlerData"):
        """Writes all of the relevant data contained within the AdlerData object to a timestamped row in a SQLite database.

        Parameters
        -----------
        filepath : path-like object
            Filepath with the location of the output SQL database.

        table_name : str, optiona
            String containing the table name to write the data to. Default is "AdlerData".

        """

        con = self._get_database_connection(filepath, create_new=True)

        row_data, required_columns = self._get_row_data_and_columns()
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
    modelFitMjd: float = np.nan


class PhaseParameterOutput:
    """Empty convenience class so that the output of AdlerData.get_phase_parameters_in_filter is an object."""

    pass
