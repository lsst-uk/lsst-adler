import numpy as np
import sys
from lsst.rsp import get_tap_service


class DataSchema:
    """Parent class for Observations (a join of DiaSource and SSSource), MPCORB
    and SSObject data classes. Largely a collection of common methods. Should never be instantiated
    by itself.
    """

    def populate(self, population_location, sql_query, sql_filename):
        """Populates the DataSchema object, either from the RSP or a SQL table. Note that this calls the methods
        get_RSP_table() and get_SQL_table(), which must exist in the child classes.

        Parameters
        -----------
        population_location : str
            String delineating source of data. Should be "RSP" for Rubin Science Platform or "SQL" for a SQL table.
        sql_query: str
            SQL query to retrieve data from database.

        """

        if population_location == "RSP":  # pragma: no cover
            self.get_RSP_table(sql_query)
        elif population_location == "SQL":
            self.get_SQL_table(sql_query, sql_filename)
        else:
            sys.exit(
                "Population source not recognised. Please supply either 'RSP' or 'SQL' for population_location argument."
            )

        self.populate_from_table()

    def get_RSP_table(self, sql_query):  # pragma: no cover
        """Retrieves the table of data from the RSP. Populates the data_table class variable.

        Parameters
        -----------
        sql_query : str
            SQL query to be sent to the RSP tap service.

        """
        self.sql_query = sql_query
        self.service = get_tap_service("ssotap")

        self.data_table = self.service.search(sql_query).to_table()

    def get_SQL_table(self, sql_query, sql_filename):
        pass

    def get_from_table(self, column_name, type):
        """Retrieves information from the data_table class variable and forces it to be a specified type.

        Parameters
        -----------
        column_name : str
            Column name under which the data of interest is stored.
        type : str
            String delineating data type. Should be "str", "float", "int" or "array".

        Returns
        -----------
        data : any type
            The data requested from the table cast to the type required.

        """
        try:
            if type == "str":
                return str(self.data_table[column_name][0])
            elif type == "float":
                return float(self.data_table[column_name][0])
            elif type == "int":
                return int(self.data_table[column_name][0])
            elif type == "array":
                return np.array(self.data_table[column_name])
            else:
                print("Type not recognised.")
        except ValueError:
            sys.exit("Could not cast column name to type.")


class Observations(DataSchema):
    """This is a SQL join of DiaSource and SSSource which contains all of the
    observations of the object. Inherits from DataSchema.
    """

    def __init__(self, ssObjectId, population_location, sql_query, sql_filename=None):
        """Initiates the Observations object.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.
        population_location : str
            String delineating source of data. Should be "RSP" for Rubin Science Platform, "SQL" for a SQL table,
            or "arguments" for arguments.
        sql_query : str
            SQL query to retrieve data from database.
        sql_filename: str, optional
            Location of local SQL database, if using.

        """

        self.ssObjectId = ssObjectId
        self.population_location = population_location
        self.populate(self.population_location, sql_query, sql_filename)
        self.calculate_reduced_mag()

    def populate_from_table(self):
        """Populates the Observations object from the data_table class variable created on initialisation."""

        self.mag = self.get_from_table("mag", "array")
        self.magErr = self.get_from_table("magErr", "array")
        self.mjd = self.get_from_table("mjd", "array")
        self.ra = self.get_from_table("ra", "array")
        self.dec = self.get_from_table("dec", "array")
        self.phaseAngle = self.get_from_table("phaseAngle", "array")
        self.topocentricDist = self.get_from_table("topocentricDist", "array")
        self.heliocentricDist = self.get_from_table("heliocentricDist", "array")

    def calculate_reduced_mag(self):
        """
        Calculates the reduced magnitude column.
        """
        self.reduced_mag = self.mag - 5 * np.log10(self.topocentricDist * self.heliocentricDist)


class MPCORB(DataSchema):
    """Grabs information from MPCORB."""

    def __init__(self, ssObjectId, population_location, sql_query, sql_filename):
        """Initiates the MPCORB object.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.
        population_location : str
            String delineating source of data. Should be "RSP" for Rubin Science Platform or "SQL" for a SQL table.
        sql_query : str
            SQL query to retrieve data from database.
        sql_filename: str, optional
            Location of local SQL database, if using.
        """

        self.ssObjectId = ssObjectId
        self.population_location = population_location
        self.populate(self.population_location, sql_query, sql_filename)

    def populate_from_table(self):
        """Populates the MPCORB object from the data_table class variable created on initialisation."""

        self.mpcDesignation = self.get_from_table("mpcDesignation", "str")
        self.mpcNumber = self.get_from_table("mpcNumber", "str")
        self.mpcH = self.get_from_table("mpcH", "float")
        self.mpcG = self.get_from_table("mpcH", "float")
        self.epoch = self.get_from_table("epoch", "float")
        self.peri = self.get_from_table("peri", "float")
        self.node = self.get_from_table("node", "float")
        self.incl = self.get_from_table("incl", "float")
        self.e = self.get_from_table("e", "float")
        self.n = self.get_from_table("n", "float")
        self.q = self.get_from_table("q", "float")
        self.uncertaintyParameter = self.get_from_table("uncertaintyParameter", "str")
        self.flags = self.get_from_table("flags", "str")


class SSObject(DataSchema):
    """Grabs information from SSObject."""

    def __init__(self, ssObjectId, population_location, sql_query, sql_filename):
        """Initiates the SSObject object.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.
        population_location : str
            String delineating source of data. Should be "RSP" for Rubin Science Platform or "SQL" for a SQL table.
        sql_query : str
            SQL query to retrieve data from database.
        sql_filename: str, optional
            Location of local SQL database, if using.
        """

        self.ssObjectId = ssObjectId
        self.population_location = population_location
        self.populate(self.population_location, sql_query, sql_filename)

    def populate_from_table(self):
        """Populates the SSObject object from the data_table class variable created on initialisation."""

        self.discoverySubmissionDate = self.get_from_table("discoverySubmissionDate", "float")
        self.firstObservationDate = self.get_from_table("firstObservationDate", "float")
        self.arc = self.get_from_table("arc", "float")
        self.numObs = self.get_from_table("numObs", "int")
        self.r_H = self.get_from_table("r_H", "float")
        self.r_G12 = self.get_from_table("r_G12", "float")
        self.r_Herr = self.get_from_table("r_Herr", "float")
        self.r_G12Err = self.get_from_table("r_G12err", "float")
        self.r_nData = self.get_from_table("r_nData", "int")
        self.maxExtendedness = self.get_from_table("maxExtendedness", "float")
        self.minExtendedness = self.get_from_table("minExtendedness", "float")
        self.medianExtendedness = self.get_from_table("medianExtendedness", "float")
