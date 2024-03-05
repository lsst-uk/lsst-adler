import numpy as np
import sys
from dataclasses import dataclass, field
from lsst.rsp import get_tap_service


class DataSchema:
    """Parent class for Observations (a join of DiaSource and SSSource), MPCORB
    and SSObject data classes. Largely a collection of common methods. Should never be instantiated
    by itself.
    """

    def populate(self, population_location, sql_query, sql_filename):
        """Populates the DataSchema object, either from the RSP or a SQL table. Note that this calls the method
        populate_from_table(), which must exist in the child classes.

        Parameters
        -----------
        population_location : str
            String delineating source of data. Should be "RSP" for Rubin Science Platform or "SQL" for a SQL table.
        sql_query: str
            SQL query to retrieve data from database.

        """

        if population_location == "RSP":  # pragma: no cover
            data_table = self.get_RSP_table(sql_query)
        elif population_location == "SQL":
            data_table = self.get_SQL_table(sql_query, sql_filename)
        else:
            sys.exit(
                "Population source not recognised. Please supply either 'RSP' or 'SQL' for population_location argument."
            )

        self.populate_from_table(data_table)

    def get_RSP_table(self, sql_query):  # pragma: no cover
        """Retrieves the table of data from the RSP. Populates the data_table class variable.

        Parameters
        -----------
        sql_query : str
            SQL query to be sent to the RSP tap service.

        """

        self.service = get_tap_service("ssotap")

        return self.service.search(sql_query).to_table()

    def get_SQL_table(self, sql_query, sql_filename):
        pass

    def get_from_table(self, data_table, column_name, data_type):
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
            if data_type == "str":
                return str(data_table[column_name][0])
            elif data_type == "float":
                return float(data_table[column_name][0])
            elif data_type == "int":
                return int(data_table[column_name][0])
            elif data_type == "array":
                return np.array(data_table[column_name])
            else:
                print("Type not recognised.")
        except ValueError:
            sys.exit("Could not cast column name to type.")


@dataclass
class Observations(DataSchema):
    """This is a SQL join of DiaSource and SSSource which contains all of the
    observations of the object. Inherits from DataSchema. All attributes carry
    the same names as the column names from the DiaSource and SSSource tables.

    Attributes:
    -----------

    ssObjectId: str
        Id of the ssObject this source was associated with, if any. If not, it is set to NULL.

    mag: array_like
        Magnitude. This is a placeholder and will be replaced by flux.

    magErr: array_like
        Magnitude error. This is a placeholder and will be replaced by flux error.

    midpointMjdTai: array_like
        Effective mid-visit time for this diaSource, expressed as Modified Julian Date, International Atomic Time.

    ra: array_like
        Right ascension coordinate of the center of this diaSource.

    dec: array_like
        Declination coordinate of the center of this diaSource.

    phaseAngle: array_like
        Phase angle.

    topocentricDist: array_like
        Topocentric distance.

    heliocentricDist: array_like
        Heliocentric distance.

    reduced_mag: array_like
        The reduced magnitude.

    """

    ssObjectId: str = ""
    mag: np.ndarray = field(default_factory=np.zeros(0))
    magErr: np.ndarray = field(default_factory=np.zeros(0))
    midpointMjdTai: np.ndarray = field(default_factory=np.zeros(0))
    ra: np.ndarray = field(default_factory=np.zeros(0))
    dec: np.ndarray = field(default_factory=np.zeros(0))
    phaseAngle: np.ndarray = field(default_factory=np.zeros(0))
    topocentricDist: np.ndarray = field(default_factory=np.zeros(0))
    heliocentricDist: np.ndarray = field(default_factory=np.zeros(0))
    reduced_mag: np.ndarray = field(default_factory=np.zeros(0))

    def __init__(self, ssObjectId, population_location, sql_query, sql_filename=None):
        """Initialises the Observations object.

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

    def populate_from_table(self, data_table):
        """Populates the Observations object from the data_table class variable created on initialisation."""

        self.mag = self.get_from_table(data_table, "mag", "array")
        self.magErr = self.get_from_table(data_table, "magErr", "array")
        self.midpointMjdTai = self.get_from_table(data_table, "midpointMjdTai", "array")
        self.ra = self.get_from_table(data_table, "ra", "array")
        self.dec = self.get_from_table(data_table, "dec", "array")
        self.phaseAngle = self.get_from_table(data_table, "phaseAngle", "array")
        self.topocentricDist = self.get_from_table(data_table, "topocentricDist", "array")
        self.heliocentricDist = self.get_from_table(data_table, "heliocentricDist", "array")

    def calculate_reduced_mag(self):
        """
        Calculates the reduced magnitude column.
        """
        self.reduced_mag = self.mag - 5 * np.log10(self.topocentricDist * self.heliocentricDist)


@dataclass
class MPCORB(DataSchema):
    """Grabs information from MPCORB. All attributes carry the same names as the column names from the MPCORB table.

    Attributes:
    -----------

    ssObjectId: str
        LSST unique identifier (if observed by LSST)

    mpcDesignation: str
        Number or provisional designation (in packed form)

    mpcNumber: int
        MPC number (if the asteroid has been numbered; NULL otherwise). Provided for convenience.

    mpcH: float
        Absolute magnitude, H

    mpcG: float
        Slope parameter, G

    epoch: float
        Epoch (in MJD, .0 TT)

    peri: float
        Argument of perihelion, J2000.0 (degrees)

    node: float
        Longitude of the ascending node, J2000.0 (degrees)

    incl: float
        Inclination to the ecliptic, J2000.0 (degrees)

    e: float
        Orbital eccentricity

    n: float
        Mean daily motion (degrees per day)

    q: float
        Perihelion distance (AU)

    uncertaintyParameter: str
        Uncertainty parameter, U

    flags: str
        4-hexdigit flags. See https://minorplanetcenter.net//iau/info/MPOrbitFormat.html for details

    """

    ssObjectId: str = ""
    mpcDesignation: str = ""
    mpcNumber: int = 0
    mpcH: float = 0.0
    mpcG: float = 0.0
    epoch: float = 0.0
    peri: float = 0.0
    node: float = 0.0
    incl: float = 0.0
    e: float = 0.0
    n: float = 0.0
    q: float = 0.0
    uncertaintyParameter: str = ""
    flags: str = ""

    def __init__(self, ssObjectId, population_location, sql_query, sql_filename):
        """Initialises the MPCORB object.

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

    def populate_from_table(self, data_table):
        """Populates the MPCORB object from the data_table class variable created on initialisation."""

        self.mpcDesignation = self.get_from_table(data_table, "mpcDesignation", "str")
        self.mpcNumber = self.get_from_table(data_table, "mpcNumber", "int")
        self.mpcH = self.get_from_table(data_table, "mpcH", "float")
        self.mpcG = self.get_from_table(data_table, "mpcH", "float")
        self.epoch = self.get_from_table(data_table, "epoch", "float")
        self.peri = self.get_from_table(data_table, "peri", "float")
        self.node = self.get_from_table(data_table, "node", "float")
        self.incl = self.get_from_table(data_table, "incl", "float")
        self.e = self.get_from_table(data_table, "e", "float")
        self.n = self.get_from_table(data_table, "n", "float")
        self.q = self.get_from_table(data_table, "q", "float")
        self.uncertaintyParameter = self.get_from_table(data_table, "uncertaintyParameter", "str")
        self.flags = self.get_from_table(data_table, "flags", "str")


@dataclass
class SSObject(DataSchema):
    """Grabs information from SSObject. All attributes carry the same names as the column names from the SSObject table.

    Attributes:
    -----------

    ssObjectId: str
        LSST unique identifier (if observed by LSST)

    discoverySubmissionDate : float
        The date the LSST first linked and submitted the discovery observations to the MPC. May be NULL if not an LSST discovery. The date format will follow general LSST conventions (MJD TAI, at the moment).

    firstObservationDate: float
        The time of the first LSST observation of this object (could be precovered)

    arc: float
        Arc of LSST observations

    numObs: int
        Number of LSST observations of this object

    r_H: float
        Best fit absolute magnitude (r band)

    r_G12: float
        Best fit G12 slope parameter (r band)

    r_Herr:
        Uncertainty of H (r band)

    r_G12Err:
        Uncertainty of G12 (r band)

    r_nData: int
        The number of data points used to fit the phase curve (r band)

    maxExtendedness: float
        maximum `extendedness` value from the DIASource

    minExtendedness: float
        minimum `extendedness` value from the DIASource

    medianExtendedness: float
        median `extendedness` value from the DIASource

    """

    ssObjectId: str = ""
    discoverySubmissionDate: float = 0.0
    firstObservationDate: float = 0.0
    arc: float = 0.0
    numObs: int = 0
    r_H: float = 0.0
    r_G12: float = 0
    r_Herr: float = 0.0
    r_G12_err: float = 0.0
    r_nData: int = 0
    maxExtendedness: float = 0.0
    minExtendedness: float = 0.0
    medianExtendedness: float = 0.0

    def __init__(self, ssObjectId, population_location, sql_query, sql_filename):
        """Initialises the SSObject object.

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

    def populate_from_table(self, data_table):
        """Populates the SSObject object from the data_table class variable created on initialisation."""

        self.discoverySubmissionDate = self.get_from_table(data_table, "discoverySubmissionDate", "float")
        self.firstObservationDate = self.get_from_table(data_table, "firstObservationDate", "float")
        self.arc = self.get_from_table(data_table, "arc", "float")
        self.numObs = self.get_from_table(data_table, "numObs", "int")
        self.r_H = self.get_from_table(data_table, "r_H", "float")
        self.r_G12 = self.get_from_table(data_table, "r_G12", "float")
        self.r_Herr = self.get_from_table(data_table, "r_Herr", "float")
        self.r_G12Err = self.get_from_table(data_table, "r_G12err", "float")
        self.r_nData = self.get_from_table(data_table, "r_nData", "int")
        self.maxExtendedness = self.get_from_table(data_table, "maxExtendedness", "float")
        self.minExtendedness = self.get_from_table(data_table, "minExtendedness", "float")
        self.medianExtendedness = self.get_from_table(data_table, "medianExtendedness", "float")
