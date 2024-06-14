from dataclasses import dataclass, field
import numpy as np

from adler.dataclasses.dataclass_utilities import get_from_table

OBSERVATIONS_KEYS = {
    "mag": np.ndarray,
    "magErr": np.ndarray,
    "midPointMjdTai": np.ndarray,
    "ra": np.ndarray,
    "dec": np.ndarray,
    "phaseAngle": np.ndarray,
    "topocentricDist": np.ndarray,
    "heliocentricDist": np.ndarray,
    "heliocentricX": np.ndarray,
    "heliocentricY": np.ndarray,
    "heliocentricZ": np.ndarray,
    "topocentricX": np.ndarray,
    "topocentricY": np.ndarray,
    "topocentricZ": np.ndarray,
    "eclipticLambda": np.ndarray,
    "eclipticBeta": np.ndarray,
}


@dataclass
class Observations:
    """A SQL join of DiaSource and SSSource which contains all of the
    observations of the object in a select filter. All attributes carry
    the same names as the column names from the DiaSource and SSSource tables.

    Attributes:
    -----------
    ssObjectId: str
        Id of the ssObject this source was associated with, if any. If not, it is set to NULL.

    filter_name : str
        Filter of the observations.

    mag: array_like of floats
        Magnitude. This is a placeholder and will be replaced by flux.

    magErr: array_like of floats
        Magnitude error. This is a placeholder and will be replaced by flux error.

    midPointMjdTai: array_like of floats
        Effective mid-visit time for this diaSource, expressed as Modified Julian Date, International Atomic Time.

    ra: array_like of floats
        Right ascension coordinate of the center of this diaSource.

    dec: array_like of floats
        Declination coordinate of the center of this diaSource.

    phaseAngle: array_like of floats
        Phase angle.

    topocentricDist: array_like of floats
        Topocentric distance.

    heliocentricDist: array_like of floats
        Heliocentric distance.

    heliocentricX: array_like of floats
        x-axis component of the heliocentric distance.

    heliocentricY: array_like of floats
        y-axis component of the heliocentric distance.

    heliocentricZ: array_like of floats
        z-axis component of the heliocentric distance.

    topocentricX: array_like of floats
        x-axis component of the topocentric distance.

    topocentricY: array_like of floats
        y-axis component of the topocentric distance.

    topocentricZ: array_like of floats
        z-axis component of the topocentric distance.

    eclipticLambda: array_like of floats
        The ecliptic longitude.

    eclipticBeta: array_like of floats
        The ecliptic latitude.

    reduced_mag: array_like of floats
        The reduced magnitude.

    num_obs : int
        The number of observations contained in this structure.

    """

    ssObjectId: str = ""
    filter_name: str = ""
    mag: np.ndarray = field(default_factory=lambda: np.zeros(0))
    magErr: np.ndarray = field(default_factory=lambda: np.zeros(0))
    midPointMjdTai: np.ndarray = field(default_factory=lambda: np.zeros(0))
    ra: np.ndarray = field(default_factory=lambda: np.zeros(0))
    dec: np.ndarray = field(default_factory=lambda: np.zeros(0))
    phaseAngle: np.ndarray = field(default_factory=lambda: np.zeros(0))
    topocentricDist: np.ndarray = field(default_factory=lambda: np.zeros(0))
    heliocentricDist: np.ndarray = field(default_factory=lambda: np.zeros(0))
    heliocentricX: np.ndarray = field(default_factory=lambda: np.zeros(0))
    heliocentricY: np.ndarray = field(default_factory=lambda: np.zeros(0))
    heliocentricZ: np.ndarray = field(default_factory=lambda: np.zeros(0))
    topocentricX: np.ndarray = field(default_factory=lambda: np.zeros(0))
    topocentricY: np.ndarray = field(default_factory=lambda: np.zeros(0))
    topocentricZ: np.ndarray = field(default_factory=lambda: np.zeros(0))
    eclipticLambda: np.ndarray = field(default_factory=lambda: np.zeros(0))
    eclipticBeta: np.ndarray = field(default_factory=lambda: np.zeros(0))
    reduced_mag: np.ndarray = field(default_factory=lambda: np.zeros(0))
    num_obs: int = 0

    @classmethod
    def construct_from_data_table(cls, ssObjectId, filter_name, data_table):
        """Initialises the Observations object from a table of data.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_name : str
            String of the filter the observations are taken in,

        data_table : table-like object
            Table of data from which attributes shoud be populated.

        Returns
        -----------
        Observations object
            Observations object with class attributes populated from data_table.

        """

        obs_dict = {"ssObjectId": ssObjectId, "filter_name": filter_name, "num_obs": len(data_table)}

        for obs_key, obs_type in OBSERVATIONS_KEYS.items():
            obs_dict[obs_key] = get_from_table(data_table, obs_key, obs_type, "SSSource/DIASource")

        obs_dict["reduced_mag"] = cls.calculate_reduced_mag(
            cls, obs_dict["mag"], obs_dict["topocentricDist"], obs_dict["heliocentricDist"]
        )

        return cls(**obs_dict)

    def calculate_reduced_mag(self, mag, topocentric_dist, heliocentric_dist):
        """
        Calculates the reduced magnitude column.


        Parameters
        -----------
        mag : array_like of floats
            Magnitude.

        topocentric_dist : array_like of floats
            Topocentric distance.

        heliocentric_dist: array_like of floats
            Heliocentric distance.


        Returns
        -----------
        array_like of floats
            The reduced magnitude.


        """
        return mag - 5 * np.log10(topocentric_dist * heliocentric_dist)
