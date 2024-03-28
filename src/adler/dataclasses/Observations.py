from dataclasses import dataclass, field
import numpy as np

from adler.dataclasses.dataclass_utilities import get_from_table


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

    midpointMjdTai: array_like of floats
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

    reduced_mag: array_like of floats
        The reduced magnitude.

    num_obs : int
        The number of observations contained in this structure.

    """

    ssObjectId: str = ""
    filter_name: str = ""
    mag: np.ndarray = field(default_factory=lambda: np.zeros(0))
    magErr: np.ndarray = field(default_factory=lambda: np.zeros(0))
    midpointMjdTai: np.ndarray = field(default_factory=lambda: np.zeros(0))
    ra: np.ndarray = field(default_factory=lambda: np.zeros(0))
    dec: np.ndarray = field(default_factory=lambda: np.zeros(0))
    phaseAngle: np.ndarray = field(default_factory=lambda: np.zeros(0))
    topocentricDist: np.ndarray = field(default_factory=lambda: np.zeros(0))
    heliocentricDist: np.ndarray = field(default_factory=lambda: np.zeros(0))
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

        mag = get_from_table(data_table, "mag", "array")
        magErr = get_from_table(data_table, "magErr", "array")
        midpointMjdTai = get_from_table(data_table, "midPointMjdTai", "array")
        ra = get_from_table(data_table, "ra", "array")
        dec = get_from_table(data_table, "dec", "array")
        phaseAngle = get_from_table(data_table, "phaseAngle", "array")
        topocentricDist = get_from_table(data_table, "topocentricDist", "array")
        heliocentricDist = get_from_table(data_table, "heliocentricDist", "array")

        reduced_mag = cls.calculate_reduced_mag(cls, mag, topocentricDist, heliocentricDist)

        return cls(
            ssObjectId,
            filter_name,
            mag,
            magErr,
            midpointMjdTai,
            ra,
            dec,
            phaseAngle,
            topocentricDist,
            heliocentricDist,
            reduced_mag,
            len(data_table),
        )

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
