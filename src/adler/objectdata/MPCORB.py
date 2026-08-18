import numpy as np
from dataclasses import dataclass

from adler.objectdata.objectdata_utilities import get_from_table, get_from_dictionary

MPCORB_KEYS = {
    "mpcDesignation": str,
    "fullDesignation": str,
    "mpcNumber": int,
    "mpcH": float,
    "mpcG": float,
    "epoch": float,
    "tperi": float,
    "peri": float,
    "node": float,
    "incl": float,
    "e": float,
    "n": float,
    "q": float,
    "uncertaintyParameter": str,
}


@dataclass
class MPCORB:
    """Object information from MPCORB. All attributes carry the same names as the column names from the MPCORB table.

    Attributes:
    -----------
    ssObjectId: str
        LSST unique identifier (if observed by LSST)

    mpcDesignation: str
        Number or provisional designation (in packed form)

    fullDesignation: str
        Number or provisional designation (in readable form)

    mpcNumber: int
        MPC number (if the asteroid has been numbered; NULL otherwise). Provided for convenience.

    mpcH: float
        Absolute magnitude, H

    mpcG: float
        Slope parameter, G

    epoch: float
        Epoch (in MJD, .0 TT)

    tperi: float
        MJD of pericentric passage

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

    """

    ssObjectId: str = ""
    mpcDesignation: str = ""
    fullDesignation: str = ""
    mpcNumber: int = 0
    mpcH: float = 0.0
    mpcG: float = 0.0
    epoch: float = 0.0
    tperi: float = 0.0
    peri: float = 0.0
    node: float = 0.0
    incl: float = 0.0
    e: float = 0.0
    n: float = 0.0
    q: float = 0.0
    uncertaintyParameter: str = ""

    # TODO: can these defaults be used to more efficiently set missing values? See SSObject FilterDependentSSO

    @classmethod
    def construct_from_data_table(cls, ssObjectId, data_table):
        """Initialises the MPCORB object from a table of data.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        data_table : table-like object
            Table of data from which attributes should be populated.

        Returns
        -----------
        MPCORB object
            MPCORB object with class attributes populated from data_table.

        """

        mpcorb_dict = {"ssObjectId": ssObjectId}

        for mpcorb_key, mpcorb_type in MPCORB_KEYS.items():
            # add null values if they don't exist
            if mpcorb_key in data_table.to_table().colnames:
                mpcorb_dict[mpcorb_key] = get_from_table(data_table, mpcorb_key, mpcorb_type, "MPCORB")
            else:
                if mpcorb_type == str:
                    mpcorb_dict[mpcorb_key] = ""  # blank string
                elif mpcorb_type == int:
                    mpcorb_dict[mpcorb_key] = (
                        0  # use zero to avoid compatiblity between np.nan and int columns
                    )
                else:
                    mpcorb_dict[mpcorb_key] = np.nan  # set everything else (floats) to np.nan

                # TODO: uncertaintyParameter is ending up as int?

        return cls(**mpcorb_dict)

    @classmethod
    def construct_from_dictionary(cls, ssObjectId, data_dict):
        """Initialises the MPCORB object from a dictionary of data.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        data_dict : dict or dict-like object
            Dictionary of data from which attributes shoud be populated.

        Returns
        -----------
        MPCORB object
            MPCORB object with class attributes populated from data_table.

        """
        mpcorb_dict = {"ssObjectId": ssObjectId}

        for mpcorb_key, mpcorb_type in MPCORB_KEYS.items():
            mpcorb_dict[mpcorb_key] = get_from_dictionary(
                data_dict, mpcorb_key.casefold(), mpcorb_type, "MPCORB"
            )

        return cls(**mpcorb_dict)
