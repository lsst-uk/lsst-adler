from dataclasses import dataclass, field
import numpy as np

from adler.objectdata.objectdata_utilities import get_from_table, get_from_dictionary

# TODO: SSO_KEYS needs editing for DP1
# NB Filter keys not here?
SSO_KEYS = {
    "discoverySubmissionDate": float,
    "firstObservationDate": float,
    "arc": float,
    "numObs": int,
    "maxExtendedness": float,
    "minExtendedness": float,
    "medianExtendedness": float,
}
SSO_FILT_KEYS = {
    "H": float,
    "G12": float,
    "HErr": float,
    "G12Err": float,
    "Ndata": int,
}


@dataclass
class SSObject:
    """Object information from SSObject. All attributes carry the same names as the column names from the SSObject table.

    Attributes:
    -----------

    ssObjectId: str
        LSST unique identifier

    filter_list : list of str
        A comma-separated list of the filters of interest.

    discoverySubmissionDate : float
        The date the LSST first linked and submitted the discovery observations to the MPC. May be NULL if not an LSST discovery. The date format will follow general LSST conventions (MJD TAI, at the moment).

    firstObservationDate: float
        The time of the first LSST observation of this object (could be precovered)

    arc: float
        Arc of LSST observations

    numObs: int
        Number of LSST observations of this object

    filter_dependent_values: list of FilterDependentSSO objects
        A list of FilterDependentSSO objects storing the filter-dependent values H, Herr, G12, G12err and nData,
        in same order as filter_list. See documentation for FilterDependentSSO object for descriptions of these variables.

    maxExtendedness: float
        maximum `extendedness` value from the DIASource

    minExtendedness: float
        minimum `extendedness` value from the DIASource

    medianExtendedness: float
        median `extendedness` value from the DIASource

    """

    ssObjectId: str = ""
    filter_list: list = field(default_factory=list)
    discoverySubmissionDate: float = 0.0
    firstObservationDate: float = 0.0
    arc: float = 0.0
    numObs: int = 0
    filter_dependent_values: list = field(default_factory=list)
    maxExtendedness: float = 0.0
    minExtendedness: float = 0.0
    medianExtendedness: float = 0.0

    # TODO: can these defaults be used to more efficiently set missing values? See SSObject FilterDependentSSO

    @classmethod
    def construct_from_data_table(cls, ssObjectId, filter_list, data_table):
        """Initialises the SSObject object from a table of data.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        data_table : table-like object
            Table of data from which attributes shoud be populated.

        Returns
        -----------
        SSObject object
            SSObject object with class attributes populated from data_table.

        """
        sso_dict = {"ssObjectId": ssObjectId, "filter_list": filter_list, "filter_dependent_values": []}

        for sso_key, sso_type in SSO_KEYS.items():

            # add null values if they don't exist
            if sso_key in data_table.colnames:
                sso_dict[sso_key] = get_from_table(data_table, sso_key, sso_type, "SSObject")
            else:
                if sso_type == str:
                    sso_dict[sso_key] = ""  # blank string
                elif sso_type == int:
                    sso_dict[sso_key] = 0  # use zero to avoid compatiblity between np.nan and int columns
                else:
                    sso_dict[sso_key] = np.nan  # set everything else (floats) to np.nan

        for filter_name in filter_list:

            filt_obj_dict = {"filter_name": filter_name}
            for sso_key, sso_type in SSO_FILT_KEYS.items():

                # add as a key to the dict if it exists
                if (filter_name + "_" + sso_key) in data_table.colnames:
                    filt_obj_dict[sso_key] = get_from_table(
                        data_table, filter_name + "_" + sso_key, sso_type, "SSObject"
                    )

            # create filter object from dict, missing values will go to default
            filter_dept_object = FilterDependentSSO(**filt_obj_dict)
            sso_dict["filter_dependent_values"].append(filter_dept_object)

        return cls(**sso_dict)

    @classmethod
    def construct_from_dictionary(cls, ssObjectId, filter_list, data_dict):
        """Initialises the SSObject object from a dictionary of data.

        Parameters
        -----------
        ssObjectId : str
            ssObjectId of the object of interest.

        filter_list : list of str
            A comma-separated list of the filters of interest.

        data_dict : dict or dict-like object
            Ditcionary of data from which attributes shoud be populated.

        Returns
        -----------
        SSObject object
            SSObject object with class attributes populated from data_dict.

        """
        sso_dict = {"ssObjectId": ssObjectId, "filter_list": filter_list, "filter_dependent_values": []}

        for sso_key, sso_type in SSO_KEYS.items():
            sso_dict[sso_key] = get_from_dictionary(data_dict, sso_key.casefold(), sso_type, "SSObject")

        for i, filter_name in enumerate(filter_list):
            filter_dept_object = FilterDependentSSO(
                filter_name=filter_name,
                H=get_from_dictionary(data_dict, (filter_name + "_H").casefold(), float, "SSObject"),
                G12=get_from_dictionary(data_dict, (filter_name + "_G12").casefold(), float, "SSObject"),
                Herr=get_from_dictionary(data_dict, (filter_name + "_HErr").casefold(), float, "SSObject"),
                G12err=get_from_dictionary(
                    data_dict, (filter_name + "_G12Err").casefold(), float, "SSObject"
                ),
                nData=get_from_dictionary(data_dict, (filter_name + "_Ndata").casefold(), float, "SSObject"),
            )

            sso_dict["filter_dependent_values"].append(filter_dept_object)

        return cls(**sso_dict)


@dataclass
class FilterDependentSSO:
    """Filter-dependent object information from SSObject. All attributes carry the same names as the column names from the SSObject table.

    Attributes:
    -----------
    filter_name : str
        Single-letter name of the filter for which these values are relevant.

    H : float
        Best fit absolute magnitude in filter.

    G12: float
        Best fit G12 slope parameters in filter.

    Herr : float
        Uncertainty of H.

    G12Err : float
        Uncertainty of G12.

    nData: int
        The number of data points used to fit the phase curve in this filter.
    """

    filter_name: str
    H: float = np.nan
    G12: float = np.nan
    HErr: float = np.nan
    G12Err: float = np.nan
    Ndata: int = 0

    # TODO: make names here consistent with adler_schema_map.csv?
