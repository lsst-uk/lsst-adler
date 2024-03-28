from dataclasses import dataclass, field
import numpy as np

from adler.dataclasses.dataclass_utilities import get_from_table


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

    @classmethod
    def construct_from_data_table(cls, ssObjectId, filter_list, data_table):
        discoverySubmissionDate = get_from_table(data_table, "discoverySubmissionDate", "float")
        firstObservationDate = get_from_table(data_table, "firstObservationDate", "float")
        arc = get_from_table(data_table, "arc", "float")
        numObs = get_from_table(data_table, "numObs", "int")

        H = np.zeros(len(filter_list))
        G12 = np.zeros(len(filter_list))
        Herr = np.zeros(len(filter_list))
        G12err = np.zeros(len(filter_list))
        nData = np.zeros(len(filter_list))

        filter_dependent_values = []

        for i, filter_name in enumerate(filter_list):
            filter_dept_object = FilterDependentSSO(
                filter_name=filter_name,
                H=get_from_table(data_table, filter_name + "_H", "float"),
                G12=get_from_table(data_table, filter_name + "_G12", "float"),
                Herr=get_from_table(data_table, filter_name + "_HErr", "float"),
                G12err=get_from_table(data_table, filter_name + "_G12Err", "float"),
                nData=get_from_table(data_table, filter_name + "_Ndata", "int"),
            )

            filter_dependent_values.append(filter_dept_object)

        maxExtendedness = get_from_table(data_table, "maxExtendedness", "float")
        minExtendedness = get_from_table(data_table, "minExtendedness", "float")
        medianExtendedness = get_from_table(data_table, "medianExtendedness", "float")

        return cls(
            ssObjectId,
            filter_list,
            discoverySubmissionDate,
            firstObservationDate,
            arc,
            numObs,
            filter_dependent_values,
            maxExtendedness,
            minExtendedness,
            medianExtendedness,
        )


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
    H: float
    G12: float
    Herr: float
    G12err: float
    nData: int = 0
