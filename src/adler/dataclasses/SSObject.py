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

    H: array_like of floats
        Best fit absolute magnitudes per-filter, in same order as filter_list.

    G12: array_like of floats
        Best fit G12 slope parameters per-filter, in same order as filter_list.

    Herr : array_like of floats
        Uncertainty of H per-filter, in same order as filter_list.

    G12Err : array_like of floats
        Uncertainty of G12 per-filter, in same order as filter_list.

    nData: array_like of ints
        The number of data points used to fit the phase curve per-filter, in same order as filter_list.

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
    H: np.ndarray = field(default_factory=lambda: np.zeros(0))
    G12: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Herr: np.ndarray = field(default_factory=lambda: np.zeros(0))
    G12err: np.ndarray = field(default_factory=lambda: np.zeros(0))
    nData: np.ndarray = field(default_factory=lambda: np.zeros(0))
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

        for i, filter in enumerate(filter_list):
            H[i] = get_from_table(data_table, filter + "_H", "float")
            G12[i] = get_from_table(data_table, filter + "_G12", "float")
            Herr[i] = get_from_table(data_table, filter + "_HErr", "float")
            G12err[i] = get_from_table(data_table, filter + "_G12Err", "float")
            nData[i] = get_from_table(data_table, filter + "_Ndata", "int")

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
            H,
            G12,
            Herr,
            G12err,
            nData,
            maxExtendedness,
            minExtendedness,
            medianExtendedness,
        )
