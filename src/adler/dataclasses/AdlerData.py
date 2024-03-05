from dataclasses import dataclass, field
import numpy as np


@dataclass
class AdlerData:
    """
    Class for storing Adler-calculated values.

    Note that for all per-filter attributes, the type is an array in order u, g, r, i, z, y.

    Attributes:
    -----------
    phaseAngle_min_adler : array_like
        Minimum phase angle of observations used in fitting model (degrees)

    phaseAngle_range_adler : array_like
        Max minus min phase angle range of observations used in fitting model (degrees)

    nobs_adler : array_like
        Number of observations used in fitting model

    arc_adler : array_like
        Observational arc used to fit model (days)

    H_P16_adler : array_like
        Absolute magnitude of the Penttila et al. 2016 model (mag)

    G12_P16_adler : array_like
        Phase parameter of the Penttila et al. 2016 model

    HErr_P16_adler : array_like
        Uncertainty in absolute magnitude (mag)

    G12Err_P16_adler : array_like
        Uncertainty in phase parameter

    """

    phaseAngle_min_adler: field(default_factory=np.zeros(6))
    phaseAngle_range_adler: field(default_factory=np.zeros(6))
    nobs_adler: field(default_factory=np.zeros(6))
    arc_adler: field(default_factory=np.zeros(6))
    H_P16_adler: field(default_factory=np.zeros(6))
    G12_P16_adler: field(default_factory=np.zeros(6))
    HErr_P16_adler: field(default_factory=np.zeros(6))
    G12Err_P16_adler: field(default_factory=np.zeros(6))
