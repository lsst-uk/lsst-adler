from dataclasses import dataclass, field
import numpy as np


@dataclass
class AdlerData:
    """
    Class for storing Adler-calculated values.

    Note that for all per-filter attributes, the type is an array in order u, g, r, i, z, y.

    """

    phaseAngle_min_adler: field(default_factory=np.zeros(6))
    phaseAngle_range_adler: field(default_factory=np.zeros(6))
    nobs_adler: field(default_factory=np.zeros(6))
    arc_adler: field(default_factory=np.zeros(6))
    H_P16_adler: field(default_factory=np.zeros(6))
    G12_P16_adler: field(default_factory=np.zeros(6))
    HErr_P16_adler: field(default_factory=np.zeros(6))
    G12Err_P16_adler: field(default_factory=np.zeros(6))
