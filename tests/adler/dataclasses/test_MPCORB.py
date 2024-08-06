from numpy.testing import assert_almost_equal
import numpy as np

from adler.dataclasses.MPCORB import MPCORB
from adler.dataclasses.dataclass_utilities import get_data_table
from adler.utilities.tests_utilities import get_test_data_filepath


def test_construct_MPCORB_from_data_table():
    ssoid = 8268570668335894776
    test_db_path = get_test_data_filepath("testing_database.db")

    test_query = f"""
                SELECT
                    ssObjectId, mpcDesignation, fullDesignation, mpcNumber, mpcH, mpcG, epoch, tperi, peri, node, incl, e, n, q, 
                    uncertaintyParameter, flags
                FROM
                    MPCORB
                WHERE
                    ssObjectId = {ssoid}
            """

    data_table = get_data_table(test_query, sql_filename=test_db_path)
    test_MPCORB = MPCORB.construct_from_data_table(ssoid, data_table)

    assert test_MPCORB.ssObjectId == 8268570668335894776
    assert test_MPCORB.mpcDesignation == "2014 QL4"
    assert (
        test_MPCORB.fullDesignation == "2011 2014 QL433"
    )  # N.B. that there are known DP0.3 issues with mpcDesignation and fullDesignation, https://dp0-3.lsst.io/data-products-dp0-3/data-simulation-dp0-3.html#known-issues
    assert test_MPCORB.mpcNumber == 0
    assert_almost_equal(test_MPCORB.mpcH, 19.8799991607666, decimal=6)
    assert_almost_equal(test_MPCORB.mpcG, 0.15000000596046448, decimal=6)
    assert_almost_equal(test_MPCORB.epoch, 60065.0, decimal=6)
    assert_almost_equal(test_MPCORB.peri, 260.5468204162153, decimal=6)
    assert_almost_equal(test_MPCORB.node, 322.8059, decimal=6)
    assert_almost_equal(test_MPCORB.incl, 4.427569999999975, decimal=6)
    assert_almost_equal(test_MPCORB.e, 0.7168805704972735, decimal=6)
    assert_almost_equal(test_MPCORB.e, 0.7168805704972735, decimal=6)
    assert np.isnan(test_MPCORB.n)
    assert_almost_equal(test_MPCORB.q, 0.5898291078470536, decimal=6)
    assert np.isnan(test_MPCORB.uncertaintyParameter)
    assert test_MPCORB.flags == "0"
