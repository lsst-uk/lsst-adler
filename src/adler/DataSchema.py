import numpy as np


class DataSchema:
    """Parent class for Observations (a join of DiaSource and SSSource), MPCORB
    and SSObject data classes.
    """

    def __init__(self, ssObjectId, sql_query, service, sql_filename=None):
        self.ssObjectId = ssObjectId
        self.sql_query = sql_query
        self.service = service

        if not sql_filename:
            self.data_table = self.get_RSP_table(self.sql_query)
        else:
            self.data_table = self.get_SQL_table(self.sql_query)

    def get_RSP_table(self, sql_query):
        rsp_table = self.service.search(sql_query).to_table()
        return rsp_table

    def get_SQL_table(self, sql_query, testing_filename):
        pass

    # should be one function to get whatever from the table and type accordingly
    def get_array_from_table(self, column_name):
        return np.array(self.data_table[column_name])

    def get_string_from_table(self, column_name):
        return str(self.data_table[column_name][0])

    def get_float_from_table(self, column_name):
        return float(self.data_table[column_name][0])

    def get_int_from_table(self, column_name):
        return int(self.data_table[column_name][0])


class Observations(DataSchema):
    """This is a SQL join of DiaSource and SSSource which contains all of the
    observations of the object.
    """

    def __init__(self, ssObjectId, observations_query, service, sql_filename=None):
        super().__init__(ssObjectId, observations_query, service, sql_filename)

        # This populates each of the variables with a numpy array of the specific column.
        # This should probably be moved to a constructor class method.
        self.mag = self.get_array_from_table("mag")
        self.magErr = self.get_array_from_table("magErr")
        self.mjd = self.get_array_from_table("mjd")
        self.ra = self.get_array_from_table("ra")
        self.dec = self.get_array_from_table("dec")
        self.phaseAngle = self.get_array_from_table("phaseAngle")
        self.topocentricDist = self.get_array_from_table("topocentricDist")
        self.heliocentricDist = self.get_array_from_table("heliocentricDist")

        self.reduced_mag = self.mag - 5 * np.log10(self.topocentricDist * self.heliocentricDist)


class MPCORB(DataSchema):
    """Grabs information from MPCORB."""

    def __init__(self, ssObjectId, observations_query, service, sql_filename=None):
        super().__init__(ssObjectId, observations_query, service, sql_filename)

        self.mpcDesignation = self.get_string_from_table("mpcDesignation")
        self.mpcNumber = self.get_string_from_table("mpcNumber")
        self.mpcH = self.get_float_from_table("mpcH")
        self.mpcG = self.get_float_from_table("mpcH")
        self.epoch = self.get_float_from_table("epoch")
        self.peri = self.get_float_from_table("peri")
        self.node = self.get_float_from_table("node")
        self.incl = self.get_float_from_table("incl")
        self.e = self.get_float_from_table("e")
        self.n = self.get_float_from_table("n")
        self.q = self.get_float_from_table("q")
        self.uncertaintyParameter = self.get_string_from_table("uncertaintyParameter")
        self.flags = self.get_string_from_table("flags")

        # no mean anomaly, no a in MPCORB table


class SSObject(DataSchema):
    """Grabs information from MPCORB."""

    def __init__(self, ssObjectId, observations_query, service, sql_filename=None):
        super().__init__(ssObjectId, observations_query, service, sql_filename)

        self.discoverySubmissionDate = self.get_float_from_table("discoverySubmissionDate")
        self.firstObservationDate = self.get_float_from_table("firstObservationDate")
        self.arc = self.get_float_from_table("arc")
        self.numObs = self.get_int_from_table("numObs")
        self.r_H = self.get_float_from_table("r_H")
        self.r_G12 = self.get_float_from_table("r_G12")
        self.r_Herr = self.get_float_from_table("r_Herr")
        self.r_G12Err = self.get_float_from_table("r_G12err")
        self.r_nData = self.get_int_from_table("r_nData")
        self.maxExtendedness = self.get_float_from_table("maxExtendedness")
        self.minExtendedness = self.get_float_from_table("minExtendedness")
        self.medianExtendedness = self.get_float_from_table("medianExtendedness")
