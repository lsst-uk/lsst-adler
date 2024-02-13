from lsst.rsp import get_tap_service, retrieve_query
from DataSchema import Observations, MPCORB, SSObject

class AdlerPlanetoid:

    def __init__(self, ssObjectId, sql_filename=None):
        
        self.ssObjectId = ssObjectId
        self.sql_filename = sql_filename
        # can also include date ranges at some point

        # can draw from a local SQL database
        if not sql_filename:
            self.service = get_tap_service("ssotap")
        else:
            self.service = None

        # this creates the AdlerPlanetoid.Observations, AdlerPlanetoid.MPCORB and
        # AdlerPlanetoid.SSObject objects.
        self.populate_observations()
        self.populate_MPCORB()
        self.populate_SSObject()


    def populate_observations(self):
        
        observations_sql_query = f"""
            SELECT
                ssObject.ssObjectId, mag, magErr, band, midpointMjdTai as mjd, ra, dec, phaseAngle,
                topocentricDist, heliocentricDist
            FROM
                dp03_catalogs_10yr.ssObject
                JOIN dp03_catalogs_10yr.diaSource ON dp03_catalogs_10yr.ssObject.ssObjectId   = dp03_catalogs_10yr.diaSource.ssObjectId
                JOIN dp03_catalogs_10yr.ssSource  ON dp03_catalogs_10yr.diaSource.diaSourceId = dp03_catalogs_10yr.ssSource.diaSourceId
            WHERE
                ssObject.ssObjectId = {self.ssObjectId} and band='r'
            """
        
        self.Observations = Observations(self.ssObjectId, observations_sql_query, self.service, self.sql_filename)
    
    
    def populate_MPCORB(self):
        
        MPCORB_sql_query = f"""
            SELECT
                ssObjectId, mpcDesignation, mpcNumber, mpcH, mpcG, epoch, peri, node, incl, e, n, q, 
                uncertaintyParameter, flags
            FROM
                dp03_catalogs_10yr.MPCORB
            WHERE
                ssObjectId = {self.ssObjectId}
        """
        
        self.MPCORB = MPCORB(self.ssObjectId, MPCORB_sql_query, self.service, self.sql_filename)
        
    
    def populate_SSObject(self):
        
        SSObject_sql_query = f"""
            SELECT
                discoverySubmissionDate, firstObservationDate, arc, numObs, 
                r_H, r_G12, r_Herr, r_G12err, r_nData,
                maxExtendedness, minExtendedness, medianExtendedness
            FROM
                dp03_catalogs_10yr.SSObject
            WHERE
                ssObjectId = {self.ssObjectId}
        """
        
        self.SSObject = SSObject(self.ssObjectId, SSObject_sql_query, self.service, self.sql_filename)
        