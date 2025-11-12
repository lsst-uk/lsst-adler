from lsst.rsp import get_tap_service
from lsst.rsp.service import get_siav2_service
from lsst.rsp.utils import get_pyvo_auth
from pyvo.dal.adhoc import DatalinkResults, SodaQuery
import lsst.geom as geom
from astropy import units as u
from lsst.afw.image import ExposureF
from lsst.afw.fits import MemFileManager
import os
from astropy.io import fits
from astropy import visualization as aviz
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils.aperture import CircularAperture
from astropy.wcs.utils import skycoord_to_pixel
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from reproject.mosaicking import reproject_and_coadd
from reproject import reproject_interp, reproject_exact


class diaSource_cutout:
    """
    Class to retrieve the cutout for a given difference image source, defined by diaSourceId.
    The additional data required to retrieve the unique image containing the diaSourceId (visit, detector, calib_level) can be passed or will be queried on initialisation.
    Likewise diaSourceId ra and dec required for defining the cutout centre can be passed or will otherwise be queried.
    NB diaSourceId is not required if all other unique image data is provided

    Attributes
    ----------
    radius: Quantity
        Circular radius of requested cutout image (a square image containing the radius is retrieved). Must include astropy angular units, e.g. u.deg
    ra: float
    dec: float
    visit: int
    detector: int
    radius: Quantity
        Angular radius of cutout image (astropy angular units), NB that the query will return the largest square cutout that fits with the radius
    calib_lev: int
        What type of image to retrieve: 1,2,3 - raw, visit, diff
    dataset: str
        Name of the RSP dataset to query for image information and data
    outdir: str
        Where to save the cutout images, will be created if it does not exist. Set to None to save no images
    outfile: str
        Specify a specific filename for the output, otherwise a default will be used
    """

    def __init__(
        self,
        diaSourceId=None,
        ra=None,
        dec=None,
        visit=None,
        detector=None,
        radius=0.01 * u.deg,
        calib_level=3,
        dataset="dp1",
        outdir=".",  # set to None to not save files
        outfile=None,
    ):

        # set up attributes
        self.diaSourceId = diaSourceId
        self.ra = ra
        self.dec = dec
        self.radius = radius
        self.calib_level = calib_level
        self.visit = visit
        self.detector = detector
        self.dataset = dataset
        self.outdir = outdir
        self.outfile = outfile

        # TODO: allow calib_level to be a list to get multiple images in one query?

        if self.outfile is None:
            self.outfile = "cutout-{}_{}.fits".format(self.diaSourceId, self.calib_level)

        # set up RSP services
        self.service_tap = get_tap_service("tap")
        assert self.service_tap is not None
        self.service_sia = get_siav2_service(self.dataset)

        if (self.outdir is not None) & (not os.path.exists(self.outdir)):
            os.makedirs(self.outdir)
            print("Created ", self.outdir)

        assert self.service_sia is not None

    def InitFromDict(self, alert_dict):
        """Set up a new diaSource_cutout object from a dictionary, i.e. an alert packet.
        One could do `diaSource_cutout(**dc_dict)` however any additional columns will cause an error.
        This function removes any dict keys that do not correspond to a diaSource_cutout attribute.

        Parameters
        -----------
        alert_dict : dict
           Dictionary containing the diaSource_cutout parameters, and possibly others

        Returns
        ----------

        dc : object
           The new diaSource_cutout class object

        """

        # clean the input dictionary
        del_keys = []

        # remove keys that aren't diaSource_cutout attributes
        for key, value in alert_dict.items():
            if not hasattr(self, key):
                del_keys.append(key)
        alert_dict = alert_dict.copy()  # make a copy to avoid changing the original dict
        for key in del_keys:
            alert_dict.pop(key, None)

        # initialise a new cutout class
        dc = diaSource_cutout(**alert_dict)
        return dc

    def query_diaSourceId(self):

        # define query to get unique image information for a given diaSourceId
        query = """SELECT ra,dec,visit,detector FROM {}.DiaSource
                WHERE diaSourceId={}
                """.format(
            self.dataset, self.diaSourceId
        )
        print(query)

        # run query
        job = self.service_tap.submit_job(query)
        job.run()
        job.wait(phases=["COMPLETED", "ERROR"])
        print("Job phase is", job.phase)
        if job.phase == "ERROR":
            job.raise_if_error()
        assert job.phase == "COMPLETED"

        # get results
        _df = job.fetch_result().to_table().to_pandas()
        self.ra = _df.iloc[0]["ra"]
        self.dec = _df.iloc[0]["dec"]
        self.visit = int(_df.iloc[0]["visit"])
        self.detector = int(_df.iloc[0]["detector"])

        return _df

    def query_image_url(self):

        # Get the datalink url for the unique image containing the diaSourceId
        query = """SELECT access_url FROM ivoa.ObsCore
        WHERE lsst_visit = {} AND lsst_detector = {} AND calib_level = {}
        """.format(
            self.visit, self.detector, self.calib_level
        )
        print(query)

        job = self.service_tap.submit_job(query)
        job.run()
        job.wait(phases=["COMPLETED", "ERROR"])
        print("Job phase is", job.phase)
        if job.phase == "ERROR":
            job.raise_if_error()
        assert job.phase == "COMPLETED"
        results = job.fetch_result().to_table()
        print(results, len(results))
        assert len(results) == 1

        df_visit = results.to_pandas()
        self.datalink_url = df_visit.iloc[0]["access_url"]

        return self.datalink_url

    def sodaCutout(self):
        """
        Get the image using soda and save it as a fits file.
        """

        # make sure we have a visit and detector id etc, i.e. all unique image data
        if (self.ra is None) or (self.dec is None) or (self.visit is None) or (self.detector is None):
            self.query_diaSourceId()

        # get the image url
        self.query_image_url()

        # find the image on the RSP
        dl_result = DatalinkResults.from_result_url(self.datalink_url, session=get_pyvo_auth())
        sq = SodaQuery.from_resource(
            dl_result, dl_result.get_adhocservice_by_id("cutout-sync"), session=get_pyvo_auth()
        )

        # define the cutout geometry
        # TODO: allow different shapes to be passed to main class?
        spherePoint = geom.SpherePoint(self.ra * geom.degrees, self.dec * geom.degrees)
        sq.circle = (
            spherePoint.getRa().asDegrees() * u.deg,
            spherePoint.getDec().asDegrees() * u.deg,
            self.radius,
        )

        # retrieve the image data for only the area covered by the cutout
        cutout_bytes = sq.execute_stream().read()
        sq.raise_if_error()
        mem = MemFileManager(len(cutout_bytes))
        mem.setData(cutout_bytes, len(cutout_bytes))
        exposure = ExposureF(mem)

        if self.outdir is not None:

            # save the cutout data to file
            self.cutout_file = os.path.join(self.outdir, self.outfile)
            with open(self.cutout_file, "bw") as f:
                f.write(sq.execute_stream().read())
                print("save {}".format(self.cutout_file))

        return exposure

    def plot_exp(self, exposure):
        """
        Plot the lsst exposure object using pyplot.
        """

        # Get image data and header
        # TODO: get WCS working when plotting directly from exp (or use RSP firefly as in tutorials?)
        img = exposure.image

        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0, 0])

        norm = aviz.ImageNormalize(img, interval=aviz.ZScaleInterval())
        s1 = ax1.imshow(img, norm=norm, origin="lower")
        cbar = plt.colorbar(s1)

        return fig

    def plot_img(self, ihdu=1, plot_wcs=False, plot_r=50, wcs_reproj=False, shape_reproj=False):
        """
        Open the saved fits image and plot it using pyplot.

        Parameters
        -----------
        ihdu : int
           Index of fits image containing the image to be plotted
        plot_wcs: Bool
            Flag to use the wcs projection (wcs is read from the fits header, which is only an approximation of the full LSST GBDES WCS model)
        plot_r: float
            Radius of circle to highlight target location (pixels?)
        wcs_reproj: astropy.wcs.wcs.WCS
            Optional, the WCS object to reproject this image to
        shape_reproj: tuple
            Optional, dimensions of the image this one will be reprojected to

        Returns
        ----------
        fig : matplotlib.figure.Figure
           The matplotlib figure object
        """

        # Get image data and header
        hdu = fits.open(self.cutout_file)
        img = hdu[ihdu].data
        hdr = hdu[ihdu].header
        wcs = WCS(hdr)

        if wcs_reproj and shape_reproj:  # TODO: also accept list of hdus and run find_optimal_celestial_wcs
            # reproject the image
            # wcs_reproj, shape_reproj = find_optimal_celestial_wcs(hdu[ihdu])
            data_out, footprint = reproject_and_coadd(
                [hdu[ihdu]],
                wcs_reproj,
                shape_out=shape_reproj,
                reproject_function=reproject_interp,
                # reproject_function=reproject_exact
            )
            # optional, set all zero pixels to nan?
            data_out[data_out == 0] = np.nan

            wcs = wcs_reproj
            img = data_out

        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1)
        if plot_wcs:
            ax1 = plt.subplot(gs[0, 0], projection=wcs)
            ax1.coords.grid(color="white", alpha=0.5, linestyle="solid")
        else:
            ax1 = plt.subplot(gs[0, 0])

        norm = aviz.ImageNormalize(img, interval=aviz.ZScaleInterval())
        s1 = ax1.imshow(img, norm=norm, origin="lower")
        cbar = plt.colorbar(s1)

        # plot target coords
        # TODO: check units here
        c = SkyCoord(self.ra, self.dec, unit=(u.deg, u.deg))
        pos = skycoord_to_pixel(c, wcs)

        # plot the detected sources as apertures with fwhm size
        aperture = CircularAperture(
            positions=np.array(pos),
            # r = df.iloc[i]["trailLength"] / pixscale,
            r=plot_r,
        )
        aperture.plot(color="r")

        plt.title("{}".format(self.visit))

        return fig
