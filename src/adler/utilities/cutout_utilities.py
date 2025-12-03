from lsst.rsp import get_tap_service
from lsst.rsp.service import get_siav2_service
from lsst.rsp.utils import get_pyvo_auth
from pyvo.dal.adhoc import DatalinkResults, SodaQuery
import lsst.geom as geom
from astropy import units as u
from lsst.afw.image import ExposureF
from lsst.afw.fits import MemFileManager
import lsst.afw.display as afwDisplay

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

CALIB_DICT = {1: "raw", 2: "visit", 3: "diff"}


class Cutout:
    """
    Class to retrieve the cutout for a given source, defined by either diaSourceId or SourceId.
    The additional data required to retrieve the unique image containing the source (visit, detector, calib_level) can be passed or will be queried on initialisation.
    Likewise the source's ra and dec required for defining the cutout centre can be passed or will otherwise be queried.
    NB a source Id is not required if all other unique image data is provided.

    Attributes
    ----------
    Id: int
        Identifier value for the source
    IdTable: str
        The table from which the source is queried, e.g. DiaSource or Source. We assume the Id column name is IdTable+"Id".
    ra: float
    dec: float
    visit: int
    detector: int
    radius: Quantity
        Angular radius of cutout image (astropy angular units), NB that the query will return the largest square cutout that fits with the radius
    calib_lev: list(int)
        List of what type of images to retrieve: 1,2,3 - raw, visit, diff. NB, this arg will accept a single int
    dataset: str
        Name of the RSP dataset to query for image information and data
    outdir: str
        Where to save the cutout images, will be created if it does not exist. Set to None to save no images
    outfile: str
        Specify a specific filename for the output, otherwise a default will be used
    """

    def __init__(
        self,
        Id=None,
        IdTable="DiaSource",
        ra=None,
        dec=None,
        visit=None,
        detector=None,
        radius=0.01 * u.deg,
        calib_level=[3],
        dataset="dp1",
        outdir=".",  # set to None to not save files
        outfile=None,
    ):

        # set up attributes
        self.Id = Id
        self.IdTable = IdTable
        self.ra = ra
        self.dec = dec
        self.visit = visit
        self.detector = detector
        self.radius = radius
        self.calib_level = calib_level
        self.dataset = dataset
        self.outdir = outdir
        self.outfile = outfile

        # define the Id column name
        self.IdCol = self.IdTable + "Id"

        # Enusure calib_level is a list
        if type(self.calib_level) is int:
            self.calib_level = [self.calib_level]

        if self.outfile is None:
            if self.Id is None:
                self.outfile = "cutout-{}-{}_ra{}_dec{}".format(self.visit, self.detector, self.ra, self.dec)
            else:
                self.outfile = "cutout-{}".format(self.Id)

        # set up RSP services
        self.service_tap = get_tap_service("tap")
        assert self.service_tap is not None
        self.service_sia = get_siav2_service(self.dataset)

        if (self.outdir is not None) & (not os.path.exists(self.outdir)):
            os.makedirs(self.outdir)
            print("Created ", self.outdir)

        assert self.service_sia is not None

    def InitFromDict(self, alert_dict):
        """Set up a new Cutout object from a dictionary, i.e. an alert packet.
        One could do `Cutout(**co_dict)` however any additional columns will cause an error.
        This function removes any dict keys that do not correspond to a Cutout attribute.

        Parameters
        -----------
        alert_dict : dict
           Dictionary containing the Cutout parameters, and possibly others

        Returns
        ----------

        co : object
           The new Cutout class object

        """

        # get any previously set values
        co_dict = self.__dict__

        # remove non-arg keys # TODO: repeat for alert_dict?
        del_keys = []
        co_args = self.__init__.__code__.co_varnames
        for key, value in co_dict.items():
            if key not in co_args:
                del_keys.append(key)
        for key in del_keys:
            co_dict.pop(key, None)

        # Force overwrite of outfile, unless it is passed in alert_dict?
        co_dict["outfile"] = None

        # Overwrite initial Cutout attributes if they are in alert_dict
        for key, value in alert_dict.items():
            if key in co_dict:
                co_dict[key] = value

        # initialise a new cutout class
        co = Cutout(**co_dict)

        return co

    def query_from_Id(self):

        # define query to get unique image information for a given source Id
        query = """SELECT ra,dec,visit,detector FROM {}.{}
                WHERE {}={}
                """.format(
            self.dataset, self.IdTable, self.IdCol, self.Id
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

        # Get the datalink url for the unique image containing the Id
        if len(self.calib_level) == 1:
            query = """SELECT access_url, calib_level FROM ivoa.ObsCore
            WHERE lsst_visit = {} AND lsst_detector = {} AND calib_level = {}
            """.format(
                self.visit, self.detector, self.calib_level[0]
            )
        else:
            query = """SELECT access_url, calib_level FROM ivoa.ObsCore
            WHERE lsst_visit = {} AND lsst_detector = {} AND calib_level IN {}
            """.format(
                self.visit, self.detector, tuple(self.calib_level)
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
        # assert len(results) == 1
        assert len(results) == len(self.calib_level)

        df_visit = results.to_pandas()
        self.datalink_url = df_visit

        return self.datalink_url

    def sodaCutout(self, small_fits=True):
        """
        Get the image using soda and save it as a fits file.

        Parameters
        -----------
        small_fits : Bool
           Flag to save only the image component to fits file (set False to save full exposure to fits file)

        Returns
        ----------

        exposure : ExposureF
           The new Cutout class object
        """

        # make sure we have a visit and detector id etc, i.e. all unique image data
        if (self.ra is None) or (self.dec is None) or (self.visit is None) or (self.detector is None):
            self.query_from_Id()

        # get the image url
        self.query_image_url()

        for i in range(len(self.datalink_url)):  # TODO: Parallelise here?

            url = self.datalink_url.iloc[i]["access_url"]
            cal_lev = self.datalink_url.iloc[i]["calib_level"]
            print(url, cal_lev)

            # find the image on the RSP
            dl_result = DatalinkResults.from_result_url(url, session=get_pyvo_auth())
            print(f"Datalink status: {dl_result.status}. Datalink service url: {url}")
            sq = SodaQuery.from_resource(
                dl_result,
                dl_result.get_adhocservice_by_id("cutout-sync"),
                session=get_pyvo_auth(),  # TODO: query fails for cal_lev = 1, raw images?
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

                # save the cutout data to file - add the calib_level to the filename
                cutout_file = os.path.join(self.outdir, self.outfile) + "_" + str(cal_lev) + ".fits"
                setattr(
                    self, "cutout_file_" + CALIB_DICT[cal_lev], cutout_file
                )  # store the filename as an attribute

                if small_fits:
                    # Save the minimum information to keep file size down, see afwDisplay.writeFitsImage in https://dp1.lsst.io/tutorials/notebook/103/notebook-103-4.html
                    afwDisplay.writeFitsImage(
                        cutout_file,  # see, https://dp1.lsst.io/tutorials/notebook/103/notebook-103-4.html
                        # exposure.maskedImage, # slightly smaller than full exposure
                        exposure.image,  # small as possible, no masks
                        wcs=exposure.wcs,
                    )
                else:
                    with open(cutout_file, "bw") as f:
                        f.write(sq.execute_stream().read())
                        print("save {}".format(cutout_file))

            # store exposure as a class attribute
            setattr(self, "exp_" + CALIB_DICT[cal_lev], exposure)

        return exposure

    def plot_img(
        self, cal_lev, plot_fits=True, ihdu=0, plot_wcs=False, plot_r=50, wcs_reproj=False, shape_reproj=False
    ):
        """
        Open the saved fits image and plot it using pyplot.
        Use the arguments wcs_reproj and shape_reproj to reproject the image to a specific wcs and shape (e.g. get North up and East pointing left). NB beware of flux conservation when reprojecting
        Find the best wcs and shape from a list of hdus using something like:
        `wcs_reproj, shape_reproj = find_optimal_celestial_wcs(hdu_list, hdu_in = "IMAGE")`

        Parameters
        -----------
        cal_lev: int
            Calibration level to be plotted
        fits: Bool
            Flag to plot from the cutout fits file
        ihdu : int or str
           Index (or hdu label) of fits image containing the image to be plotted
        plot_wcs: Bool
            Flag to use the wcs projection (wcs is read from the fits header, which is only an approximation of the full LSST GBDES WCS model)
        plot_r: float
            Radius of circle to highlight target location (pixels?)
        wcs_reproj: astropy.wcs.wcs.WCS
            Optional, a specific WCS object to reproject this image to
        shape_reproj: tuple
            Optional, specific dimensions of the image this one will be reprojected to

        Returns
        ----------
        fig : matplotlib.figure.Figure
           The matplotlib figure object
        """

        # get the ExposureF object
        exp = getattr(self, "exp_" + CALIB_DICT[cal_lev])

        if plot_fits:
            # Get image data and header
            cutout_file = getattr(self, "cutout_file_" + CALIB_DICT[cal_lev])
            hdu = fits.open(cutout_file)
            img = hdu[ihdu].data
            hdr = hdu[ihdu].header
            wcs = WCS(hdr)
        else:
            # get the fits equivalent data
            img = exp.getImage().getArray()  # image data
            # hdr = Header(exp.getMetadata().toDict()) # fits (primary) header info (without WCS) - https://community.lsst.org/t/obtain-only-image-metadata-with-butler-get/5922
            wcs = WCS(
                exp.getWcs().getFitsMetadata()
            )  # ExposureF WCS converted to fits WCS. NB Will require a shift to be applied to calculated pixel coords, getXY0(), https://community.lsst.org/t/visualizing-images-in-sky-coordinates-using-wcs-in-a-notebook/4210/8

        if wcs_reproj and shape_reproj:
            if not plot_fits:
                # TODO: log an error here
                print(
                    "ERROR: Reprojection requires plotting from the fits file WCS representation (set plot_fits=True)"
                )
            else:
                # reproject the image, using shape and wcs determined from list of hdus and find_optimal_celestial_wcs
                data_out, footprint = reproject_and_coadd(
                    [hdu[ihdu]],
                    wcs_reproj,
                    shape_out=shape_reproj,
                    reproject_function=reproject_interp,
                    # reproject_function=reproject_exact
                )
                # Set all zero pixels to nan to make plot nicer (TODO: make optional?)
                data_out[data_out == 0] = np.nan

                wcs = wcs_reproj
                img = data_out

        # create the matplotlib figure
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1)
        if plot_wcs:  # plot the image using the Fits approximate WCS
            ax1 = plt.subplot(gs[0, 0], projection=wcs)
            ax1.coords.grid(color="white", alpha=0.5, linestyle="solid")
        else:  # just plot in pixel scale
            ax1 = plt.subplot(gs[0, 0])

        # Plot the image with norm and scaling
        norm = aviz.ImageNormalize(img, interval=aviz.ZScaleInterval())
        s1 = ax1.imshow(img, norm=norm, origin="lower")
        cbar = plt.colorbar(s1)

        # plot target coords using the INEXACT fits approximation WCS
        c = SkyCoord(self.ra, self.dec, unit=(u.deg, u.deg))
        pos = skycoord_to_pixel(c, wcs)
        if not plot_fits:
            print("correct pixel pos for ExposureF WCS to fits")
            pos = pos - np.array(exp.image.getXY0())
        print("Fits WCS pos = {}".format(pos))
        # plot the detected sources as apertures (scale with fwhm/trail size?)
        aperture = CircularAperture(
            positions=np.array(pos),
            # r = df.iloc[i]["trailLength"] / pixscale,
            r=plot_r,
        )
        aperture.plot(color="r")

        # Plot the exact position with the ExposureF WCS
        if wcs_reproj and shape_reproj:
            # TODO: log a warning here
            print(
                "WARNING: Reprojection requires the fits file WCS representation and does not work with ExposureF WCS"
            )
        coord = geom.SpherePoint(self.ra * geom.degrees, self.dec * geom.degrees)
        pos = exp.wcs.skyToPixel(coord)  # pos.x, pos.y
        pos = (
            pos - exp.image.getXY0()
        )  # Get the shifted position to correct for the WCS - https://community.lsst.org/t/how-to-use-wcss-in-dp1-and-commissioning-processing/10769
        print("ExposureF WCS pos = {}".format(pos))
        aperture = CircularAperture(
            positions=np.array(pos),
            r=plot_r,
        )
        aperture.plot(color="w")

        plt.title("{}\n{} {}".format(self.Id, self.visit, CALIB_DICT[cal_lev]))

        return fig
