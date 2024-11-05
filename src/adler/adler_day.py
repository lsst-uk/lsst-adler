from adler.objectdata.AdlerPlanetoid import AdlerPlanetoid
from adler.science.PhaseCurve import PhaseCurve
import os
import pandas as pd
import astropy.units as u
import numpy as np


def run_adler_day(ssobject_list):
    """Runs the Adler day process for a list of SSObjectIds.

    Parameters
    -----------
    ssobject_list : list of str
        List of strings of the SSObjectIDs we want to investigate.

    """

    phase_model = "HG12_Pen16"

    for ssoid in ssobject_list:

        # date range is start of survey to "today's" date, these are dummy dates
        # this should also pull in the most recent alert packet, which should be in Cassandra by now?
        planetoid = AdlerPlanetoid.construct_from_cassandra(ssoid, date_range=[60290.0, 61902])

        try:
            dirname = os.path.dirname(__file__)
            database_path = os.path.join(dirname, "../../tests/data/adler-day-testing.db")
            planetoid.attach_previous_adler_data(database_path)
            previousAdlerData = True  # this could be a flag on the AdlerPlanetoid object
        except ValueError:
            print("No data found for this object in AdlerDatabase.")
            previousAdlerData = False

        for filt in planetoid.filter_list:
            obs = planetoid.observations_in_filter(filt)
            df_obs = pd.DataFrame(obs.__dict__)

            ######## DO WE FIT? ########

            # here we assume that there have been enough new points to justify refitting the phase curve
            # ie: that triage has already been done in adler-night
            # whether we refit depends on some complex considerations (n_obs, phase angle range, minimum phase angle)
            # however: could do a basic check: five new data points?

            ######## FIT NEW PHASE CURVE ########

            if previousAdlerData:
                params = planetoid.PreviousAdlerData.get_phase_parameters_in_filter(filt, phase_model)
                H = params.H
                G = params.phase_parameter_1
            else:
                sso = planetoid.SSObject_in_filter(filt)
                H = sso.H
                G = sso.G12

            if np.isnan(H) or np.isnan(G):
                print(
                    "No H and/or phase parameter value found for this object and filter in SSObject or Adler databases."
                )
                print("Using guess values. Fit may take longer.")
                H = 15.0
                G = 0.5

            pc = PhaseCurve(
                H=H * u.mag,
                phase_parameter_1=G,
                model_name=phase_model,
            )

            pc_fit = pc.FitModel(
                np.array(df_obs["phaseAngle"]) * u.deg,
                np.array(df_obs["reduced_mag"]) * u.mag,
                np.array(df_obs["magErr"]) * u.mag,
            )

            pc_fit = pc.InitModelSbpy(pc_fit)
            pc_fit.__dict__["H"] = pc_fit.__dict__["H"].value

            planetoid.AdlerData.populate_phase_parameters(filt, **pc_fit.__dict__)

            ######## OUTLIER DETECTION: ONE METHOD ########

            # can do phase curve fit including all but last point for outlier detection
            # use H and G from alert if first time, otherwise pull from AdlerData

            # check for that phase curve model what the predicted values are, check alert point against those within threshold
            # can later compare to LSST-predicted mag

            # of course this only checks the most recent point - what if the outlier was ten alerts ago but we didn't know?
            # what if it was an alert we triaged out?

        # write new row to Adler database
        # planetoid.AdlerData.write_row_to_database("adler-day-testing.db")


if __name__ == "__main__":

    # an assumption that these are the ssObjectIds of interest
    # in the future, there will be a function which reads from the Kafka topic of triaged alerts from adler-night
    # and returns this list.
    # why are some of them negative? Ken doesn't know
    ssobject_list = [
        "685835657114313227",
        "-6411717652251119994",
        "3881603647193383121",
        "6711448832195578320",
        "2959207197656523004",
        "-8724852777827460529",
        "-7394109909111663666",
        "2154070653759973672",
        "-7805034407650925178",
        "-4490542680786866656",
    ]

    run_adler_day(ssobject_list)
