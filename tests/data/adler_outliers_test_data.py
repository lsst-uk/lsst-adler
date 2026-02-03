from adler.adler_outliers_day import main
from adler.utilities.tests_utilities import get_test_data_path


def generate():
    output_path = get_test_data_path()

    # MPC run with median model
    argv1 = [
        "--process-mjd",
        "60799.5",
        "--input-sql-file",
        "mpc_obs_sbn_testing_database.sqlite",
        "--schema",
        "MPC",
        "--filter-list",
        "g",
        "r",
        "--model-name",
        "median",
        "--data-timespan",
        "7",
        "--n-new-nights",
        "3",
        "--diff-cut",
        "1.5",
        "--std-cut",
        "5.0",
        "--sig-clip-val",
        "3.0",
        "--output-dir",
        f"{output_path}",
        "--logs-dir",
        f"{output_path}",
    ]

    main(argv1)

    # DP0.3 run with phase curve model
    argv2 = [
        "--process-mjd",
        "63335.5",
        "--input-sql-file",
        "testing_database.db",
        "--schema",
        "dp03_catalogs_10yr",
        "--filter-list",
        "r",
        "i",
        "--model-name",
        "HG12_Pen16",
        "--data-timespan",
        "400",
        "--n-new-nights",
        "31",
        "--diff-cut",
        "1.5",
        "--std-cut",
        "5.0",
        "--sig-clip-val",
        "3.0",
        "--output-dir",
        f"{output_path}",
        "--logs-dir",
        f"{output_path}",
    ]

    main(argv2)


if __name__ == "__main__":
    generate()
