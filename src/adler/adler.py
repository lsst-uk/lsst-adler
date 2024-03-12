import argparse

from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid


def runAdler(args):
    planetoid = AdlerPlanetoid.construct_from_RSP(args.ssoid, args.filter_list, args.date_range)

    planetoid.do_pretend_science()


def main():
    parser = argparse.ArgumentParser(description="Runs Adler for a select planetoid and given user input.")

    parser.add_argument("-s", "--ssoid", help="SSObject ID of planetoid.", type=str, required=True)
    parser.add_argument(
        "-f", "--filters", help="Comma-separated list of filters required.", type=str, default="u,g,r,i,z,y"
    )
    parser.add_argument(
        "-d",
        "--date_range",
        help="Minimum and maximum MJD(TAI) of required observations. Default is to pull all observations.",
        nargs=2,
        type=float,
        default=[60000.0, 67300.0],
    )

    args = parser.parse_args()

    args.filter_list = args.filters.split(",")

    runAdler(args)


if __name__ == "__main__":
    main()
