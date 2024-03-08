import argparse

from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid


def runAdler(args):
    planetoid = AdlerPlanetoid.construct_from_RSP(args.ssoid, args.filter_list)

    planetoid.do_pretend_science()


def main():
    parser = argparse.ArgumentParser(description="Runs Adler for a select planetoid and given user input.")

    parser.add_argument("-s", "--ssoid", help="SSObject ID of planetoid.", type=str, required=True)
    parser.add_argument(
        "-f", "--filters", help="Comma-separated list of filters required.", type=str, default="u,g,r,i,z,y"
    )

    # can add arguments to specify a date range etc later

    args = parser.parse_args()

    args.filter_list = args.filters.split(",")

    runAdler(args)


if __name__ == "__main__":
    main()
