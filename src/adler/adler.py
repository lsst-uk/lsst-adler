import argparse

from adler.dataclasses.AdlerPlanetoid import AdlerPlanetoid
from adler.science.PhaseCurve import PhaseCurve


def runAdler(args):
    planetoid = AdlerPlanetoid(args.ssoid)

    planetoid.do_pretend_science()

    # now let's do some phase curves!
    pc = PhaseCurve(abs_mag=planetoid.SSObject.r_H, phase_param=0.2, model_name="HG")
    print(pc)


def main():
    parser = argparse.ArgumentParser(description="Runs Adler for a select planetoid and given user input.")

    parser.add_argument("-s", "--ssoid", help="SSObject ID of planetoid.", type=str, required=True)

    # can add arguments to specify a date range etc later
    # alternatively we may start using a config file

    args = parser.parse_args()

    runAdler(args)


if __name__ == "__main__":
    main()
