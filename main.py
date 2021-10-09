import pandas as pd
import argparse
from pathlib import Path


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--threshold",
        "-t",
        type=float,
        default=1.0,
        help="Proportion of genomes a locus must be present in to be included",
    )

    parser.add_argument("pangenome", type=Path)


def main():
    pass


if __name__ == "__main__":
    main()
