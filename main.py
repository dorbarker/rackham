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


def threshold_count_calculate(gffs_dir: Path, threshold: float) -> int:

    n_genomes = len(list(gffs_dir.glob("*.gff")))

    return round(threshold * n_genomes)


def loci_over_threshold(gene_families: Path, threshold_count: int) -> set[str]:

    gene_families = pd.read_csv(gene_families, sep=",")

    over_threshold = gene_families["number_genomes"] > threshold_count

    loci = gene_families["gene_family"].loc[over_threshold]

    return set(loci)


def paralogs_get(paralog_list: Path) -> set[str]:

    return set(paralog_list.read_text().splitlines())


def loci_paralogs_remove(loci, paralogs):

    return loci - paralogs


def main():
    pass


if __name__ == "__main__":
    main()
