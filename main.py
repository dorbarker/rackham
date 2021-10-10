import pandas as pd
import argparse
from pathlib import Path
from shutil import copyfile


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


def copy_alleles(
    filtered_loci: set[str], sequence_directory: Path, allele_directory: Path
) -> None:

    allele_directory.mkdir(parents=True, exist_ok=True)

    feature_sequences = sequence_directory.glob("*.fasta")

    for seq in feature_sequences:

        if seq.replace(".", "_").split("_")[0] in filtered_loci:

            dst = allele_directory.joinpath(seq.name)

            copyfile(seq, dst)

    # check that all is well
    if len(list(allele_directory.glob("*.fasta"))) != len(filtered_loci):
        raise IncompleteCopyError


def main():
    pass


class IncompleteCopyError(Exception):
    pass


if __name__ == "__main__":
    main()
