import pandas as pd
import argparse
from pathlib import Path
from shutil import copyfile
from Bio import SeqIO
import sys
import logging

logging.basicConfig(
    format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S", level=logging.INFO
)


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--carriage-threshold",
        "-t",
        type=float,
        default=1.0,
        help="Proportion of genomes a locus must be present in to be included",
    )

    parser.add_argument(
        "--length-tolerance",
        "-l",
        type=float,
        default=0.01,
        help="Maximum allowed value for (min_length / max_length) for alleles of a given locus",
    )

    parser.add_argument(
        "pangenome", type=Path, help="Pangenome directory created by PIRATE"
    )

    parser.add_argument(
        "--output",
        "-o",
        required=True,
        type=Path,
        help="Output directory containing cgMLST alleles",
    )

    args = parser.parse_args()

    if not args.pangenome.is_dir():
        msg = f"{args.pangenome} is not a directory"
        print(msg, file=sys.stderr)
        sys.exit(1)

    return args


def threshold_count_calculate(gffs_dir: Path, threshold: float) -> int:

    n_genomes = len(list(gffs_dir.glob("*.gff")))
    logging.info(f"n_genomes: {n_genomes}")

    cutoff_count = round(threshold * n_genomes)

    logging.info(f"cutoff_count: {cutoff_count}")

    return cutoff_count


def loci_over_threshold(gene_families: Path, threshold_count: int) -> set[str]:

    gene_families = pd.read_csv(gene_families, sep="\t")

    over_threshold = gene_families["number_genomes"] > threshold_count

    loci = gene_families["gene_family"].loc[over_threshold]

    logging.info(f"loci over threshold: {len(loci)}")

    return set(loci)


def paralogs_get(paralog_list: Path) -> set[str]:

    return set(paralog_list.read_text().splitlines())


def loci_paralogs_remove(loci, paralogs):

    return loci - paralogs


def is_length_variable(locus: Path, length_tolerance: float):

    with locus.open("r") as f:
        lengths = [len(rec.seq) for rec in SeqIO.parse(f, "fasta")]

    minimum, maximum = min(lengths), max(lengths)

    is_variable = (1.0 - (maximum / minimum)) > length_tolerance

    logging.info(f"{minimum} {maximum} {is_variable}")

    return is_variable

def copy_alleles(
    length_tolerance: float,
    filtered_loci: set[str],
    sequence_directory: Path,
    allele_directory: Path,
) -> None:

    allele_directory.mkdir(parents=True, exist_ok=True)

    feature_sequences = sequence_directory.glob("*.fasta")

    for seq in feature_sequences:

        bn = seq.stem.replace(".", "_").split("_")[0]
        logging.info(bn)

        if bn in filtered_loci:

            if not is_length_variable(seq, length_tolerance):

                logging.info(f"Copying {seq}")

                dst = allele_directory.joinpath(seq.name)

                copyfile(seq, dst)


def filtered_loci_get(
    carriage_threshold: float,
    length_tolerance: float,
    gffs_dir: Path,
    paralog_clusters: Path,
    gene_families: Path,
    feature_sequences_directory: Path,
    alleles_directory: Path,
) -> None:

    threshold_count = threshold_count_calculate(gffs_dir, carriage_threshold)

    high_carriage = loci_over_threshold(gene_families, threshold_count)

    paralogs = paralogs_get(paralog_clusters)

    high_carriage_orthologs = loci_paralogs_remove(high_carriage, paralogs)

    copy_alleles(
        length_tolerance=length_tolerance,
        filtered_loci=high_carriage_orthologs,
        sequence_directory=feature_sequences_directory,
        allele_directory=alleles_directory,
    )


def main():

    args = arguments()

    filtered_loci_get(
        carriage_threshold=args.carriage_threshold,
        length_tolerance=args.length_tolerance,
        gffs_dir=args.pangenome / "modified_gffs/",
        paralog_clusters=args.pangenome / "paralog_clusters.tab",
        gene_families=args.pangenome / "PIRATE.gene_families.tsv",
        feature_sequences_directory=args.pangenome / "feature_sequences/",
        alleles_directory=args.output,
    )


class IncompleteCopyError(Exception):
    pass


if __name__ == "__main__":
    main()
