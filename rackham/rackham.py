import numpy as np
import pandas as pd
import argparse
from pathlib import Path
from Bio import SeqIO
import sys
import logging

from rackham import __version__

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
        help="Proportion of genomes a locus must be present in to be included [1.0]",
    )

    parser.add_argument(
        "--length-tolerance",
        "-l",
        type=float,
        default=0.00,
        help="Maximum allowed value for (min_length / max_length) for alleles of a given locus [0.00]",
    )

    parser.add_argument(
        "--max-gaps",
        "-g",
        type=int,
        default=0,
        help="Maximum allowable gaps in a locus [0]",
    )

    parser.add_argument(
        "--mode",
        choices=("mlst", "fsac"),
        default="fsac",
        help="Allele header format to use. mlst = '>aspA_1', fsac = '>1' [fsac]",
    )

    parser.add_argument(
        "--output",
        "-o",
        required=True,
        type=Path,
        help="Output directory containing cgMLST alleles",
    )

    parser.add_argument(
        "-v", "--version", action="version", version=f"{parser.prog} {__version__}"
    )

    parser.add_argument(
        "pangenome", type=Path, help="Pangenome directory created by PIRATE"
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


def is_length_variable(locus: Path, length_tolerance: float, max_gaps: int):

    lengths = []
    gaps = []
    with locus.open("r") as f:
        for record in SeqIO.parse(f, "fasta"):
            lengths.append(len(record.seq))
            gaps.append(str(record.seq).count("-"))

    minimum, maximum = min(lengths), max(lengths)

    is_variable_length = (1.0 - (maximum / minimum)) > length_tolerance
    is_excessively_gapped = max(gaps) > max_gaps

    is_variable = is_variable_length or is_excessively_gapped

    return is_variable


def copy_alleles(
    length_tolerance: float,
    filtered_loci: set[str],
    sequence_directory: Path,
    allele_directory: Path,
    mode: str,
    max_gaps: int,
) -> dict:

    lookups = {}

    allele_directory.mkdir(parents=True, exist_ok=True)

    feature_sequences = sequence_directory.glob("*.fasta")

    for seq in feature_sequences:

        basename = seq.stem.replace(".", "_").split("_")[0]
        logging.info(basename)

        if basename in filtered_loci:

            if not is_length_variable(seq, length_tolerance, max_gaps):

                logging.info(f"Copying {seq}")

                dst = allele_directory.joinpath(seq.name)

                if mode == "mlst":
                    dst = dst.with_suffix(".tfa")

                seq_names = cluster_locus(seq)
                allele_lookup = name_alleles(seq_names)

                lookups[basename] = allele_lookup

                fasta = format_renamed_fasta(seq_names, allele_lookup, basename, mode)

                dst.write_text(fasta)

    return lookups


def convert_call_table(lookups, gene_families):

    selected_loci = lookups.keys()

    loci = pd.read_csv(gene_families, sep="\t")
    selected_rows = loci["gene_family"].isin(selected_loci)
    index = loci["gene_family"].loc[selected_rows]

    selected_columns = list(range(20, len(loci.columns) - 1))

    loci = loci.loc[selected_rows].iloc[:, selected_columns]
    loci.index = index
    loci = loci.transpose()

    for locus, values in loci.iteritems():
        updated_values = values.map(
            lambda call: 0 if pd.isnull(call) else lookups[locus][call]
        )
        loci[locus] = updated_values

    return loci


def cluster_locus(filepath: Path) -> dict[str, list[str]]:

    alleles = {}

    with filepath.open("r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq = str(record.seq)
            name = record.id

            try:
                alleles[seq].append(name)
            except KeyError:
                alleles[seq] = [name]

    return alleles


def name_alleles(seq_names: dict[str, list[str]]) -> dict[str, int]:

    allele_number_lookup = {}

    for i, names in enumerate(seq_names.values(), 1):
        for name in names:
            allele_number_lookup[name] = i

    return allele_number_lookup


def format_renamed_fasta(
    seq_names: dict[str, list[str]],
    allele_number_lookup: dict[str, int],
    basename: str,
    mode: str,
):

    leader = f"{basename}_" if mode == "mlst" else ""
    multifasta = []
    for sequence, names in seq_names.items():
        allele_num = allele_number_lookup[names[0]]

        multifasta.append((allele_num, sequence))

    formatted_multifasta = "\n".join(
        [f">{leader}{i}\n{seq}" for i, seq in sorted(multifasta)]
    )

    return formatted_multifasta


def create_st_table(calls_table: pd.DataFrame) -> pd.DataFrame:

    column_order = ["ST"] + [str(col) for col in calls_table.columns]

    sequence_types = calls_table.drop_duplicates()
    sequence_types["ST"] = np.arange(1, len(sequence_types) + 1)

    sequence_types = sequence_types.reset_index(drop=True)[column_order]

    return sequence_types


def filtered_loci_get(
    carriage_threshold: float,
    length_tolerance: float,
    gffs_dir: Path,
    paralog_clusters: Path,
    gene_families: Path,
    feature_sequences_directory: Path,
    alleles_directory: Path,
    mode: str,
    max_gaps: int,
) -> pd.DataFrame:

    threshold_count = threshold_count_calculate(gffs_dir, carriage_threshold)

    high_carriage = loci_over_threshold(gene_families, threshold_count)

    paralogs = paralogs_get(paralog_clusters)

    high_carriage_orthologs = loci_paralogs_remove(high_carriage, paralogs)

    lookups = copy_alleles(
        length_tolerance=length_tolerance,
        filtered_loci=high_carriage_orthologs,
        sequence_directory=feature_sequences_directory,
        allele_directory=alleles_directory,
        mode=mode,
        max_gaps=max_gaps,
    )

    calls_table = convert_call_table(lookups, gene_families)

    return calls_table


def main():

    args = arguments()

    calls_table = filtered_loci_get(
        carriage_threshold=args.carriage_threshold,
        length_tolerance=args.length_tolerance,
        gffs_dir=args.pangenome / "modified_gffs/",
        paralog_clusters=args.pangenome / "paralog_clusters.tab",
        gene_families=args.pangenome / "PIRATE.gene_families.tsv",
        feature_sequences_directory=args.pangenome / "feature_sequences/",
        alleles_directory=args.output / "alleles",
        mode=args.mode,
        max_gaps=args.max_gaps,
    )

    calls_table.to_csv(args.output / "calls.tsv", sep="\t", index_label="genome")

    sequence_types = create_st_table(calls_table)
    sequence_types.to_csv(args.output / "sts.txt", sep="\t", index=False)


class IncompleteCopyError(Exception):
    pass


if __name__ == "__main__":
    main()
