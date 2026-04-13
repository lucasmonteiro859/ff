#!/opt/fastfilter/venv/bin/python

"""
fastfilter.py - version 1.0a1
Author: Gil Poiares-Oliveira <gpo@ciencias.ulisboa.pt>
PI: Margarida Gama-Carvalho <mhcarvalho@ciencias.ulisboa.pt>

RNA Systems Biology Lab
BioISI - Biosystems and Integrative Sciences Institute
Department of Chemistry and Biochemistry
Faculty of Sciences, University of Lisbon

(C) 2022-2023
"""

import argparse
import csv
import multiprocessing
import sys
import time
from pathlib import Path

import inquirer
import pandas as pd
from Bio import SeqIO, SeqRecord
from Bio.SeqIO.QualityIO import FastqPhredIterator

from tqdm import tqdm

# Defaults
# Can you edit these values directly in code to fit your needs? Yes!
# Should you? No! Supply the values as arguments on runtime, please!
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
PROJECTS_DIR_DEFAULT = Path("/data/working_directory/projects")

# These variables are meant to be overwritten if the user specifies the
# parameters manually at runtime
min_seq_len = MIN_LENGTH_DEFAULT
homopolymer_coeff = HOMOPOLYMER_COEFF_DEFAULT
min_score = MIN_SCORE_DEFAULT
single_end = False

# For debugging, currently
dryrun = False


# Parse arguments
def parse_arguments() -> None:
    """Parse arguments supplied on runtime and create global variables to
    store those values. Any unsupplied parameters will be set to the default
    values retrived from MIN_LENGTH_DEFAULT and HOMOPOLYMER_COEFF_DEFAULT
    """
    parser = argparse.ArgumentParser(description="Filter FASTQ files.")
    parser.add_argument(
        "-l",
        "--minlen",
        metavar="50",
        type=int,
        default=MIN_LENGTH_DEFAULT,
        help="sequence length threshold",
    )
    parser.add_argument(
        "-p",
        "--homopolymerlen",
        metavar="3",
        type=int,
        default=HOMOPOLYMER_COEFF_DEFAULT,
        help="homopolymer length",
    )
    parser.add_argument(
        "-d",
        "--dryrun",
        action=argparse.BooleanOptionalAction,
        help="Do a dry run for debbuging.",
    )
    parser.add_argument("-s", "--min-score", type=int, metavar="30")
    parser.add_argument(
        "-i",
        "--sequences-dir",
        type=str,
        metavar="/data/working_directory/project/cutadapt",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        metavar="/data/working_directory/project/fastfilter",
    )
    parser.add_argument("-j", "--cpus", type=int, default=1, help="Number of parallel CPUs")
    args = parser.parse_args()
    global min_seq_len
    global homopolymer_coeff
    global dryrun
    global min_score
    global seq_dir
    global output_dir
    global num_cpus
    min_seq_len = args.minlen
    homopolymer_coeff = args.homopolymerlen
    dryrun = args.dryrun
    min_score = args.min_score
    seq_dir = Path(args.sequences_dir) if args.sequences_dir else None
    output_dir = Path(args.output_dir) if args.output_dir else None
    num_cpus = args.cpus


def find_homopolymers(*seq) -> dict:
    """Find homopolymers in a sequence. If no sequence is given, all
    homopolymer counts are returned as 0.

    Keyword arguments:
    seq -- The nucleotide sequence of the readout
    """
    global homopolymer_coeff
    if len(seq) == 0:
        """If no sequence is provided, dict gets initialized
        will 0 in all counts"""
        a_count = 0
        t_count = 0
        g_count = 0
        c_count = 0
    else:
        seq = seq[0]
        a_count = int(seq.count(homopolymer_coeff * "A") > 0)
        t_count = int(seq.count(homopolymer_coeff * "T") > 0)
        g_count = int(seq.count(homopolymer_coeff * "G") > 0)
        c_count = int(seq.count(homopolymer_coeff * "C") > 0)

    return {
        homopolymer_coeff * "A": a_count,
        homopolymer_coeff * "T": t_count,
        homopolymer_coeff * "G": g_count,
        homopolymer_coeff * "C": c_count,
    }


def export_length_frequencies(
    seq_lens: list[int], filename_no_extension: str, output_dir: Path
) -> None:
    """Export ordered length frequencies to CSV file"""
    lenData = (
        pd.Series(seq_lens)
        .value_counts()
        .reset_index()
        .sort_values("index")
        .reset_index(drop=True)
    )
    lenData.columns = ["Length", "Frequency"]
    output_path = output_dir / (filename_no_extension + "_FREQ.csv")
    lenData.to_csv(output_path, index=False)


def analyze_sequence(record: SeqRecord) -> dict:
    """Main sequence analysis function. Scans one readout for
    homopolymers (by calling find_polymers())
    """
    qual_values = record.letter_annotations["phred_quality"]
    try:
        qual_score = sum(qual_values) / len(qual_values)
    except ZeroDivisionError:
        qual_score = 0
    n_count = record.seq.count("N")
    dot_count = record.seq.count(".")

    global min_seq_len
    seq_len = len(record.seq)
    should_include = True
    too_short = False
    found_homopolymer = False
    low_score = False
    # Homopolymers dict with counts each of the different types
    # of homopolymers
    homopolymers = find_homopolymers(record.seq)
    # homopolymers_exist set to True if any homopolymer is found
    # in the sequence
    homopolymers_exist = any(x > 0 for x in homopolymers.values())
    if seq_len < min_seq_len:
        should_include = False
        too_short = True
    if n_count > 0:
        should_include = False
    if dot_count > 0:
        should_include = False
    if homopolymers_exist:
        should_include = False
        found_homopolymer = True
    if qual_score < min_score:
        should_include = False
        low_score = True
    seq_qualities = {
        "meets_criteria": should_include,
        "length": seq_len,
        "qual_score": qual_score,
        "n_count": n_count,
        "dot_count": dot_count,
        "homopolymers": homopolymers,
        "too_short": too_short,
        "should_include": should_include,
        "found_homopolymer": found_homopolymer,
        "low_score": low_score,
    }

    return seq_qualities


def parse_file_single_end(r1_filename: Path) -> dict:
    filename1_no_extension = r1_filename.stem
    output_dir = r1_filename.parent.parent / "fastfilter"
    table1_path = Path(output_dir / (filename1_no_extension + "_TABLE.csv"))

    # Initialize empty variables
    total_seqs = 0
    good_reads_r1 = 0
    homopolymers = find_homopolymers()
    seq_lens1 = []
    exclusion_reasons_r1 = {
        "too_short": 0,
        "found_dot": 0,
        "found_n": 0,
        "found_homopolymer": 0,
        "low_score": 0,
    }

    with open(table1_path, "x", newline="") as csv1_file:
        header = [
            "id",
            "sequence",
            "score",
            "size",
            "n_count",
            "dot_count",
            "homopolymers",
        ]
        table1_writer = csv.writer(csv1_file, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        table1_writer.writerow(header)

        f_r1 = open(r1_filename, "r")

        filtered1 = open(output_dir / (filename1_no_extension + "_FILTERED.fastq"), "x")

        r1_iterator = FastqPhredIterator(f_r1)

        good_sequences = 0

        for record1 in r1_iterator:
            total_seqs += 1  # Counter of total sequences in input files

            seq1 = analyze_sequence(record1)

            seq_lens1.append(seq1["length"])

            homopolymers = {
                x: homopolymers[x] + seq1["homopolymers"][x]
                for x in set(homopolymers).intersection(seq1["homopolymers"])
            }

            good_reads_r1 += seq1["meets_criteria"]

            if seq1["meets_criteria"]:
                # If both sequences meet criteria, write to output FASTQs
                SeqIO.write(record1, filtered1, "fastq")
                good_sequences += 1

            else:
                # Otherwise, quantify all the exclusion_reasons
                # The following lines of code take advantage of the fact that
                # Python considers bool(x) == True, if x != 0

                # R1
                exclusion_reasons_r1["too_short"] += bool(seq1["too_short"])
                exclusion_reasons_r1["found_dot"] += bool(seq1["dot_count"])
                exclusion_reasons_r1["found_n"] += bool(seq1["n_count"])
                exclusion_reasons_r1["found_homopolymer"] += bool(
                    seq1["found_homopolymer"]
                )
                exclusion_reasons_r1["low_score"] += seq1["low_score"]

            # Add sequence to output list
            if not dryrun:
                # Add row to the corresponding "TABLE" file
                table1_writer.writerow([record1.id, record1.seq] + list(seq1.items()))

        # Export file with sequence length frequencies
        export_length_frequencies(seq_lens1, filename1_no_extension, output_dir)

        # Close opened files
        f_r1.close()
        filtered1.close()

        # Return dict with the paired end's cleaning stats to be compiled in
        # fastfilter_overview.csv file
        return {
            "file": filename1_no_extension,
            "total_seqs": total_seqs,
            "good_sequences_count": good_sequences,
            "homopolymers": homopolymers,
            "exclusion_reasons_r1": exclusion_reasons_r1,
            "good_reads_r1": good_reads_r1,
        }


def parse_file(r1_filename: Path, r2_filename: Path, output_dir: Path, position: int = 0) -> dict:
    """Run the parsing workflow for a pair of R1/R2 FASTQ files.
    This function should be called for each pair of paired-end readouts.
    """
    filename1_no_extension = r1_filename.stem
    filename2_no_extension = r2_filename.stem
    # output_dir = r1_filename.parent.parent / "fastfilter"
    table1_path = Path(output_dir / (filename1_no_extension + "_TABLE.csv"))
    table2_path = Path(output_dir / (filename2_no_extension + "_TABLE.csv"))

    # Initialize empty variables
    total_seqs = 0
    good_reads_r1 = 0
    good_reads_r2 = 0
    homopolymers = find_homopolymers()
    seq_lens1 = []
    seq_lens2 = []
    exclusion_reasons_r1 = {
        "too_short": 0,
        "found_dot": 0,
        "found_n": 0,
        "found_homopolymer": 0,
        "low_score": 0,
    }

    exclusion_reasons_r2 = {
        "too_short": 0,
        "found_dot": 0,
        "found_n": 0,
        "found_homopolymer": 0,
        "low_score": 0,
    }

    with open(table1_path, "x", newline="") as csv1_file, open(
        table2_path, "x", newline=""
    ) as csv2_file:
        header = [
            "id",
            "sequence",
            "score",
            "size",
            "n_count",
            "dot_count",
            "homopolymers",
        ]
        table1_writer = csv.writer(csv1_file, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        table2_writer = csv.writer(csv2_file, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        # table1_writer.writerow(header)
        # table2_writer.writerow(header)

        f_r1 = open(r1_filename, "r")
        f_r2 = open(r2_filename, "r")

        filtered1 = open(output_dir / (filename1_no_extension + "_FILTERED.fastq"), "x")
        filtered2 = open(output_dir / (filename2_no_extension + "_FILTERED.fastq"), "x")

        r1_iterator = FastqPhredIterator(f_r1)
        r2_iterator = FastqPhredIterator(f_r2)

        good_sequences = 0

        # --- NEW: streaming sync of R1/R2 by read ID ---
        r1_next = next(r1_iterator, None)
        r2_next = next(r2_iterator, None)

        pbar = tqdm(desc=f"Processing {filename1_no_extension}",
                    unit=" reads",
                    position=position,
                    leave=True,
                    bar_format="{l_bar}{bar} | {n:,} reads [{elapsed}, {rate_fmt}]")

        while r1_next is not None and r2_next is not None:

            id1 = r1_next.id.split()[0].rstrip("/12")
            id2 = r2_next.id.split()[0].rstrip("/12")

            # Skip orphan reads until alignment is restored
            while id1 < id2:
                r1_next = next(r1_iterator, None)
                if r1_next is None:
                    break
                id1 = r1_next.id.split()[0].rstrip("/12")

            while id2 < id1:
                r2_next = next(r2_iterator, None)
                if r2_next is None:
                    break
                id2 = r2_next.id.split()[0].rstrip("/12")

            if r1_next is None or r2_next is None or id1 != id2:
                # Cannot align, skip to next
                continue

            record1 = r1_next
            record2 = r2_next

            total_seqs += 1

            seq1 = analyze_sequence(record1)
            seq2 = analyze_sequence(record2)

            seq_lens1.append(seq1["length"])
            seq_lens2.append(seq2["length"])

            good_reads_r1 += seq1["meets_criteria"]
            good_reads_r2 += seq2["meets_criteria"]

            if seq1["meets_criteria"] and seq2["meets_criteria"]:
                SeqIO.write(record1, filtered1, "fastq")
                SeqIO.write(record2, filtered2, "fastq")
                good_sequences += 1
            else:
                exclusion_reasons_r1["too_short"] += seq1["too_short"]
                exclusion_reasons_r1["found_dot"] += bool(seq1["dot_count"])
                exclusion_reasons_r1["found_n"] += bool(seq1["n_count"])
                exclusion_reasons_r1["found_homopolymer"] += seq1["found_homopolymer"]
                exclusion_reasons_r1["low_score"] += seq1["low_score"]

                exclusion_reasons_r2["too_short"] += seq2["too_short"]
                exclusion_reasons_r2["found_dot"] += bool(seq2["dot_count"])
                exclusion_reasons_r2["found_n"] += bool(seq2["n_count"])
                exclusion_reasons_r2["found_homopolymer"] += seq2["found_homopolymer"]
                exclusion_reasons_r2["low_score"] += seq2["low_score"]

            # Advance both iterators
            r1_next = next(r1_iterator, None)
            r2_next = next(r2_iterator, None)

            pbar.update(1)

        pbar.close()
        # --- END NEW ---

        # Export file with sequence length frequencies
        export_length_frequencies(seq_lens1, filename1_no_extension, output_dir)
        export_length_frequencies(seq_lens2, filename2_no_extension, output_dir)

        # Close opened files
        f_r1.close()
        f_r2.close()
        filtered1.close()
        filtered2.close()

        print(f"{r1_filename.name}: Finished. {good_sequences}/{total_seqs} passed.")

        # Return dict with the paired end's cleaning stats to be compiled in
        # fastfilter_overview.csv file
        return {
            "file": filename1_no_extension.replace("_R1", ""),
            "total_seqs": total_seqs,
            "good_sequences_count": good_sequences,
            "homopolymers": homopolymers,
            "exclusion_reasons_r1": exclusion_reasons_r1,
            "good_reads_r1": good_reads_r1,
            "exclusion_reasons_r2": exclusion_reasons_r2,
            "good_reads_r2": good_reads_r2,
        }


def query_seqDir() -> Path:
    """Interactively query the user for the folder that contains the FASTQ
    files
    """
    try:
        locations = [dir for dir in projects_dir.iterdir() if dir.is_dir()]
    except FileNotFoundError:
        print(
            f"Error! Couldn't find {projects_dir}. Please manually specifiy "
            "the sequences folder using the -i flag. \n\nUse: \n"
            "fastfilter.py -i /path/to/sequences/folder"
        )
        sys.exit()
    print("Which folder contains the input FASTQ files?")
    print(
        "Use your up/down arrows to find the folder you want, and select"
        " it with ENTER."
    )
    questions = [inquirer.List("folder", message="Folder", choices=locations)]
    answer = inquirer.prompt(questions)["folder"]
    return Path(answer)


def generate_reports(result) -> None:
    """Generate 2 report files:
    - fastfilter_overview.csv, stats on filtered sequences per paired-end
    - fastfilter_overview.txt, with processing time and input parameters
    """
    global single_end
    global elapsed_min
    with open(output_dir / "fastfilter_overview.csv", "x", newline="") as f:
        w = csv.writer(f, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        header = [
            "file",
            "total_seqs",
            "good_sequences_count",
            "homopolymers",
            "R1_too_short",
            "R1_found_dot",
            "R1_found_n",
            "R1_found_homopolymer",
            "R1_low_score",
            "R1_good_reads",
        ]
        if not single_end:
            header = header + [
                "R2_too_short",
                "R2_found_dot",
                "R2_found_n",
                "R2_found_homopolymer",
                "R2_low_score",
                "R2_good_reads",
            ]
        w.writerow(header)
        for r in result:
            entry = [
                r["file"],
                r["total_seqs"],
                r["good_sequences_count"],
                r["homopolymers"],
                r["exclusion_reasons_r1"]["too_short"],
                r["exclusion_reasons_r1"]["found_dot"],
                r["exclusion_reasons_r1"]["found_n"],
                r["exclusion_reasons_r1"]["found_homopolymer"],
                r["exclusion_reasons_r1"]["low_score"],
                r["good_reads_r1"],
            ]
            if not single_end:
                entry = entry + [
                    r["exclusion_reasons_r2"]["too_short"],
                    r["exclusion_reasons_r2"]["found_dot"],
                    r["exclusion_reasons_r2"]["found_n"],
                    r["exclusion_reasons_r2"]["found_homopolymer"],
                    r["exclusion_reasons_r2"]["low_score"],
                    r["good_reads_r2"],
                ]
            w.writerow(entry)
    with open(output_dir / "fastfilter_overview.txt", "x", newline="") as f:
        f.write(f"Time elapsed: {elapsed_min}\n")
        f.write(f"Sequence length threshold: {min_seq_len}\n")
        f.write(f"Min. of nucleotides in a homopolymer: {homopolymer_coeff}\n")
        f.write(f"Score threshold: {min_score}\n")


def main():
    multiprocessing.set_start_method("spawn")
    # Save the current time when the script is started
    start = time.time()
    global projects_dir
    # Create a multiple-core processing pool
    projects_dir = PROJECTS_DIR_DEFAULT
    # Parse the arguments from the user
    parse_arguments()

    global output_dir
    global seq_dir
    global single_end
    global elapsed_min

    if seq_dir is None:
        project_dir = query_seqDir()
        seq_dir = project_dir / "cutadapt"
        if seq_dir.is_dir():
            print("Cutadapt folder found! Using the sequences in it")
        else:
            raise Exception(
                "No cutadapt folder found, create one and put the "
                "sequences in it. Or specify a sequences dir "
                "using the -d argument."
            )
    if output_dir is None:
        output_dir = seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Gather all the files ending in '.fastq'
    fastq_r1_files = sorted(seq_dir.glob("*_R1_*.fastq"))
    fastq_r2_files = sorted(seq_dir.glob("*_R2_*.fastq"))

    pool = multiprocessing.Pool(processes=num_cpus)

    print("Running filtering. This may take hours." "Please wait...")

    if len(fastq_r1_files) + len(fastq_r2_files) > 2:
        # List of tuples for matching paired-end reads (R1.fastq, R2.fastq)
        paired_ends = list(zip(fastq_r1_files, fastq_r2_files))
        positions = list(range(len(paired_ends))) # NEW
        # Run the parse_file() function for each of the pairs of FASTQ files
        args_for_pool = [(r1, r2, output_dir, pos) for (r1, r2), pos in zip(paired_ends, positions)] # NEW
        return_values = pool.starmap(parse_file, args_for_pool) # CHANGED
    else:
        single_end = True
        fastq_files = list(seq_dir.glob("*.fastq"))
        return_values = pool.map(parse_file_single_end, fastq_files)
    finish = time.time()
    # Print the time elapsed
    elapsed_min = (finish - start) / 60
    print("Ran successfully in", elapsed_min, "min")
    generate_reports(return_values)


if __name__ == "__main__":
    main()
