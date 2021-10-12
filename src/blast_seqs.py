# 1. BLAST IF against reference genome
# 2. BLAST IF against unicycler assembly

import csv
import logging
import os
import numpy as np
import pandas as pd
import re

from argparse import ArgumentParser
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

BLAST_DB_FILES = [".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"]


def create_blast_db(fasta):
    for ext in BLAST_DB_FILES:
        if not os.path.exists(fasta + ext):
            cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=fasta)
            logging.info(f"creating BLAST database for {fasta}")
            cline()
            break

 # class Bio.Blast.Record.HSP
 #    Bases: object
 #    Stores information about one hsp in an alignment hit.
 #    Members:
 #            score BLAST score of hit. (float)
 #            bits Number of bits for that score. (float)
 #            expect Expect value. (float)
 #            num_alignments Number of alignments for same subject. (int)
 #            identities Number of identities (int) if using the XML parser. Tuple of number of identities/total aligned (int, int) if using the (obsolete) plain text parser.
 #            positives Number of positives (int) if using the XML parser. Tuple of number of positives/total aligned (int, int) if using the (obsolete) plain text parser.
 #            gaps Number of gaps (int) if using the XML parser. Tuple of number of gaps/total aligned (int, int) if using the (obsolete) plain text parser.
 #            align_length Length of the alignment. (int)
 #            strand Tuple of (query, target) strand.
 #            frame Tuple of 1 or 2 frame shifts, depending on the flavor.
 #            query The query sequence.
 #            query_start The start residue for the query sequence. (1-based)
 #            query_end The end residue for the query sequence. (1-based)
 #            match The match sequence.
 #            sbjct The sbjct sequence.
 #            sbjct_start The start residue for the sbjct sequence. (1-based)
 #            sbjct_end The end residue for the sbjct sequence. (1-based)


def get_blast_output(sample_id, genome, genome_fasta, query_fasta):
    print(str(query_fasta) + " vs. " + genome_fasta)
    seq_id = os.path.basename(os.path.splitext(query_fasta)[0])
    cline = NcbiblastnCommandline(
        query=query_fasta,
        db=genome_fasta,
        evalue=0.001,
        outfmt=5,
        out=f"build/{seq_id}_{genome}.xml",
    )
    logging.info(f"BLASTing {seq_id} against {genome}")
    cline()

    handle = open(f"build/{seq_id}_{genome}.xml")
    record = NCBIXML.read(handle)
    # single return in nested for loops... just returns first best (as requested)?
    # may want to return more
    #   multiples HSPs for same query/database pair
    #   qc
    for algn in record.alignments:
        for hsp in algn.hsps:
            hsp_dict = {
                "sample": sample_id,
                "genome": genome,
                "seq_id": seq_id,
                "s_start": hsp.sbjct_start,
                "s_end": hsp.sbjct_end,
                "expect": hsp.expect,
                "alignments": len(record.alignments),
                "a1_hsps": len((algn.hsps))
            }
            return hsp_dict
    return {"sample": sample_id, "genome": genome, "seq_id": seq_id, "alignments": 0, "a1_hsps": 0}


def remove_float(i):
    if i is None or pd.isnull(i) or pd.isna(i):
        return None
    return str(i).replace(".0", "")


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "reference_genomes", help="CSV containing master sample ID -> reference genome"
    )
    parser.add_argument("output", help="Path to CSV output")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    else:
        logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    # Read in sample->genome map
    with open(args.reference_genomes, "r") as f:
        reader = csv.reader(f)
        next(reader)
        sample_genome = {row[0]: row[1] for row in reader}

    reference_df = pd.DataFrame()
    assembly_df = pd.DataFrame()
    for sample_id, genome in sample_genome.items():
        sample_dir = f"data/iarpa/TE/{sample_id}"
        seq_dir = os.path.join(sample_dir, "target_sequence")
        if not os.path.exists(seq_dir):
            # logging.info(f"sample {sample_id} contains no target sequences")
            continue
        seq_files = [os.path.join(seq_dir, x) for x in os.listdir(seq_dir)]

        # Find the FASTA for the reference genome
        genome_fasta = None
        for f in os.listdir("data/reference_genomes"):
            if re.match(rf"^{genome}.(fa|fasta|fna)$", f):
                genome_fasta = os.path.join("data/reference_genomes", f)
                break
        if not genome_fasta:
            logging.error(f"unable to find FASTA for {genome} in reference genomes")
            continue

        # Check if we need to create the BLAST database for the ref genome
        create_blast_db(genome_fasta)

        assembly_fasta = os.path.join(sample_dir, "assembly/unicycler_assembly.fasta")
        if not os.path.exists(assembly_fasta) and not os.path.islink(assembly_fasta):
            logging.info(f"no unicycler assembly for sample {sample_id} at " + assembly_fasta)
            continue

        # Check if we need to create the BLAST database for the assembly
        create_blast_db(assembly_fasta)

        for sf in seq_files:
            # BLAST against reference genome
            reference_df = reference_df.append(
                get_blast_output(sample_id, genome, genome_fasta, sf), ignore_index=True
            )
            # print(reference_df)
            assembly_df = assembly_df.append(
                get_blast_output(sample_id, "unicycler", assembly_fasta, sf), ignore_index=True
            )
            # print(assembly_df)
            # TODO: if not aligned to assembly, blast against nt database
            # but we will get lots of 100% matches, how do we filter them?

    output_df = pd.merge(reference_df, assembly_df, on=["sample", "seq_id"], how="outer")
    raw_list = output_df.columns
    working_list = raw_list
    working_list = [re.sub(r'(.*)_x$', r'reference_\1', entry) for entry in raw_list]
    working_list = [re.sub(r'(.*)_y$', r'assembly_\1', entry) for entry in working_list]
    rename_dict = {raw_list[i]: working_list[i] for i in range(len(raw_list))}

    output_df = output_df.rename(
        columns=rename_dict
    )

    # output_df = output_df.rename(
    #     columns={
    #         "genome_x": "ref_genome",
    #         "start_x": "ref_start",
    #         "end_x": "ref_end",
    #         "start_y": "assembly_start",
    #         "end_y": "assembly_end",
    #     }
    # )

    # output_df = output_df[
    #     [
    #         "sample",
    #         "ref_genome",
    #         "seq_id",
    #         "ref_start",
    #         "ref_end",
    #         "assembly_start",
    #         "assembly_end",
    #     ]
    # ]
    # output_df["ref_start"] = output_df["ref_start"].apply(remove_float)
    # output_df["ref_end"] = output_df["ref_end"].apply(remove_float)
    # output_df["assembly_start"] = output_df["assembly_start"].apply(remove_float)
    # output_df["assembly_end"] = output_df["assembly_end"].apply(remove_float)
    output_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
