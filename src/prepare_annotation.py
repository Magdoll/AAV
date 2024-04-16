#!/usr/bin/env python3
"""Format vector annotations and other references as PacBio annotation.txt.

The output annotation.txt is a bespoke format required by the PacBio scripts. It looks
like::

    NAME=myVector;TYPE=vector;REGION=1795-6553;
    NAME=myCapRep;TYPE=repcap;REGION=1895-5987;
    NAME=myHelper;TYPE=helper;
    NAME=chr1;TYPE=host;
    NAME=chr2;TYPE=host;

Specifically:

  - NAME is a sequence name/ID used in the reference sequence files that reads will be
    mapped to.
  - TYPE is the source type of the reference sequence -- one of 'vector', 'helper',
    'repcap', 'host', or 'lambda' (not used here).
  - There must be exactly one line where TYPE is 'vector'.
  - The 'vector' line also has a REGION field with the 1-based start and end positions.
    In AAV, this is the payload region including the ITRs.
"""

import argparse
import sys
from collections import namedtuple

AnnRow = namedtuple("AnnRow", "seq_name source_type start1 end")

# From summarize_AAV_alignment.py
ANNOT_TYPES = {"vector", "repcap", "helper", "lambda", "host"}


def read_annotation_bed(fname):
    """Read a UCSC BED4 file and extract the coordinates of a row labeled 'vector'.

    Take 'vector' information from a UCSC BED4 file of construct annotations (0-based
    coordinates), ignoring all rows aside from the one labeled 'vector', and format it
    for annotation.txt.

    Check that the input contains exactly one 'vector' row. Ignore all other BED rows.
    """
    out_row = None
    with open(fname) as infile:
        for line in infile:
            # Require BED4 or more
            seq_name, start0, end, label = line.rstrip().split("\t")[:4]
            if label == "vector":
                if out_row is not None:
                    raise RuntimeError(
                        f"Input {args.annotation_bed} contains more than one row "
                        "labeled 'vector'."
                    )
                out_row = AnnRow(seq_name, label, int(start0) + 1, end)
    if out_row is None:
        raise RuntimeError(
            f"Input {args.annotation_bed} must contain a row labeled 'vector'."
        )
    return out_row


def read_reference_names(fname):
    """Read a 2-column TSV of reference sequence names and source types."""
    with open(fname) as infile:
        for line in infile:
            seq_name, source_type = line.split()
            if source_type not in ANNOT_TYPES:
                raise RuntimeError(
                    f"{source_type} must be one of: vector, repcap, helper, host"
                )
            yield AnnRow(seq_name, source_type, None, None)


def write_annotation_txt(out_fname, vector_row, other_rows):
    """Write PacBio-style annotations to `out_fname`.

    Take the vector annotations and non-'vector' sequence names and source types, format
    it for annotation.txt, and append it to the same output file.

    Skip any duplicate 'vector' sequence label appearing in the other references, and
    catch if a sequence name is reused across multiple source types.
    """
    with open(out_fname, "w+") as outf:
        outf.write("NAME={};TYPE={};REGION={}-{};\n".format(*vector_row))
        seen_seq_names_and_sources = {vector_row.seq_name: vector_row.source_type}
        for orow in other_rows:
            if orow.source_type == "vector":
                if orow.seq_name != vector_row.seq_name:
                    raise RuntimeError(
                        "Source type 'vector' listed in additional "
                        f"reference names, but sequence name {orow.seq_name} does not "
                        f"match the previously given {vector_row.seq_name}"
                    )
                continue
            if orow.seq_name in seen_seq_names_and_sources:
                prev_type = seen_seq_names_and_sources[orow.seq_name]
                if orow.source_type != prev_type:
                    raise RuntimeError(
                        f"Sequence name {orow.seq_name} listed with "
                        f"different source types: first {prev_type}, then "
                        f"{orow.source_type}"
                    )
                continue
            outf.write("NAME={seq_name};TYPE={source_type};\n".format(**orow._asdict()))


if __name__ == "__main__":
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument(
        "annotation_bed",
        help="Vector sequence annotations, as UCSC BED4 with 0-based coordinates.",
    )
    AP.add_argument(
        "reference_names",
        help="Reference sequence names and their sources, in 2 columns.",
    )
    AP.add_argument("-o", "--output", help="Output filename (*.txt).")
    args = AP.parse_args()

    try:
        vec_row = read_annotation_bed(args.annotation_bed)
        otr_rows = read_reference_names(args.reference_names)
        write_annotation_txt(args.output, vec_row, otr_rows)
    except RuntimeError as exc:
        sys.exit(str(exc))
