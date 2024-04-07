#!/usr/bin/env python3
import gzip
import os
import pdb
import random
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from csv import DictReader, DictWriter
from multiprocessing import Process

import pysam

CIGAR_DICT = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X",
    9: "B",
}

SUMMARY_FIELDS = [
    "read_id",
    "read_len",
    "is_mapped",
    "is_supp",
    "map_name",
    "map_start0",
    "map_end1",
    "map_len",
    "map_iden",
    "map_type",
    "map_subtype",
]

PER_READ_FIELDS = [
    "read_id",
    "read_len",
    "has_primary",
    "has_supp",
    "assigned_type",
    "assigned_subtype",
    "effective_count",
]

NONMATCH_FIELDS = ["read_id", "pos0", "type", "type_len"]

name_map_scAAV = {
    "full": "full",
    "left-partial": "left-partial",  #'wtITR-partial',
    "right-partial": "right-partial",  #'mITR-partial',
    "partial": "partial",
    "backbone": "backbone",
    "vector+backbone": "vector+backbone",
}

annot_rex = re.compile("NAME=(\S+);TYPE=([a-zA-Z]+);(REGION=\d+\-\d+){0,1}")
ccs_rex = re.compile("\S+\/\d+\/ccs(\/fwd|\/rev)?")
ANNOT_TYPE_PRIORITIES = {"vector": 1, "repcap": 2, "helper": 3, "lambda": 4, "host": 5}

MAX_DIFF_W_REF = 100
TARGET_GAP_THRESHOLD = 200  # skipping through the on-target region for more than this is considered "full-gap"
DEBUG_GLOBAL_FLAG = False


def subset_sam_by_readname_list(
    in_bam,
    out_bam,
    per_read_csv,
    wanted_types,
    wanted_subtypes,
    max_count=None,
    exclude_subtype=False,
    exclude_type=False,
):
    qname_list = {}  # qname --> (a_type, a_subtype)
    for r in DictReader(open(per_read_csv), delimiter="\t"):
        # pdb.set_trace()
        if (
            wanted_types is None
            or (not exclude_type and r["assigned_type"] in wanted_types)
            or (exclude_type and r["assigned_type"] not in wanted_types)
        ) and (
            wanted_subtypes is None
            or (not exclude_subtype and (r["assigned_subtype"] in wanted_subtypes))
            or (exclude_subtype and (r["assigned_subtype"] not in wanted_subtypes))
        ):
            qname_list[r["read_id"]] = (r["assigned_type"], r["assigned_subtype"])

    cur_count = 0
    reader = pysam.AlignmentFile(in_bam, "rb", check_sq=False)
    writer = pysam.AlignmentFile(out_bam, "wb", header=reader.header)
    for r in reader:
        if r.qname in qname_list:
            d = add_assigned_types_to_record(r, *qname_list[r.qname])
            writer.write(pysam.AlignedSegment.from_dict(d, reader.header))
            cur_count += 1
            if max_count is not None and cur_count >= max_count:
                break
    reader.close()
    writer.close()


def iter_cigar(rec):
    # first we exclude cigar front/end that is hard-clipped (due to supp alignment)
    cigar_list = rec.cigar
    if CIGAR_DICT[cigar_list[0][0]] == "H":
        cigar_list = cigar_list[1:]
    if CIGAR_DICT[cigar_list[-1][0]] == "H":
        cigar_list = cigar_list[:-1]

    # warning_printed = False
    for _type, _count in cigar_list:
        x = CIGAR_DICT[_type]
        if x in ("M", "=", "X", "I", "D", "S", "N"):
            for i in range(_count):
                yield x, _count
        else:
            raise Exception("Unexpected cigar {0}{1} seen! Abort!".format(_count, x))


def iter_cigar_w_aligned_pair(rec, writer):
    ii = iter_cigar(rec)
    prev_cigar_type = None
    prev_r_pos = 0
    total_err = 0
    total_len = 0
    for q_pos, r_pos in rec.get_aligned_pairs():
        cigar_type, cigar_count = next(ii)
        if cigar_type == "S":  # nothing to do if soft-clipped, r_pos must be None
            assert r_pos is None
            continue
        total_len += cigar_count
        if cigar_type != prev_cigar_type:
            if cigar_type in ("I", "D", "X", "N"):
                total_err += cigar_count
                info = {
                    "read_id": rec.qname,
                    "pos0": r_pos if cigar_type != "I" else prev_r_pos,
                    "type": cigar_type,
                    "type_len": cigar_count,
                }
                writer.writerow(info)
            prev_cigar_type = cigar_type
        if r_pos is not None:
            prev_r_pos = r_pos
    return total_err, total_len


def read_annotation_file(annot_filename):
    """
    example
    NAME=chr1;TYPE=host;
    NAME=chr2;TYPE=host;
    NAME=myVector;TYPE=vector;REGION=1795-6553;
    NAME=myCapRep;TYPE=repcap;REGION=1895-5987;
    NAME=myHelper;TYPE=helper;

    :param annot_filename: Annotation file following the format indicated above. Only "vector" is required. Others optional.
    :return:
    """
    d = {}
    for line in open(annot_filename):
        stuff = line.strip()
        m = annot_rex.match(stuff)
        if m is None:
            raise Exception(
                "{0} is not a valid annotation line! Should follow format `NAME=xxxx;TYPE=xxxx;REGION=xxxx;`. Abort!".format(
                    stuff
                )
            )
            sys.exit(-1)

        _name = m.group(1)
        _type = m.group(2)
        _region = (
            None
            if m.group(3) is None
            else tuple(map(int, m.group(3).split("=")[1].split("-")))
        )
        if _type in d:
            raise Exception(
                "Annotation file has multiple {0} types. Abort!".format(_type)
            )
            sys.exit(-1)
        elif _type not in ANNOT_TYPE_PRIORITIES:
            raise Exception(
                "{0} is not a valid type (host, repcap, vector, helper). Abort!".format(
                    _type
                )
            )
            sys.exit(-1)
        else:
            d[_name] = {"type": _type, "region": _region}
    return d


def is_on_target(r, valid_start, valid_end):
    """
    Possible assign types:
     - full (within valid_start, valid_end)
     - backbone (outside valid_start, valid_end)
     - left-partial (incomplete only on the 3'/right end)
     - right-partial (incomplete only on the 5'/left end)
     - partial (incomplete both on the 5' and 3' end)
     - vector+backbone (any amount that crosses over on target and backbone)

    NOTE: at the calling of this method we don't know if it's scAAV/ssAAV yet
    So later we will further split subtype assignment
    ssAAV: full|left-partial|right-partial|partial
    //scAAV: full|wtITR-partial|mITR-partial|partial
    """
    if r.reference_end < valid_start or r.reference_start > valid_end:
        return "backbone"

    diff_start = r.reference_start - valid_start
    diff_end = valid_end - r.reference_end
    # this is what a true "full length" (without large deletions) size would be
    valid_len = valid_end - valid_start

    if abs(diff_start) <= MAX_DIFF_W_REF:  # complete 5' start/left
        if abs(diff_end) <= MAX_DIFF_W_REF:  # complete 3' end/right
            for cigar_type, num in iter_cigar(r):
                if cigar_type == "N" and num >= TARGET_GAP_THRESHOLD:
                    # pdb.set_trace()
                    return "full-gap"
            return "full"
        elif diff_end > MAX_DIFF_W_REF:  # left-partial (incomplete on right/3')
            return "left-partial"
        elif diff_end < -MAX_DIFF_W_REF:  # into backbone
            return "vector+backbone"
    elif diff_start > MAX_DIFF_W_REF:  # right-partial (incomplete on left/5')
        if abs(diff_end) <= MAX_DIFF_W_REF:
            return "right-partial"
        elif diff_end > MAX_DIFF_W_REF:
            return "partial"
        elif diff_end < -MAX_DIFF_W_REF:
            return "vector+backbone"
    else:  # diff < -MAX_DIFF_W_REF, into backbone
        return "vector+backbone"


def assign_read_type(r, annotation):
    """
    :param read_dict: dict of {'supp', 'primary'}
    :return: assigned_type, which could be (scAAV, ssAAV, unclassified) + (super, full, partial, unclassified)

    <assigned_type: ssAAV, scAAV, backbone, helper, repcap, host, can use “+” sign>,
    <assigned_subtype: full, partial, nonAAV>
    <map_stat: unmapped | fully_aligned | partial_aligned | chimeric_aligned>,
    <map to: comma-delimited list of [chr:start-end]>,
    <comma-delimited list of unmapped portion, if any>,
    """
    _type = annotation[r.reference_name]["type"]
    if annotation[r.reference_name]["region"] is not None:
        return _type, is_on_target(r, *annotation[r.reference_name]["region"])
    else:
        return _type, "NA"


def process_alignment_bam(
    sorted_sam_filename,
    annotation,
    output_prefix,
    starting_readname=None,
    ending_readname=None,
):
    """
    :param sorted_sam_filename: Sorted (by read name) SAM filename
    :param annotation:
    :param output_prefix:
    """
    f1 = open(output_prefix + ".summary.csv", "w")
    f2 = open(output_prefix + ".nonmatch_stat.csv", "w")
    f3 = open(output_prefix + ".per_read.csv", "w")

    writer1 = DictWriter(f1, SUMMARY_FIELDS, delimiter="\t")
    writer2 = DictWriter(f2, NONMATCH_FIELDS, delimiter="\t")
    writer3 = DictWriter(f3, PER_READ_FIELDS, delimiter="\t")
    writer1.writeheader()
    writer2.writeheader()
    writer3.writeheader()

    debug_count = 0

    reader = pysam.AlignmentFile(sorted_sam_filename, check_sq=False)
    bam_writer = pysam.AlignmentFile(
        output_prefix + ".tagged.bam", "wb", header=reader.header
    )

    if starting_readname is not None:
        # progress forward until we get to the read
        while True:
            try:
                cur_r = next(reader)
                if cur_r.qname == starting_readname:
                    records = [cur_r]
                    break
            except StopIteration:
                break
    else:
        records = [
            next(reader)
        ]  # records will hold all the multiple alignment records of the same read

    while True:
        try:
            cur_r = next(reader)
            if ending_readname is not None and cur_r.qname == ending_readname:
                break
            if cur_r.qname != records[-1].qname:
                process_alignment_records_for_a_read(
                    records, annotation, writer1, writer2, writer3, bam_writer
                )
                records = [cur_r]
            else:
                records.append(cur_r)
        except StopIteration:  # finished reading the SAM file
            break

    process_alignment_records_for_a_read(
        records, annotation, writer1, writer2, writer3, bam_writer
    )
    bam_writer.close()
    f1.close()
    f2.close()
    f3.close()
    return f3.name, output_prefix + ".tagged.bam"


MIN_PRIM_SUPP_COV = 0.8  # at minimum the total of prim + main supp should cover this much of the original sequence


def find_companion_supp_to_primary(prim, supps):
    """
    Return the most likely companion supp to the primary
    :param prim: the primary info
    :param supps: the list of supp info
    :return: return the most likely companion supp to the primary
    """

    def get_true_start_end(rec, true_qlen):
        # first we need to look at the strand
        # then also look at clipping
        cigartype, cigarlen = rec.cigartuples[0]
        offset = cigarlen if CIGAR_DICT[cigartype] == "H" else 0
        if rec.is_reverse:  # on - strand
            # we need to know the true length
            return true_qlen - (rec.qend + offset), true_qlen - (rec.qstart + offset)
        else:  # on + strand
            # just need to look at clipping
            return rec.qstart + offset, rec.qend + offset

    # if prim['rec'].qname=='m64011_220616_211638/9503552/ccsfwd':
    # pdb.set_trace()
    supp_should_be_rev = not prim["rec"].is_reverse
    # first look for a +/- supp
    for supp in supps:
        if supp["rec"].is_reverse == supp_should_be_rev:
            prim_start, prim_end = get_true_start_end(prim["rec"], prim["read_len"])
            supp_start, supp_end = get_true_start_end(supp["rec"], supp["read_len"])
            min_start = min(prim_start, supp_start)
            max_end = max(prim_end, supp_end)
            if (max_end - min_start) >= prim["read_len"] * MIN_PRIM_SUPP_COV:
                return supp, "+/-"

    # if that didn't work, check if there's a +/+ supp
    for supp in supps:
        if supp["rec"].is_reverse == prim["rec"].is_reverse:
            prim_start, prim_end = get_true_start_end(prim["rec"], prim["read_len"])
            supp_start, supp_end = get_true_start_end(supp["rec"], supp["read_len"])
            min_start = min(prim_start, supp_start)
            max_end = max(prim_end, supp_end)
            if (max_end - min_start) >= prim["read_len"] * MIN_PRIM_SUPP_COV:
                return supp, "+/+"
    return None, None


def add_assigned_types_to_record(r, a_type, a_subtype):
    """
    Add BAM tags
    AT tag <type:scAAV|ssAAV|unclassified>
    AS tag <type:>
    AX tag which is "AT-AX"
    """
    d = r.to_dict()
    d["tags"].append("AT:Z:" + a_type)
    d["tags"].append("AS:Z:" + a_subtype)
    d["tags"].append("AX:Z:" + a_type + "-" + a_subtype)
    return d


def process_alignment_records_for_a_read(
    records, annotation, writer1, writer2, writer3, bam_writer
):
    """
    For each, find the most probable assignment, prioritizing vector > rep/cap > helper > host

    :param records: list of alignment records for the same read
    :return:
    """
    read_tally = {"primary": None, "supp": []}
    for r in records:
        # check ccs id format is <movie>/<zmw>/ccs[/rev or /fwd]
        if ccs_rex.fullmatch(r.qname) is None:
            print(
                "WARNING: sequence ID does not follow format movie/zmw/ccs[/rev or /fwd]. Might undercount ssAAV!"
            )

        info = {
            "read_id": r.qname,
            "read_len": r.query_length,
            "is_mapped": "N" if r.is_unmapped else "Y",
            "is_supp": "NA",
            "rec": r,  # we won't write this out later, it's a holder here for processing prim v supp
            "map_name": "NA",
            "map_start0": "NA",
            "map_end1": "NA",
            "map_len": "NA",
            "map_iden": "NA",
            "map_type": "NA",
            "map_subtype": "NA",
        }
        if r.is_unmapped:
            read_tally["primary"] = info
        else:
            cigar_list = r.cigar
            seq_len = r.query_length
            if CIGAR_DICT[cigar_list[0][0]] == "H":
                seq_len += cigar_list[0][1]
            if CIGAR_DICT[cigar_list[-1][0]] == "H":
                seq_len += cigar_list[-1][1]

            info["map_name"] = r.reference_name
            info["read_len"] = seq_len
            info["is_supp"] = "Y" if r.is_supplementary else "N"
            info["map_start0"] = r.reference_start
            info["map_end1"] = r.reference_end
            info["map_len"] = r.reference_end - r.reference_start
            total_err, total_len = iter_cigar_w_aligned_pair(r, writer2)
            info["map_iden"] = 1 - total_err * 1.0 / total_len

            a_type, a_subtype = assign_read_type(r, annotation)
            info["map_type"] = a_type
            info["map_subtype"] = a_subtype
            if DEBUG_GLOBAL_FLAG:
                print(r.qname, a_type, a_subtype)
            # pdb.set_trace()

            if r.is_supplementary:
                read_tally["supp"].append(info)
            else:
                assert read_tally["primary"] is None
                read_tally["primary"] = info
        # writer1.writerow(info) # not writing here -- writing later when we rule out non-compatible subs

    # summarize it per read, now that all relevant alignments have been processed
    prim = read_tally["primary"]
    supps = read_tally["supp"]

    if len(supps) == 0:
        supp = None
        supp_orientation = None
    elif len(supps) >= 1:  # there's multiple supp, find the companion matching supp
        supp, supp_orientation = find_companion_supp_to_primary(prim, supps)
        # supp could be None, in which case there is best matching supp!
        # in the case supp is None we wanna see if this is a weird read (ex: mapped twice to + strand)

    # write the assigned type / subtype to the new BAM output
    bam_writer.write(
        pysam.AlignedSegment.from_dict(
            add_assigned_types_to_record(
                prim["rec"], prim["map_type"], prim["map_subtype"]
            ),
            prim["rec"].header,
        )
    )
    del prim["rec"]
    writer1.writerow(prim)
    if supp is not None:
        bam_writer.write(
            pysam.AlignedSegment.from_dict(
                add_assigned_types_to_record(
                    supp["rec"], supp["map_type"], supp["map_subtype"]
                ),
                supp["rec"].header,
            )
        )
        del supp["rec"]
        writer1.writerow(supp)

    sum_info = {
        "read_id": prim["read_id"],
        "read_len": prim["read_len"],
        "has_primary": prim["is_mapped"],
        "has_supp": "Y" if supp is not None else "N",
        "assigned_type": "NA",
        "assigned_subtype": "NA",
        "effective_count": 1,
    }
    if sum_info["has_primary"] == "Y":
        if prim["map_type"] == "vector":
            if supp is None:  # special case: primary only, maps to vector --> is ssAAV
                # double check the special case where there was supp candidates but no companion
                if len(supps) > 0:
                    sum_type = "unclassified"  # might be a weird case ex: a read covers the region twice as on + strand
                    sum_subtype = prim["map_subtype"]
                else:  # never had any supp candidates, def ssAAV
                    sum_type = "ssAAV"
                    sum_subtype = prim["map_subtype"]
            else:
                if supp["map_type"] == "vector":
                    if supp_orientation == "+/-":
                        # special case, primary+ supp, maps to vector, --> is scAAV
                        sum_type = "scAAV"
                    else:
                        assert supp_orientation == "+/+"
                        sum_type = "tandem"
                    if supp["map_subtype"] == prim["map_subtype"]:
                        if (
                            sum_type == "scAAV"
                        ):  # special case, rename subtype for scAAV
                            sum_subtype = name_map_scAAV[prim["map_subtype"]]
                        else:
                            sum_subtype = prim["map_subtype"]
                    else:
                        sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
                else:  # primary is in vector, supp not in vector
                    sum_type = prim["map_type"] + "|" + supp["map_type"]
                    sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
        else:  # mapping to non-AAV vector region
            if supp is None:
                sum_type = prim["map_type"]
                sum_subtype = prim["map_subtype"]
            elif supp["map_type"] == prim["map_type"]:
                sum_type = prim["map_type"]
                if supp["map_subtype"] == prim["map_subtype"]:
                    sum_subtype = prim["map_subtype"]
                else:
                    sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
            else:
                sum_type = prim["map_type"] + "|" + supp["map_type"]
                sum_subtype = prim["map_subtype"] + "|" + supp["map_subtype"]
        sum_info["assigned_type"] = sum_type
        sum_info["assigned_subtype"] = sum_subtype

    # ###############################################################
    # 'effective_count' - now look at whether this is an ssAAV type
    # ex: <movie>/<zmw>/ccs   means potential two species (effective count of 2)
    # ex: <movie>/<zmw>/ccs/fwd or rev   is just one
    # NOTE: this can still be undercounting becuz not considering thing that are not ssAAV
    # ###############################################################
    if sum_info["assigned_type"] == "ssAAV":
        if sum_info["read_id"].endswith("/ccs"):
            sum_info["effective_count"] = 2
        # elif sum_info['read_id'].endswith('/ccs/fwd') or sum_info['read_id'].endswith('/ccs/rev'):
        #    sum_info['effective_count'] = 1  # not needed, default is to 1

    writer3.writerow(sum_info)
    if DEBUG_GLOBAL_FLAG:
        print(sum_info)
    # pdb.set_trace()


def run_processing_parallel(sorted_sam_filename, d, output_prefix, num_chunks=1):
    reader = pysam.AlignmentFile(open(sorted_sam_filename), check_sq=False)
    readname_list = [next(reader).qname]
    for r in reader:
        if r.qname != readname_list[-1]:
            readname_list.append(r.qname)

    total_num_reads = len(readname_list)
    chunk_size = (total_num_reads // num_chunks) + 1
    print(
        f"Total {total_num_reads} reads, dividing into {num_chunks} chunks of size {chunk_size}..."
    )

    pool = []
    for i in range(num_chunks):
        p = Process(
            target=process_alignment_bam,
            args=(
                sorted_sam_filename,
                d,
                output_prefix + "." + str(i + 1),
                readname_list[i * chunk_size],
                None
                if (i + 1) * chunk_size > total_num_reads
                else readname_list[(i + 1) * chunk_size],
            ),
        )
        p.start()
        pool.append(p)
        print("Going from {0} to {1}".format(i * chunk_size, (i + 1) * chunk_size))
    for i, p in enumerate(pool):
        if DEBUG_GLOBAL_FLAG:
            print(f"DEBUG: Waiting for {i}th pool to finish.")
        p.join()

    # combine the data together for
    # *.nonmatch_stat.csv, *.per_read.csv, *.summary.csv, *.tagged.bam

    # copy the first chunk over
    o = output_prefix + ".1"

    # shutil.copy(o + '.nonmatch_stat.csv', output_prefix + '.nonmatch_stat.csv')
    # shutil.copy(o + '.per_read.csv', output_prefix + '.per_read.csv')
    # shutil.copy(o + '.summary.csv', output_prefix + '.summary.csv')
    # shutil.copy(o + '.tagged.bam', output_prefix + '.tagged.bam')

    f1 = gzip.open(output_prefix + ".nonmatch_stat.csv.gz", "wb")
    f2 = open(output_prefix + ".per_read.csv", "w")
    f3 = open(output_prefix + ".summary.csv", "w")

    with open(o + ".nonmatch_stat.csv") as h:
        for line in h:
            f1.write(line.encode())
    with open(o + ".per_read.csv") as h:
        for line in h:
            f2.write(line)
    with open(o + ".summary.csv") as h:
        for line in h:
            f3.write(line)
    reader = pysam.AlignmentFile(o + ".tagged.bam", "rb", check_sq=False)
    f4 = pysam.AlignmentFile(output_prefix + ".tagged.bam", "wb", template=reader)
    for r in reader:
        f4.write(r)

    if DEBUG_GLOBAL_FLAG:
        print("Combining chunk data...")
    for i in range(1, num_chunks):
        o = output_prefix + "." + str(i + 1)
        with open(o + ".nonmatch_stat.csv") as h:
            h.readline()  # ignore the header
            f1.write(h.read().encode())
        with open(o + ".per_read.csv") as h:
            h.readline()  # ignore the header
            f2.write(h.read())
        with open(o + ".summary.csv") as h:
            h.readline()  # ignore the header
            f3.write(h.read())
        for r in pysam.AlignmentFile(open(o + ".tagged.bam"), "rb", check_sq=False):
            f4.write(r)

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    reader.close()
    # delete the chunk data
    if DEBUG_GLOBAL_FLAG:
        print("Data combining complete. Deleting chunk data.")
    for i in range(num_chunks):
        o = output_prefix + "." + str(i + 1)
        os.remove(o + ".nonmatch_stat.csv")
        os.remove(o + ".per_read.csv")
        os.remove(o + ".summary.csv")
        os.remove(o + ".tagged.bam")

    return output_prefix + ".per_read.csv", output_prefix + ".tagged.bam"


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("sam_filename", help="Sorted by read name SAM file")
    parser.add_argument("annotation_txt", help="Annotation file")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument(
        "--max_allowed_missing_flanking",
        default=100,
        type=int,
        help="Maximum allowed missing flanking bp to be still considered 'full' (default:100)",
    )
    parser.add_argument(
        "--cpus", type=int, default=1, help="Number of CPUs (default: 1)"
    )
    parser.add_argument("--debug", action="store_true", default=False)
    # parser.add_argument("-f", "--random_frac", default=1., type=float, help="default: off. Random fraction of alignments to subsample.")
    # parser.add_argument("-m", "--max_reads", type=int, default=None, \
    #                    help="default: off. Maximum number of records to process. Can use in conjunction with --random_frac")

    args = parser.parse_args()

    if args.debug:
        DEBUG_GLOBAL_FLAG = True

    MAX_DIFF_W_REF = args.max_allowed_missing_flanking

    d = read_annotation_file(args.annotation_txt)
    if args.cpus == 1:
        per_read_csv, full_out_bam = process_alignment_bam(
            args.sam_filename, d, args.output_prefix
        )
    else:
        per_read_csv, full_out_bam = run_processing_parallel(
            args.sam_filename, d, args.output_prefix, num_chunks=args.cpus
        )

    # subset BAM files into major categories for ease of loading into IGV for viewing
    # subset_sam_by_readname_list(in_bam, out_bam, per_read_csv, wanted_types, wanted_subtypes)
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".scAAV-full.tagged.bam",
        per_read_csv,
        ["scAAV"],
        ["full"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".scAAV-partials.tagged.bam",
        per_read_csv,
        ["scAAV"],
        ["partial", "left-partial", "right-partial"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".scAAV-other.tagged.bam",
        per_read_csv,
        ["scAAV"],
        ["partial", "left-partial", "right-partial", "full"],
        exclude_subtype=True,
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".ssAAV-full.tagged.bam",
        per_read_csv,
        ["ssAAV"],
        ["full"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".ssAAV-partials.tagged.bam",
        per_read_csv,
        ["ssAAV"],
        ["partial", "left-partial", "right-partial"],
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".ssAAV-other.tagged.bam",
        per_read_csv,
        ["ssAAV"],
        ["partial", "left-partial", "right-partial", "full"],
        exclude_subtype=True,
    )
    subset_sam_by_readname_list(
        full_out_bam,
        args.output_prefix + ".others.tagged.bam",
        per_read_csv,
        ["ssAAV", "scAAV"],
        None,
        exclude_type=True,
    )

    # samtools sort/index the above files
    try:
        subprocess.check_call("samtools --help > /dev/null", shell=True)
    except:
        print("WARNING: unable to call samtools to sort the output BAM files. End.")
        sys.exit(-1)

    o = args.output_prefix
    files = [
        o + ".scAAV-full",
        o + ".scAAV-partials",
        o + ".scAAV-other",
        o + ".ssAAV-full",
        o + ".ssAAV-partials",
        o + ".ssAAV-other",
        o + ".others",
    ]
    for p in files:
        subprocess.check_call(
            f"samtools sort {p}.tagged.bam > {p}.tagged.sorted.bam", shell=True
        )
        subprocess.check_call(f"samtools index {p}.tagged.sorted.bam", shell=True)
