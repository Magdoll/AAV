#!/usr/bin/env python3
import os, sys
import random
from csv import DictReader, DictWriter
from collections import defaultdict
import pysam

CIGAR_DICT = {0: 'M',
              1: 'I',
              2: 'D',
              3: 'N',
              4: 'S',
              5: 'H',
              6: 'P',
              7: '=',
              8: 'X',
              9: 'B'}

SUMMARY_FIELDS = ['read_id',
                  'read_len',
                  'is_mapped',
                  'is_supp',
                  'map_start0',
                  'map_end1',
                  'map_len',
                  'map_iden']

PER_READ_FIELDS = ['read_id',
                   'read_len',
                   'has_primary',
                   'has_supp',
                   'assigned_type',
                   'assigned_subtype']

NONMATCH_FIELDS = ['read_id',
                   'pos0',
                   'type',
                   'type_len']

def iter_cigar(rec):
    # first we exclude cigar front/end that is hard-clipped (due to supp alignment)
    cigar_list = rec.cigar
    if CIGAR_DICT[cigar_list[0][0]] == 'H': cigar_list = cigar_list[1:]
    if CIGAR_DICT[cigar_list[-1][0]] == 'H': cigar_list = cigar_list[:-1]

    for _type, _count in cigar_list:
        x = CIGAR_DICT[_type]
        if x in ('M', '=', 'X', 'I', 'D', 'S', 'N'):
            for i in range(_count): yield x, _count
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
        if cigar_type == 'S': # nothing to do if soft-clipped, r_pos must be None
            assert r_pos is None
            continue
        total_len += cigar_count
        if cigar_type != prev_cigar_type:
            if cigar_type in ('I', 'D', 'X'):
                total_err += cigar_count
                info = {'read_id': rec.qname,
                        'pos0': r_pos if cigar_type!='I' else prev_r_pos,
                        'type': cigar_type,
                        'type_len': cigar_count}
                writer.writerow(info)
            prev_cigar_type = cigar_type
        if r_pos is not None: prev_r_pos = r_pos
    return total_err, total_len

MAX_DIFF_W_REF = 100
def assign_read_type(read_dict, valid_ref_start, valid_ref_end):
    """
    :param read_dict: dict of {'supp', 'primary'}
    :return: assigned_type, which could be (scAAV, ssAAV, unknown) + (super, full, partial, unknown)
    """
    r1, r2 = read_dict['primary'], read_dict['supp']
    subtype1, subtype2 = None, None
    if r1 is not None:
        diff_start = r1.reference_start - valid_ref_start
        diff_end = r1.reference_end - valid_ref_end
        if r1.reference_end < valid_ref_start or r1.reference_end > valid_ref_end:
            subtype1 = 'unknown'
        elif (diff_start < -MAX_DIFF_W_REF) or \
                (diff_end > MAX_DIFF_W_REF):
            subtype1 = 'super'
        elif (diff_start > MAX_DIFF_W_REF and diff_end <= 0):
            subtype1 = 'partial'
        else:
            subtype1 = 'full'
    if r2 is not None:
        diff_start = r2.reference_start - valid_ref_start
        diff_end = r2.reference_end - valid_ref_end
        if r2.reference_end < valid_ref_start or r2.reference_end > valid_ref_end:
            subtype2 = 'unknown'
        elif (diff_start < -MAX_DIFF_W_REF) or \
                (diff_end > MAX_DIFF_W_REF):
            subtype2 = 'super'
        elif (diff_start > MAX_DIFF_W_REF and diff_end <= 0):
            subtype2 = 'partial'
        else:
            subtype2 = 'full'

    if subtype1 is None:
        if subtype2 is None: return 'unknown', 'unknown'
        else: return 'unknown', subtype2
    else:
        if subtype2 is None: return 'ssAAV', subtype1
        else: return 'scAAV', subtype1+'-'+subtype2


def process_alignment_bam(bam_filename, output_prefix, valid_ref_start, valid_ref_end, random_frac=None, max_reads=None):
    """
    :param bam_filename: Aligned BAM file
    :param seq_len_dict: dict of seq id --> seq length (we need this becuz hard-clipped lengths aren't retained)
    :param random_frac: (optional) if given a fraction of (0, 1], will subsample the reads, useful for not reading the whole file
    :param max_reads: (optional) if given, will stop after max number of reads have been processed; can be used in conjunction w random_frac
    """
    f1 = open(output_prefix+'.summary.csv', 'w')
    f2 = open(output_prefix+'.nonmatch_stat.csv', 'w')
    f3 = open(output_prefix+'.per_read.csv', 'w')

    writer1 = DictWriter(f1, SUMMARY_FIELDS, delimiter='\t')
    writer2 = DictWriter(f2, NONMATCH_FIELDS, delimiter='\t')
    writer3 = DictWriter(f3, PER_READ_FIELDS, delimiter='\t')
    writer1.writeheader()
    writer2.writeheader()
    writer3.writeheader()

    debug_count = 0
    read_tally = defaultdict(lambda: {'primary': None, 'supp': None})

    for r in pysam.AlignmentFile(bam_filename, check_sq=False):
        if max_reads is not None and debug_count > max_reads: break
        if random_frac is not None and random.random() > random_frac: continue
        debug_count += 1
        info = {'read_id': r.qname,\
                'read_len': r.query_length,\
                'is_mapped': 'N' if r.is_unmapped else 'Y',\
                'is_supp': 'NA',
                'map_start0': 'NA',
                'map_end1': 'NA',
                'map_len': 'NA',
                'map_iden': 'NA'
                }
        if not r.is_unmapped:
            cigar_list = r.cigar
            seq_len = r.query_length
            if CIGAR_DICT[cigar_list[0][0]] == 'H': seq_len += cigar_list[0][1]
            if CIGAR_DICT[cigar_list[-1][0]] == 'H': seq_len += cigar_list[-1][1]

            if r.is_supplementary: read_tally[r.qname]['supp'] = r
            else: read_tally[r.qname]['primary'] = r
            read_tally[r.qname]['read_len'] = seq_len

            info['read_len'] = seq_len
            info['is_supp'] = 'Y' if r.is_supplementary else 'N'
            info['map_start0'] = r.reference_start
            info['map_end1'] = r.reference_end
            info['map_len'] = r.reference_end - r.reference_start
            total_err, total_len = iter_cigar_w_aligned_pair(r, writer2)
            info['map_iden'] = 1 - total_err*1./total_len
        writer1.writerow(info)

    f1.close()
    f2.close()
    # now go through read tallies
    for read_id, d in read_tally.items():
        a_type, a_subtype = assign_read_type(d, valid_ref_start, valid_ref_end)
        info = {'read_id': read_id,
                'read_len': d['read_len'],
                'has_primary': 'Y' if d['primary'] is not None else 'N',
                'has_supp': 'Y' if d['supp'] is not None else 'N',
                'assigned_type': a_type,
                'assigned_subtype': a_subtype}
        writer3.writerow(info)
    f3.close()

    print("Processed {0} records from {1}.".format(debug_count, bam_filename))



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("bam_filename", help="Aligned BAM file")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("--valid_ref_range", default=None, help="Expected sequencing range, ex: (60,2500). Default: off")
    parser.add_argument("-f", "--random_frac", default=1., type=float, help="default: off. Random fraction of alignments to subsample.")
    parser.add_argument("-m", "--max_reads", type=int, default=None, \
                        help="default: off. Maximum number of records to process. Can use in conjunction with --random_frac")

    args = parser.parse_args()
    if args.valid_ref_range is None:
        valid_ref_start, valid_ref_end = 0, 99999999
    else:
        valid_ref_start, valid_ref_end = eval(args.valid_ref_range)
    process_alignment_bam(args.bam_filename, args.output_prefix, valid_ref_start, valid_ref_end, args.random_frac, args.max_reads)
