import parasail
import sys
import pysam
from csv import DictReader
from Bio import SeqIO

"""
Get ITR flip flop configurations
Must have already run `summarize_AAV_alignment.py` to get a .tagged.BAM file!

Required:
"""

SW_SCORE_MATRIX = parasail.matrix_create("ACGT", 2, -5)

SEQ_LEFT_FLIP='ttggccactccctctctgcgcgctcgctcgctcactgaggccgggcgaccaaaggtcgcccgacgcccgggctttgcccgggcggcctcagtgagcgagcgagcgcgcagagagggagtggccaactccatcactaggggttcct'.upper()
SEQ_LEFT_FLOP='TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT'
SEQ_RIGHT_FLIP='AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA'
SEQ_RIGHT_FLOP='AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA'

#POS_LEFT_FLIP=0 # start position
#POS_RIGHT_FLIP=4603 # start position

def read_flip_flop_fasta(fasta_filename):
    global SEQ_LEFT_FLIP
    global SEQ_LEFT_FLOP
    global SEQ_RIGHT_FLIP
    global SEQ_RIGHT_FLOP

    flag = 0

    print(f"Reading {fasta_filename}....")

    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        if r.id == 'SEQ_LEFT_FLIP':
            SEQ_LEFT_FLIP = str(r.seq)
            flag += 0b1000
        elif r.id == 'SEQ_LEFT_FLOP':
            SEQ_LEFT_FLOP = str(r.seq)
            flag += 0b0100
        elif r.id == 'SEQ_RIGHT_FLIP':
            SEQ_RIGHT_FLIP = str(r.seq)
            flag += 0b0010
        elif r.id == 'SEQ_RIGHT_FLOP':
            SEQ_RIGHT_FLOP = str(r.seq)
            flag += 0b0001
        else:
            print("WARNING: Sequence IDs must be SEQ_LEFT_FLIP|SEQ_LEFT_FLOP|SEQ_RIGHT_FLIP|SEQ_RIGHT_FLOP. Is {0} instead. Ignoring!".format(r.id))

    # check that all 4 needed sequence IDs are seen
    if (flag >> 3) == 0:
        print("ERROR: SEQ_LEFT_FLIP is not given. Abort!")
        sys.exit(-1)
    flag &= 0b0111
    if (flag >> 2) == 0:
        print("ERROR: SEQ_LEFT_FLOP is not given. Abort!")
        sys.exit(-1)
    flag &= 0b0011
    if (flag >> 1) == 0:
        print("ERROR: SEQ_RIGHT_FLIP is not given. Abort!")
        sys.exit(-1)
    flag &= 0b0001
    if flag == 0:
        print("ERROR: SEQ_RIGHT_FLOP is not given. Abort!")
        sys.exit(-1)


def identify_flip_flop(r):
    """
    Assume record tag:AT is vector, tag:AX can be full|left-partial|right-partial|partial
    Add back a tag 'AF' that is [flip/flop]-[flip/flop]
    """
    t = dict(r.tags)
    try:
        assert t['AX'] in ('vector-full', 'vector-left-partial', 'vector-right-partial', 'vector-partial')
    except AssertionError:
        print("Input BAM records must have a `AX` tag assigned by first running summarize_AAV_alignment.py. Abort!")
        sys.exit(-1)
    config_left, config_right = 'unclassified', 'unclassified'
    if t['AX']=='vector-partial': # ignore, since both sides are missing chunks of ITR
        return 'unclassified', 'unclassified'

    if t['AX'] in ('vector-full', 'vector-left-partial'):
        o1 = parasail.sw_trace(r.query, SEQ_LEFT_FLIP, 3, 1, SW_SCORE_MATRIX)
        o2 = parasail.sw_trace(r.query, SEQ_LEFT_FLOP, 3, 1, SW_SCORE_MATRIX)
        if o1.score > o2.score and o1.score > 250: 
            config_left = 'flip'
        elif o2.score > o1.score and o2.score > 250:
            config_left = 'flop'
        else:
            config_left = 'unclassified'
    
    if t['AX'] in ('vector-full', 'vector-right-partial'):
        o1 = parasail.sw_trace(r.query[-len(SEQ_RIGHT_FLIP)-10:], SEQ_RIGHT_FLIP, 3, 1, SW_SCORE_MATRIX)
        o2 = parasail.sw_trace(r.query[-len(SEQ_RIGHT_FLOP)-10:], SEQ_RIGHT_FLOP, 3, 1, SW_SCORE_MATRIX)
        if o1.score > o2.score and o1.score > 250:
            config_right = 'flip'
        elif o2.score > o1.score and o2.score > 250:
            config_right = 'flop'
        else:
            config_right = 'unclassified'
    return config_left, config_right

def main(per_read_csv, tagged_bam, output_prefix):

    read_info = {}
    for r in DictReader(open(per_read_csv), delimiter='\t'):
        read_info[r['read_id']] = r

    fout = open(output_prefix + '.flipflop_assignments.txt', 'w')
    fout.write("name\ttype\tsubtype\tstart\tend\tleftITR\trightITR\n")
    reader = pysam.AlignmentFile(open(tagged_bam), 'rb', check_sq=False)
    writer1 = pysam.AlignmentFile(open(output_prefix+'.vector-full-flipflop.bam', 'w'), 'wb',
                                  header=reader.header)
    writer2 = pysam.AlignmentFile(open(output_prefix+'.vector-leftpartial-flipflop.bam', 'w'), 'wb',
                                  header=reader.header)
    writer3 = pysam.AlignmentFile(open(output_prefix+'.vector-rightpartial-flipflop.bam', 'w'), 'wb',
                                  header=reader.header)
    for r in reader:
        t = dict(r.tags)
        # if r.qname=='m54278_220522_043945/5964021/ccs/rev': break
        if t['AT'] == 'vector' and t['AX'] in ('vector-full', 'vector-left-partial', 'vector-right-partial'):
            c_l, c_r = identify_flip_flop(r)
            d = r.to_dict()
            a_type = read_info[r.qname]['assigned_type']
            if a_type not in ('scAAV', 'ssAAV'): continue
            d['tags'].append('AF:Z:' + c_l + '-' + c_r)
            d['tags'].append('AG:Z:' + a_type)
            if t['AX'] == 'vector-full':
                writer = writer1
            elif t['AX'] == 'vector-right-partial':
                writer = writer2
            elif t['AX'] == 'vector-left-partial':
                writer = writer3
            writer.write(pysam.AlignedSegment.from_dict(d, r.header))
            fout.write(r.qname + '\t' + a_type + '\t' + t['AX'] + '\t' +
                       str(r.reference_start) + '\t' +
                       str(r.reference_end) + '\t' +
                       c_l + '\t' + c_r + '\n')

    writer1.close()
    writer2.close()
    writer3.close()
    fout.close()

    print("Output summmary: {0}".format(fout.name))
    print(f"Indidual BAM files written: {output_prefix}.vector- full,leftpartial,rightpartial -flipflop.bam")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("sorted_tagged_bam", help="Sorted tagged BAM file")
    parser.add_argument("per_read_csv", help="Per read CSV file")
    parser.add_argument("-o", "--output_prefix", help="Output prefix", required=True)
    parser.add_argument("--flipflop_fasta", default=None, help="(optional) flip flop fasta file (if not given, uses AAV2 default)")

    args = parser.parse_args()

    if args.flipflop_fasta is not None:
        read_flip_flop_fasta(args.flipflop_fasta)

    main(args.per_read_csv, args.sorted_tagged_bam, args.output_prefix)

