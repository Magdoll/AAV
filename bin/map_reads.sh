#!/bin/bash -ex
sample_name=$1
reads=$2
vector_fa=$3
helper_fa=$4
repcap_fa=$5
host_fa=$6

# NOTE: the sequence IDs should be free of blank spaces and symbols. Stick with numbers,
# alphabet letters, and _ and -. If necessary, rename the sequence IDs in the combined
# fasta file.

# STEPS:
# - cat the reference fasta sequences
# - [enh] fix reference fasta sequence names
# - align reads to references OUTPUT
# - BAM, aligned to all references and sorted by read name (?!) SOMEWHERE
# - extract all references' names and lengths (.genome, via faidx)
# DOWNSTREAM:
# - take vector annotations (BED or TXT)
# - append reference names and lengths to vector annotations TXT

ls -Alh

function get_reference_names {
    fasta=$1
    typename=$2
    grep '^>' "$fasta" \
        | cut -f1 -d\  | cut -f2 -d\> \
        | awk  '{ print $1"\t'$typename'" }'
}

get_reference_names "$vector_fa" "vector" > reference_names.tsv
if [ -e "$helper_fa" ]; then
    get_reference_names "$helper_fa" "helper" >> reference_names.tsv
fi
if [ -e "$repcap_fa" ]; then
    get_reference_names "$repcap_fa" "repcap" >> reference_names.tsv
fi
if [ -e "$host_fa" ]; then
    get_reference_names "$host_fa" "host" >> reference_names.tsv
fi
cat reference_names.tsv
echo

cat "$vector_fa" "$helper_fa" "$repcap_fa" "$host_fa" > all_refs.fa
grep '^>' all_refs.fa
echo

samtools fastq -n -0 reads.fq "$reads"
minimap2 --eqx -a --secondary=no -t $(nproc) all_refs.fa reads.fq > mapped.sam
samtools sort -n -O SAM mapped.sam > "$sample_name.sort_by_name.sam"

ls -Alh
