#!/bin/bash -ex
sample_name=$1
mapped_reads_sam=$2
annotation_txt=$3
flipflop_name=$4

ls -Alh

if [[ $mapped_reads_sam == *.bam ]]; then
    echo "Converting $mapped_reads_sam from BAM to SAM"
    sam_fname="${mapped_reads_sam%%.sam}"
    samtools view -b -o "$sam_fname" "$mapped_reads_sam"
    mapped_reads_sam="$sam_fname"
else
    echo "Reads $mapped_reads_sam appear to be in SAM format"
fi

echo
echo "Starting summarize_AAV_alignment"
summarize_AAV_alignment.py "$mapped_reads_sam" "$annotation_txt" "$sample_name" --cpus $(nproc)
echo "Finished summarize_AAV_alignment"
ls -Alh

if [ -n "$flipflop_name" ]; then
    # ENH: take user-provided ITR sequences other than AAV2
    echo "Assuming $flipflop_name == AAV2"

    echo
    echo "Starting get_flipflop_config"
    get_flipflop_config.py "${sample_name}.tagged.bam" "${sample_name}.per_read.csv" -o "$sample_name"
    echo "Finished get_flipflop_config"
    ls -Alh
    flipflop_assignments="${sample_name}.flipflop_assignments.txt"
else
    echo "Skipping flip/flop analysis"
    flipflop_assignments=""
fi

echo
echo "Starting create_report"
create_report.R "./${sample_name}" "$annotation_txt" "$sample_name" "$flipflop_assignments"
echo "Finished create_report"
ls -Alh
