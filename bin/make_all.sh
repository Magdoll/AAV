#!/bin/bash -ex
mapped_reads_sam=$1
annotation_txt=$2
sample_name=$3

ls -Alh /usr/local/bin

echo
echo "Starting summarize_AAV_alignment"
summarize_AAV_alignment.py "${mapped_reads_sam}" "${annotation_txt}" "${sample_name}" --cpus 4
echo "Finished summarize_AAV_alignment"
ls -Alh

echo
echo "Starting get_flipflop_config"
get_flipflop_config.py "${sample_name}.tagged.bam" "${sample_name}.per_read.csv" -o "${sample_name}"
echo "Finished get_flipflop_config"
ls -Alh

echo
echo "Starting create_report"
create_report.R ./${sample_name} ${annotation_txt} ${sample_name} ${sample_name}.flipflop_assignments.txt
echo "Finished create_report"
ls -Alh
