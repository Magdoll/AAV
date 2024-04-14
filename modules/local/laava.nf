process map_reads() {
    publishDir "$params.output", mode: "copy"

    input:
    val sample_name
    path reads
    path vector_fa
    path helper_fa
    path repcap_fa
    path host_fa

    output:
    path("${sample_name}.sort_by_name.sam"), emit: mapped_reads
    path("reference_names.tsv"), emit: reference_names
    val "${sample_name}", emit: sample_name

    script:
    """
    map_reads.sh ${sample_name} ${reads} ${vector_fa} ${helper_fa} ${repcap_fa} ${host_fa}
    """
}


process make_report() {
    publishDir "$params.output", mode: "copy"

    input:
    val sample_name
    val flipflop_name
    path mapped_reads
    path vector_annotation
    path reference_names

    output:
    // summarize alignment
    path("${sample_name}.per_read.csv"), emit: per_read_csv
    path("${sample_name}.summary.csv"), emit: summary_csv
    path("${sample_name}.nonmatch_stat.csv.gz"), emit: nonmatch_stat_csvgz
    path("${sample_name}.tagged.bam"), emit: tagged_bam
    // flip-flop
    path("${sample_name}.flipflop_assignments.txt"), emit: flipflop_assignments_txt, optional: true
    // report
    path("${sample_name}.alignments.tsv"), emit: alignments_tsv
    path("${sample_name}.readsummary.tsv"), emit: readsummary_tsv
    path("${sample_name}.sequence-error.tsv"), emit: sequence_error_tsv
    path("${sample_name}.flipflop.tsv"), emit: flipflop_tsv, optional: true
    path("${sample_name}.Rdata"), emit: rdata, optional: true
    path("${sample_name}_AAV_report.html"), emit: aav_report_html
    path("${sample_name}_AAV_report.pdf"), emit: aav_report_pdf

    script:
    """
    prepare_annotation.py "${vector_annotation}" "${reference_names}" -o annotation.txt
    make_report.sh "${sample_name}" "${mapped_reads}" annotation.txt "${flipflop_name}"
    """
}
