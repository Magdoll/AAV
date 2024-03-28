process make_all() {
    publishDir "$params.outdir", mode: "copy"
    container = "ghcr.io/formbio/laava:latest"

    input:
    path mapped_reads_sam
    path annotation_txt
    val flipflop_name
    val sample_name

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
    make_all.sh ${mapped_reads_sam} ${annotation_txt} ${sample_name} ${flipflop_name}
    """
}
