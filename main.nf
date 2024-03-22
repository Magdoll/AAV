#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process make_all() {

    input:
    path mapped_reads_sam
    path annotation_txt
    val flipflop_name
    val sample_name

    output:
    // summarize alignment
    path("*.per_read.csv"), emit: per_read_csv
    path("*.summary.csv"), emit: summary_csv
    path("*.nonmatch_stat.csv.gz"), emit: nonmatch_stat_csvgz
    path("*.tagged.bam"), emit: tagged_bam
    // flip-flop
    path("*.flipflop_assignments.txt"), emit: flipflop_assignments_txt, optional: true
    // report
    path("*.alignments.tsv"), emit: alignments_tsv
    path("*.readsummary.tsv"), emit: readsummary_tsv
    path("*.sequence-error.tsv"), emit: sequence_error_tsv
    path("*.flipflop.tsv"), emit: flipflop_tsv, optional: true
    path("*.Rdata"), emit: rdata, optional: true
    path("*_AAV_report.html"), emit: aav_report_html
    path("*_AAV_report.pdf"), emit: aav_report_pdf

    script:
    """
    make_all.sh ${mapped_reads_sam} ${annotation_txt} ${sample_name} ${flipflop_name}
    """
}


workflow {
    make_all(
        Channel.fromPath(params.mapped_reads_sam),
        Channel.fromPath(params.annotation_txt),
        params.flipflop_name,
        file(params.mapped_reads_sam).getSimpleName())
}
