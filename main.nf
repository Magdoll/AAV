#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process make_aav_report() {

    input:
    path annotation_txt
    path flipflop_assignments_txt
    path per_read_csv
    path summary_csv
    path nonmatch_stat_csvgz
    val sample_name

    output:
    path("*.alignments.tsv"), emit: alignments_tsv
    path("*.readsummary.tsv"), emit: readsummary_tsv
    path("*.sequence-error.tsv"), emit: sequence_error_tsv
    path("*.flipflop.tsv"), emit: flipflop_tsv
    //path("*.Rdata"), emit: rdata
    path("*_AAV_report.html"), emit: aav_report_html
    path("*_AAV_report.pdf"), emit: aav_report_pdf
    // e.g.
    // ssAAV_AAV_report.html
    // ssAAV_AAV_report.pdf
    // *these get (over)written twice...
    // ssAAV.alignments.tsv
    // ssAAV.sequence-error.tsv
    // ssAAV.readsummary.tsv
    // ssAAV.flipflop.tsv

    script:
    """
    ls -Alh
    which create_report.R
    create_report.R ${sample_name} ${annotation_txt} ${sample_name} ${flipflop_assignments_txt}
    """
}


process make_read_summary() {

    input:
    path mapped_reads_sam
    path annotation_txt
    val sample_name

    output:
    path("*.per_read.csv"), emit: per_read_csv
    path("*.summary.csv"), emit: summary_csv
    path("*.nonmatch_stat.csv.gz"), emit: nonmatch_stat_csvgz
    path("*.tagged.bam"), emit: tagged_bam
    // e.g.
    // ssAAV.per_read.csv
    // ssAAV.summary.csv
    // ssAAV.nonmatch_stat.csv.gz
    // ssAAV.tagged.bam

    //path("*.flipflop_assignments.txt"), emit: flipflop_assignments_txt, optional: true
    // e.g.  scAAV.flipflop_assignments.txt

    script:
    """
	summarize_AAV_alignment.py \\
        ${mapped_reads_sam} \\
        ${annotation_txt} \\
        ${sample_name}.per_read.csv \\
        --cpus 4
    """
    //get_flipflop_config.py ${tagged_bam} ${per_read_csv} -o ${sample_name}
}


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

    //make_read_summary(
    //    Channel.fromPath(params.mapped_reads_sam),
    //    Channel.fromPath(params.annotation_txt),
    //    file(params.mapped_reads_sam).getSimpleName())

    //make_aav_report(
    //       Channel.fromPath(params.annotation_txt),
    //       Channel.fromPath(params.flipflop_assignments_txt),
    //       Channel.fromPath(params.per_read_csv),
    //       Channel.fromPath(params.summary_csv),
    //       Channel.fromPath(params.nonmatch_stat_csvgz),
    //       file(params.per_read_csv).getSimpleName())

}
