#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process make_all() {
    input:
    //$(ss_inputs)/mapped.sort_by_read.sam
    path mapped_reads_sam
    //$(ss_inputs)/annotation.txt
    path annotation_txt
    val sample_name

    output:
    path("*.per_read.csv"), emit: per_read_csv
    path("*.summary.csv"), emit: summary_csv
    path("*.nonmatch_stat.csv.gz"), emit: nonmatch_stat_csvgz
    path("*.tagged.bam"), emit: tagged_bam
    path("*.flipflop_assignments.txt"), emit: flip_flop_assignments
    path("*.Rdata"), emit: rdata
    path("*_AAV_report.html"), emit: aav_report_html
    path("*_AAV_report.pdf"), emit: aav_report_pdf
    path("*.alignments.tsv"), emit: alignments_tsv
    path("*.sequence-error.tsv"), emit: sequence_error_tsv
    path("*.readsummary.tsv"), emit: readsummary_tsv
    path("*.flipflop.tsv"), emit: flipflop_tsv

    script:
    """
	summarize_AAV_alignment.py \\
        ${mapped_reads_sam} \\
        ${annotation_txt} \\
        ${sample_name} \\
        --cpus 4

    samtools index ${sample_name}.tagged.bam

	get_flipflop_config.py \\
        ${sample_name}.tagged.bam \\
        ${sample_name}.per_read.csv \\
        -o ${sample_name}

    create_report.R \\
        ./${sample_name}
        ${annotation_txt} \\
        ${sample_name}
        ${sample_name}.flipflop_assignments.txt
    """
}


workflow {
    make_all(
        Channel.fromPath(params.mapped_reads_sam),
        Channel.fromPath(params.annotation_txt),
        file(params.mapped_reads_sam).getSimpleName())
}
