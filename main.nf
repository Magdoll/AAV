#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { make_all } from './modules/laava'

workflow {
    make_all(
        Channel.fromPath(params.mapped_reads_sam),
        Channel.fromPath(params.annotation_txt),
        params.flipflop_name,
        file(params.mapped_reads_sam).getSimpleName())
}
