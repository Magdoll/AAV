#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { map_reads; make_report } from './modules/local/laava'

workflow {
    result = map_reads(
        file(params.seq_reads).getSimpleName(),
        Channel.fromPath(params.seq_reads),
        Channel.fromPath(params.vector_fa),
        params.helper_fa ? Channel.fromPath(params.helper_fa) : null,
        params.repcap_fa ? Channel.fromPath(params.repcap_fa) : null,
        params.host_fa ? Channel.fromPath(params.host_fa): null,
    )
    make_report(
        result.sample_name,
        params.flipflop_name,
        result.mapped_reads,
        Channel.fromPath(params.vector_bed),
        result.reference_names,
    )
}
