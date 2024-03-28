## Summary

This workflow can be used with PacBio long-read sequencing to
characterize adeno-associated virus (AAV) products by examining
alignment orientation and coverage over annotated sequence regions. It
produces a report of the alignment statistics in HTML and PDF formats.

## Workflow Diagram

``` mermaid
graph TD

bam(Aligned BAM) --> counts(Read counts)
anno(Annotation) --> counts
counts --> sumtables(Summary tables)
counts --> sumreport(PDF/HTML Report)
```

## Methods

This workflow performs an automated analysis of long-read sequencing
data from adeno-associated virus (AAV) products. The sequencing data
should be from the PacBio sequencer run in AAV mode, or equivalent
circular consensus sequencing (CCS) reads (Travers et al., 2010), and
then aligned to the construct, packaging, and masked host reference
sequences using Minimap2 (Li, 2018).

In this analysis, aligned sequencing reads are filtered for quality to
include primary alignments and reads with mapping quality scores greater
than 10.  The alignment coordinates and orientation of reads passing
these filters are then compared to the annotated vector region in the
reference sequence, which comprises the left and right ITRs and the
genomic region between them, to assign each read to a type and (for AAV
reads) subtype classification according to the definitions above.

Finally a report is generated with relevant quality metrics.

