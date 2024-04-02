## Inputs

### Required

- Mapped reads (BAM/SAM)
  - Sequence file from PacBio Machine run in AAV Mode.
- Vector Annotations (BED)
  - Annotated construct region coordinates including ITR, promotor and CDS regions.
  - Must have either:
	- multiple "ITR" regions (not case sensitive), e.g. "ITR-L" and "ITR-R"
	- a region named "vector" that spans ITR to ITR
- Reference sequences (FASTA)
    - Vector plasmid (Required)
    - Helper plasmid, repCap plasmid (optional)
    - Host genome (optional but recommended)

### Options

- AAV Serotype
  - Serotype name for ITR flip/flop analysis, e.g. "AAV2", or leave blank to skip.
