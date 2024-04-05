# Local tests: Run the core 3 analysis steps on a downsampled copy of the public PacBio scAAV data

small_inputs = small
out_small_pfx = build/small


.PSEUDO: fast slow clean
fast: $(out_small_pfx)_AAV_report.html

slow: $(out_small_pfx).tagged.bam $(out_small_pfx).flipflop_assignments.txt $(out_small_pfx)_AAV_report.html


clean:
	rm -vf build/*
	rm -vf $(out_small_pfx)_AAV_report.html $(out_small_pfx)_AAV_report.pdf

$(out_small_pfx)_AAV_report.html: ../src/create_report.R ../src/report.Rmd $(out_small_pfx).annotation.txt $(out_small_pfx).flipflop_assignments.txt
	$< $(out_small_pfx) $(small_inputs)/annotation.txt pAV-CMV-GFP $(out_small_pfx).flipflop_assignments.txt

$(out_small_pfx).annotation.txt: ../src/prepare_annotation.py $(small_inputs)/annotation.bed $(small_inputs)/reference_names.tsv
	$^ -o $@

$(out_small_pfx).flipflop_assignments.txt: ../src/get_flipflop_config.py $(out_small_pfx).tagged.bam $(out_small_pfx).per_read.csv
	$^ -o $(out_small_pfx)

$(out_small_pfx).tagged.bam: ../src/summarize_AAV_alignment.py $(small_inputs)/scaav-cba-egfp.bam $(small_inputs)/annotation.txt
	$^ $(out_small_pfx) --cpus 4
