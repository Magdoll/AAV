#!/usr/bin/env Rscript
library(rmarkdown)

working_dir = getwd()
message(paste("Working directory:", working_dir))

args = commandArgs(trailingOnly = TRUE)

input_prefix = args[1]
annot_filename = args[2] # ex: annotation.txt
sample_id = args[3]
flipflop_summary = ''
if (length(args) > 3) {
    flipflop_summary = args[4]
}

rmd_dir = dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[4]))
rmd_path = paste0(rmd_dir, "/report.Rmd")
message(paste("Report template location:", rmd_path, sep=" "))

out_path = paste0(input_prefix, "_AAV_report")  # Adds .html and .pdf automatically
out_dir = dirname(out_path)
out_filename = basename(out_path)
message(paste("Output location:", out_path, sep=" "))

input_params = list(
    input_prefix = input_prefix,
    annot_filename = annot_filename,
    sample_id = sample_id,
    flipflop_summary = flipflop_summary)
message("Parameters:")
print(input_params)

rmarkdown::render(
    rmd_path,
    output_format="all",
    output_file=out_filename,
    output_dir=out_dir,
    knit_root_dir=working_dir,
    params=input_params)
