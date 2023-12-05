library(rmarkdown)

working_dir = getwd()
message(paste("Working directory:", working_dir))

args = commandArgs(trailingOnly = TRUE)

input_prefix = args[1]
annot_filename = args[2] # ex: annotation.txt
sampleid = args[3]
flipflop_summary = ''
if (length(args) > 3) {
    flipflop_summary = args[4]
}

input_params = list(
    input_prefix = input_prefix,
    annot_filename = annot_filename,
    sampleid = sampleid,
    flipflop_summary = flipflop_summary)
message("Parameters:")
print(input_params)

rmd_dir = dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[4]))
rmd_path = paste0(rmd_dir, "/report.Rmd")
outfile = basename(paste0(input_prefix, "_AAV_report.html"))
message(paste("Rmd location:", rmd_path, sep=" "))
message(paste("Output location:", outfile, sep=" "))


rmarkdown::render(
    rmd_path,
    output_format="html_document",
    output_file=outfile,
    output_dir=working_dir,
    knit_root_dir=working_dir,
    params=input_params)
