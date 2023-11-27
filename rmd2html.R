library(rmarkdown)

args <- commandArgs(trailingOnly=TRUE)

infile=args[1]
outfile=args[2]

render(infile,
       output_format="html_document",
       output_file=outfile)
