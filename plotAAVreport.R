
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
input.prefix <- args[1]
print(input.prefix)
pdf.report.file <- paste0(input.prefix, "_AAV_report.pdf")

x.summary <- read.table(paste0(input.prefix, '.summary.csv'),sep='\t',header=T)
x.err <- read.table(paste0(input.prefix, '.nonmatch_stat.csv'),sep='\t',header=T)
x.read <- read.table(paste0(input.prefix, '.per_read.csv'),sep='\t',header=T)


total_err <- dim(x.err)[1]
x.err$pos0_div <- (x.err$pos0%/%10 * 10)
df.err <- x.err %>% group_by(pos0_div, type) %>% summarise(count=n())
x.err$type_len_cat <- "1-10"
x.err[x.err$type_len>10, "type_len_cat"] <- "11-100"
x.err[x.err$type_len>100, "type_len_cat"] <- "100-500"
x.err[x.err$type_len>500, "type_len_cat"] <- ">500"
x.err$type_len_cat <- ordered(x.err$type_len_cat, levels=c('1-10', '11-100', '100-500', '>500'))
df.err_len_cat <- x.err %>% group_by(type, type_len_cat) %>% summarise(count=n()) %>% mutate(freq=round(100*count/total_err, 2))


min_r_start <- min((x.summary$map_start0%/%100) * 100, na.rm=T)
if (min_r_start == 0) { min_r_start <- -100 }
max_r_end   <- max((x.summary$map_end1%/%100+1) * 100, na.rm=T)

total_read_count <- dim(x.read)[1]
df.read <- x.read %>% group_by(assigned_type, assigned_subtype) %>% summarise(count=n()) %>% mutate(freq=round(count*100/total_read_count,2))


p1.map_iden <- ggplot(x.summary, aes(map_iden*100)) + geom_histogram(fill='black', binwidth=0.01) +
               xlab("Mapping Identity (%)") + ylab("Read Count") +
               labs(title="Distribution of Mapped Identity to Reference")


p1.err_dot <- ggplot(x.err, aes(x=pos0+1, y=type_len)) + geom_point(aes(color=type), alpha=0.5) +
               xlim(c(min_r_start, max_r_end)) +
               xlab("Reference Position") + ylab("Sub/Ins/Del Length") +
               labs(title="Distribution of Non-Matches", subtitle="Each point is a non-match from a read")

p1.err_dot_close <- ggplot(x.err, aes(x=pos0+1, y=type_len)) + geom_point(aes(color=type), alpha=0.5) +
               xlim(c(min_r_start, max_r_end)) +
               ylim(c(0, 100)) +
               xlab("Reference Position") + ylab("Sub/Ins/Del Length") +
               labs(title="Distribution of Non-Matches (of sizes <100 only)", subtitle="Each point is a non-match from a read")

p1.err_sub <- ggplot(filter(df.err,type=='X'), aes(x=pos0_div, y=count)) + geom_bar(fill='darkgreen', stat='identity') +
              xlim(c(min_r_start, max_r_end)) +
              labs(title="Distribution of Non-matches by Reference Position, Substitutions",
                   subtitle="Higher bars indicate hot spots for substitutions w.r.t reference") +
               xlab("Reference Position") + ylab("Frequency")

p1.err_del <- ggplot(filter(df.err,type=='D'), aes(x=pos0_div, y=count)) + geom_bar(fill='darkred', stat='identity') +
              xlim(c(min_r_start, max_r_end)) +
              labs(title="Distribution of Non-matches by Reference Position, Deletions",
                   subtitle="Higher bars indicate hot spots for deletion w.r.t reference") +
               xlab("Reference Position") + ylab("Frequency")

p1.err_ins <- ggplot(filter(df.err,type=='I'), aes(x=pos0_div, y=count)) + geom_bar(fill='darkblue', stat='identity') +
              xlim(c(min_r_start, max_r_end)) +
              labs(title="Distribution of Non-matches by Reference Position, Insertions",
                   subtitle="Higher bars indicate hot spots for insertion w.r.t reference") +
               xlab("Reference Position") + ylab("Frequency")

p1.map_len <- ggplot(x.summary, aes(map_len)) + geom_histogram(aes(y=..count../sum(..count..))) +
               xlab("Mapped Reference Length") + ylab("Fraction of Reads") +
               labs(title="Distribution of Mapped Reference Length")

p1.map_starts <- ggplot(x.summary, aes(map_start0)) + geom_histogram(fill='orange', aes(y=..count../sum(..count..))) +
                xlim(c(min_r_start, max_r_end)) +
                xlab("Mapped Reference Start Position") + ylab("Fraction of Reads") +
                labs(title="Distribution of Mapped Reference Start Position")

p1.map_ends <- ggplot(x.summary, aes(map_end1)) + geom_histogram(fill='magenta', aes(y=..count../sum(..count..))) +
                xlim(c(min_r_start, max_r_end)) +
                xlab("Mapped Reference End Position") + ylab("Fraction of Reads") +
                labs(title="Distribution of Mapped Reference End Position")


p2.atype_violin <- ggplot(x.read, aes(x=assigned_type, y=read_len)) + geom_violin() +
                    xlab("Assigned AAV Type") + ylab("Read Length") +
                    labs(title="Distribution of Read Lengths by Assigned AAV Type")


  pdf(file=pdf.report.file, width = 6.5, height = 6.5)

  #cover
  grid.newpage()
  cover <- textGrob("AAV Report",
                    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
  grid.draw(cover)
  grid.arrange(p1.map_starts, p1.map_ends, p1.map_len, ncol=1)

  grid.arrange(p1.err_sub, p1.err_del, p1.err_ins, ncol=1)
  grid.arrange(p1.map_iden, p1.err_dot, p1.err_dot_close)


    table.err_len_cat <- tableGrob(df.err_len_cat, rows = NULL, cols = c("Err Type", "Err Length", "Count", "Frequency (%)"))
    title.err_len_cat <- textGrob("Length Distribution of Different Non-matches", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.err_len_cat <- gTree(children=gList(title.err_len_cat, table.err_len_cat))
    grid.arrange(gt.err_len_cat)

    table.atype <- tableGrob(df.read, rows = NULL, cols = c("Assigned Type", "Assigned Subtype", "Count", "Frequency (%)"))
    title.atype <- textGrob("Assigned AAV Types By Read Alignment Characteristics", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.atype <- gTree(children=gList(title.atype, table.atype))
  grid.arrange(gt.atype, p2.atype_violin)
  dev.off()
