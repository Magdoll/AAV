
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
input.prefix <- args[1]
annot.filename <- args[2] # ex: annotation.txt
print(input.prefix)
pdf.report.file <- paste0(input.prefix, "_AAV_report.pdf")

# read the annotation file to find the vector target region
#
#NAME=myHost;TYPE=host;
#NAME=myVector;TYPE=vector;REGION=1795-6553;
#NAME=mRepCap;TYPE=repcap;REGION=1895-5987;


TARGET_REGION_START <- 0
TARGET_REGION_END <- 0
TARGET_REGION_START_REPCAP <- 0
TARGET_REGION_END_REPCAP <- 0

annot <- read.table(annot.filename)
for (i in 1:dim(annot)[1]) {
    if (unlist(strsplit(annot[i,],';'))[2]=='TYPE=vector') {
        p <- unlist(strsplit(annot[i,],';'))[3];
        s_e <- as.integer(unlist(strsplit(unlist(strsplit(p, '='))[2], '-')));
        TARGET_REGION_START <- s_e[1];
        TARGET_REGION_END <- s_e[2];
    }
    else if (unlist(strsplit(annot[i,],';'))[2]=='TYPE=repcap') {
        p <- unlist(strsplit(annot[i,],';'))[3];
        s_e <- as.integer(unlist(strsplit(unlist(strsplit(p, '='))[2], '-')));
        TARGET_REGION_START_REPCAP <- s_e[1];
        TARGET_REGION_END_REPCAP <- s_e[2];
    }
}


x.all.summary <- read.table(paste0(input.prefix, '.summary.csv'),sep='\t',header=T)
x.all.err <- read.table(paste0(input.prefix, '.nonmatch_stat.csv.gz'),sep='\t',header=T)
x.all.read <- read.table(paste0(input.prefix, '.per_read.csv'),sep='\t',header=T)

x.all.err[x.all.err$type=='D',"type"] <- 'deletion'
x.all.err[x.all.err$type=='I',"type"] <- 'insertion'
x.all.err[x.all.err$type=='X',"type"] <- 'mismatch'
x.all.err[x.all.err$type=='N',"type"] <- 'gaps'

# ----------------------------------------------------
# produce stats for vector only (ssAAV or scAAV)
# ----------------------------------------------------
x.read.vector <- filter(x.all.read, assigned_type %in% c('scAAV', 'ssAAV'))
x.err.vector <- filter(x.all.err, read_id %in% x.read.vector$read_id)
x.summary.vector <- filter(x.all.summary, read_id %in% x.read.vector$read_id)

total_num_reads <- dim(x.read.vector)[1]

total_err <- dim(x.err.vector)[1]
x.err.vector$pos0_div <- (x.err.vector$pos0%/%10 * 10)
df.err.vector <- x.err.vector %>% group_by(pos0_div, type) %>% summarise(count=n())
x.err.vector$type_len_cat <- "1-10"
x.err.vector[x.err.vector$type_len>10, "type_len_cat"] <- "11-100"
x.err.vector[x.err.vector$type_len>100, "type_len_cat"] <- "100-500"
x.err.vector[x.err.vector$type_len>500, "type_len_cat"] <- ">500"
x.err.vector$type_len_cat <- ordered(x.err.vector$type_len_cat, levels=c('1-10', '11-100', '100-500', '>500'))
df.err_len_cat.vector <- x.err.vector %>% group_by(type, type_len_cat) %>% summarise(count=n()) %>% mutate(freq=round(100*count/total_err, 2))

df.read_stat_N <- filter(x.err.vector,type=='gaps') %>% group_by(read_id) %>% summarise(max_del_size=max(type_len))
num_reads_large_del <- sum(df.read_stat_N$max_del_size>=200)
freq_reads_large_del <- round(num_reads_large_del*100/total_num_reads, 2)
df.read_stat_N_summary <- data.frame(category=c("Total Reads", "Reads with gaps >200bp"),
                           value=c(total_num_reads, paste0(num_reads_large_del, " (", freq_reads_large_del, "%)")))


ERR_SAMPLE_SIZE <- 50000
x.err2.vector <- x.err.vector[sample(1:dim(x.err.vector)[1], ERR_SAMPLE_SIZE),]
p1.err_dot <- ggplot(x.err2.vector, aes(x=pos0+1, y=type_len)) + geom_point(aes(color=type), alpha=0.5) +
               xlim(c(TARGET_REGION_START, TARGET_REGION_END)) +
               xlab("Reference Position") + ylab("Sub/Ins/Del Length") +
               labs(title="Distribution of Non-Matches", subtitle="Each point is a non-match from a read, only 50k points at most")

p1.err_dot_close <- ggplot(x.err2.vector, aes(x=pos0+1, y=type_len)) + geom_point(aes(color=type), alpha=0.5) +
               xlim(c(TARGET_REGION_START, TARGET_REGION_END)) +
               ylim(c(0, 100)) +
               xlab("Reference Position") + ylab("Sub/Ins/Del Length") +
               labs(title="Distribution of Non-Matches (of sizes <100 only)", subtitle="Each point is a non-match from a read, only 50k points at most")

p1.err_sub <- ggplot(filter(df.err.vector,type=='mismatch'), aes(x=pos0_div, y=count)) + geom_bar(fill='darkgreen', stat='identity') +
              xlim(c(TARGET_REGION_START, TARGET_REGION_END)) +
              labs(title="Distribution of Non-matches by Reference Position, Substitutions",
                   subtitle="Higher bars indicate hot spots for substitutions w.r.t reference") +
               xlab("Reference Position") + ylab("Frequency")

p1.err_del <- ggplot(filter(df.err.vector,type=='deletion'), aes(x=pos0_div, y=count)) + geom_bar(fill='darkred', stat='identity') +
              xlim(c(TARGET_REGION_START, TARGET_REGION_END)) +
              labs(title="Distribution of Non-matches by Reference Position, Deletions",
                   subtitle="Higher bars indicate hot spots for deletion w.r.t reference") +
               xlab("Reference Position") + ylab("Frequency")

p1.err_ins <- ggplot(filter(df.err.vector,type=='insertion'), aes(x=pos0_div, y=count)) + geom_bar(fill='darkblue', stat='identity') +
              xlim(c(TARGET_REGION_START, TARGET_REGION_END)) +
              labs(title="Distribution of Non-matches by Reference Position, Insertions",
                   subtitle="Higher bars indicate hot spots for insertion w.r.t reference") +
               xlab("Reference Position") + ylab("Frequency")

p1.map_iden <- ggplot(x.summary.vector, aes(map_iden*100, fill=map_subtype)) + geom_histogram(binwidth=0.01) +
               xlab("Mapping Identity (%)") + ylab("Read Count") +
               labs(title="Distribution of Mapped Identity to Reference")

p1.map_len <- ggplot(x.summary.vector, aes(map_len, fill=map_subtype)) + geom_histogram(aes(y=..count../sum(..count..))) +
               xlab("Mapped Reference Length") + ylab("Fraction of Reads") +
               labs(title="Distribution of Mapped Reference Spanning Region Size")

p1.map_starts <- ggplot(x.summary.vector, aes(map_start0+1, fill=map_subtype)) +
                geom_histogram(aes(y=..count../sum(..count..))) +
                geom_vline(xintercept=TARGET_REGION_START, color='red', lty=2) +
                geom_vline(xintercept=TARGET_REGION_END, color='red', lty=2) +
                xlab("Mapped Reference Start Position") + ylab("Fraction of Reads") +
                labs(title="Distribution of Mapped Reference Start Position")

p1.map_ends <- ggplot(x.summary.vector, aes(map_end1, fill=map_subtype)) +
                geom_histogram(aes(y=..count../sum(..count..))) +
                geom_vline(xintercept=TARGET_REGION_START, color='red', lty=2) +
                geom_vline(xintercept=TARGET_REGION_END, color='red', lty=2) +
                xlab("Mapped Reference End Position") + ylab("Fraction of Reads") +
                labs(title="Distribution of Mapped Reference End Position")

x.read.vector$subtype <- x.read.vector$assigned_subtype
x.read.vector[!x.read.vector$subtype %in% c("full", "full-gap", "5-partial", "3-partial", "partial", "vector+backbone"), "subtype"] <- 'other'

p1.scAAV_len_hist <- ggplot(filter(x.read.vector, assigned_type=='scAAV'), aes(x=read_len, color=subtype)) +
                       geom_freqpoly() +
                       xlab("Read length (bp)") +
                       ylab("Count") +
                       labs(title="Distribution of read length, scAAV, by subtype")

p1.ssAAV_len_hist <- ggplot(filter(x.read.vector, assigned_type=='ssAAV'), aes(x=read_len, color=subtype)) +
                       geom_freqpoly() +
                       xlab("Read length (bp)") +
                       ylab("Count") +
                       labs(title="Distribution of read length, ssAAV, by subtype")

# ----------------------------------------------------
# produce stats for repcap (if exists)
# ----------------------------------------------------
x.read.repcap <- filter(x.all.read, assigned_type=='repcap')
#x.err.repcap <- filter(x.all.err, read_id %in% x.read.repcap$read_id)
x.summary.repcap <- filter(x.all.summary, read_id %in% x.read.repcap$read_id)

if (dim(x.read.repcap)[1] > 10) { # only plot if at least 10 reads
    p1.map_len.repcap <- ggplot(x.summary.repcap, aes(map_len, fill=map_subtype)) + geom_histogram(aes(y=..count../sum(..count..))) +
                   xlab("Mapped Reference Length") + ylab("Fraction of Reads") +
                   labs(title="Repcap: Distribution of Mapped Reference Spanning Region Size")

    p1.map_starts.repcap <- ggplot(x.summary.repcap, aes(map_start0+1, fill=map_subtype)) +
                    geom_histogram(aes(y=..count../sum(..count..))) +
                    geom_vline(xintercept=TARGET_REGION_START_REPCAP, color='red', lty=2) +
                    xlab("Mapped Reference Start Position") + ylab("Fraction of Reads") +
                    labs(title="Repcap: Distribution of Mapped Reference Start Position")

    p1.map_ends.repcap <- ggplot(x.summary.repcap, aes(map_end1, fill=map_subtype)) +
                    geom_histogram(aes(y=..count../sum(..count..))) +
                    geom_vline(xintercept=TARGET_REGION_END_REPCAP, color='red', lty=2) +
                    xlab("Mapped Reference End Position") + ylab("Fraction of Reads") +
                    labs(title="Repcap: Distribution of Mapped Reference End Position")
}

allowed_subtypes <- c('full', 'full-gap', 'vector+backbone')
p2.atype_violin <-ggplot(filter(x.read.vector, assigned_subtype %in% allowed_subtypes), aes(x=paste(assigned_type, assigned_subtype,sep='-'), y=read_len)) +
                    geom_violin() +
                    xlab("Assigned AAV Type") + ylab("Read Length") +
                    labs(title="Distribution of Read Lengths by Assigned AAV Type") +
                    theme(axis.text.x=element_text(angle = -45, hjust = 0))


p3.err_Ns <- ggplot(filter(df.err.vector,type=='gaps'), aes(x=pos0_div, y=count)) + geom_bar(fill='orange', stat='identity') +
              xlim(c(TARGET_REGION_START, TARGET_REGION_END)) +
              labs(title="Distribution of large deletion events (cigar 'N'), by position",
                   subtitle="Higher bars indicate hot spots for large deletions w.r.t reference") +
               xlab("Reference Position") + ylab("Frequency")

p3.err_size_Ns <- ggplot(df.read_stat_N, aes(max_del_size)) + geom_histogram(binwidth=100) +
                xlab("Maximum large deletion size") + ylab("Number of Reads") +
                labs(title="Distribution of biggest deletion for reads")


  pdf(file=pdf.report.file, width = 6.5, height = 6.5)

  #cover
  grid.newpage()
  cover <- textGrob("AAV Report",
                    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
  grid.draw(cover)

  valid_types <- c('ssAAV', 'scAAV', 'host', 'repcap', 'helper', 'lambda', 'unmapped', 'chimeric')
  x.all.read[is.na(x.all.read$assigned_type), "assigned_type"] <- 'unmapped'
  x.all.read[grep("|", as.character(x.all.read$assigned_type), fixed=T), "assigned_type"] <- 'chimeric'
  x.all.read[!(x.all.read$assigned_type %in% valid_types), "assigned_type"] <- 'other'

  valid_subtypes <- c('full', 'full-gap', '3-partial', '5-partial', 'partial', 'backbone', 'vector+backbone')
  x.all.read[!(x.all.read$assigned_subtype %in% valid_subtypes), "assigned_subtype"] <- 'other'

  min_show_freq <- 0.01
  total_read_count.all <- dim(x.all.read)[1]
  df.read1 <- x.all.read %>% group_by(assigned_type) %>%
           summarise(count=n()) %>% mutate(freq=round(count*100/total_read_count.all,2))
  df.read1 <- df.read1[order(-df.read1$freq),]
  df.read2 <- x.all.read %>% group_by(assigned_type, assigned_subtype) %>%
           summarise(count=n()) %>% mutate(freq=round(count*100/total_read_count.all,2))
  df.read2 <- df.read2[order(-df.read2$freq),]

  table.atype1 <- tableGrob(df.read1, rows = NULL, cols = c("Assigned Type", "Count", "Frequency (%)"))
  title.atype1 <- textGrob("Assigned Types By Read Alignment Characteristics, overview", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
  gt.atype1 <- gTree(children=gList(title.atype1, table.atype1))
  grid.arrange(gt.atype1)
#   table.atype2 <- tableGrob(filter(df.read2,freq>=min_show_freq), rows = NULL,
#                   cols = c("Assigned Type", "Assigned Subtype", "Count", "Frequency (%)"),
#                   theme=ttheme_minimal(base_size = 8))
#   title.atype2 <- textGrob("Assigned Types, detailed (>1% only)", gp=gpar(fontface="italic", fontsize=12), vjust=-18)
#   gt.atype2 <- gTree(children=gList(title.atype2, table.atype2))
#   grid.arrange(gt.atype2)


  total_read_count.vector <- dim(x.read.vector)[1]
  df.read.vector1 <- x.read.vector %>% group_by(assigned_type) %>%
         summarise(count=n()) %>% mutate(freq=round(count*100/total_read_count.vector,2))
  df.read.vector1 <- df.read.vector1[order(-df.read.vector1$freq),]
  df.read.vector2 <- x.read.vector %>% group_by(assigned_type, assigned_subtype) %>%
         summarise(count=n()) %>% mutate(freq=round(count*100/total_read_count.vector,2))
  df.read.vector2 <- df.read.vector2[order(-df.read.vector2$freq),]

  table.atype.vector1 <- tableGrob(df.read.vector1, rows = NULL, cols = c("Assigned Type",  "Count", "Frequency (%)"))
  title.atype.vector1 <- textGrob("Assigned AAV Types By Read Alignment Characteristics, overview", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
  gt.atype.vector1 <- gTree(children=gList(title.atype.vector1, table.atype.vector1))
  grid.arrange(gt.atype.vector1, p2.atype_violin)
  table.atype.vector2 <- tableGrob(df.read.vector2[1:20,], rows = NULL,
                         cols = c("Assigned Type, detailed", "Assigned Subtype", "Count", "Frequency (%)"),
                         theme=ttheme_minimal(base_size = 8))
  title.atype.vector2 <- textGrob("Assigned AAV Types, detailed (top 20 only)", gp=gpar(fontface="italic", fontsize=12), vjust=-24)
  gt.atype.vector2 <- gTree(children=gList(title.atype.vector2, table.atype.vector2))
  grid.arrange(gt.atype.vector2)


  grid.arrange(p1.scAAV_len_hist, p1.ssAAV_len_hist)

  grid.arrange(p1.map_starts, p1.map_ends, p1.map_len, ncol=1)
  if (dim(x.read.repcap)[1] > 10) { # only plot if at least 10 reads
    grid.arrange(p1.map_starts.repcap, p1.map_ends.repcap, p1.map_len.repcap, ncol=1)
  }

  grid.arrange(p1.err_sub, p1.err_del, p1.err_ins, ncol=1)
  grid.arrange(p1.map_iden, p1.err_dot, p1.err_dot_close)


  table.err_len_cat <- tableGrob(df.err_len_cat.vector, rows = NULL, cols = c("Err Type", "Err Length", "Count", "Frequency (%)"))
  title.err_len_cat <- textGrob("Length Distribution of Different Non-matches", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
  gt.err_len_cat <- gTree(children=gList(title.err_len_cat, table.err_len_cat))
  grid.arrange(gt.err_len_cat)

  #only plot the "gap/cigarN" page if it was run with --splice
  if (sum(df.err.vector$type=='gaps')>0) {
    table.err_Ns_summary <- tableGrob(df.read_stat_N_summary, rows=NULL, cols=c("Category", "Count"))
    #title.err_Ns_summary <- textGrob("Reads with large deletions (cigar 'N')", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    #gt.err_Ns_summary <- gTree(children=gList(title.err_Ns_summary, table.err_Ns_summary))
    grid.arrange(p3.err_Ns, p3.err_size_Ns, table.err_Ns_summary, ncol=1)
  }

  dev.off()
