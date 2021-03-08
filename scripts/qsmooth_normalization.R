library(qsmooth)
library(preprocessCore)

args <- commandArgs(trailingOnly=TRUE) 
count <- args[1]
group <- args[2]

prefix <- substring(count, 1, nchar(count) - 9)
out_file_name_qsmooth <- paste0(prefix, 'qsmooth.tsv')
out_file_name_quantile <- paste0(prefix, 'quantile.tsv')

count <- read.table(count, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
group  <- read.table(group, sep='\t', header=FALSE, stringsAsFactors=FALSE)

rownames(group) <- group[, 1]
group <- group[colnames(count), ]

qsmooth_model <- qsmooth(count, group[, 2])
norm <- qsmoothData(qsmooth_model)
rownames(norm) <- rownames(count)
colnames(norm) <- colnames(count)

write.table(
    norm, 
    out_file_name_qsmooth,
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)

norm <- normalize.quantiles(as.matrix(count))
rownames(norm) <- rownames(count)
colnames(norm) <- colnames(count)

write.table(
    norm, 
    out_file_name_quantile,
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)