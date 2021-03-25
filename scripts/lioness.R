library(lionessR)
library(stats)
library(SummarizedExperiment)
library(limma)
library(igraph)
library(plotrix)
library(scales)
library(biomaRt)
library(fields)

spearman <- function(x){
    return(cor(t(x), method="spearman"))
}

radian.rescale <- function(x, start=0, direction=1) {
   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

args <- commandArgs(trailingOnly=TRUE) 
ref_file <- args[1]
before_file <- args[2]
after_file <- args[3]
out_file <- args[4]

xprs <- list()
cormat <- list()
for(index in seq(3)){
    xprs[[index]] <- read.table(args[index], sep='\t', header=TRUE, row.names=1)
    
    if(index == 1){
        ensembl <- useMart(
            "ensembl", 
            dataset = "mmusculus_gene_ensembl", 
            host = "uswest.ensembl.org"
        )
        annotation <- getBM(
            attributes = c('ensembl_gene_id', 'uniprot_gn_symbol'),
            filters = 'ensembl_gene_id', 
            values = rownames(xprs[[index]]), 
            mart = ensembl
        )
        annotation <- annotation[annotation$uniprot_gn_symbol != '' & !duplicated(annotation$uniprot_gn_symbol), ]
    }

    num_all_genes <- nrow(xprs[[index]])
    xprs[[index]] <- xprs[[index]][annotation$ensembl_gene_id, ]

    message(num_all_genes - nrow(xprs[[index]]), ' genes without annotation are removed')

    rownames(xprs[[index]]) = annotation$uniprot_gn_symbol
    tmp_cormat <- lioness(as.matrix(xprs[[index]]), spearman)
    cormat[[index]] <- assays(tmp_cormat)[['lioness']]
}

message(nrow(xprs[[index]]), ' genes. ', nrow(xprs[[index]]) * (nrow(xprs[[index]]) - 1) / 2 - nrow(xprs[[index]]), ' total edges.')

edges <- t(matrix(unlist(c(strsplit(row.names(cormat[[1]]), "_"))), 2))

group <- factor(rep(c("Reference", "Target"), each=dim(cormat[[2]])[2]))
design <- model.matrix(~0+group)
cont.matrix <- makeContrasts(groupTarget-groupReference, levels = design)  

fit <- lmFit(cbind(cormat[[1]], cormat[[2]]), design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
toptable_before <- topTable(fit2e, number=nrow(cormat[[2]]), adjust="fdr")

fit <- lmFit(cbind(cormat[[1]], cormat[[3]]), design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
toptable_after <- topTable(fit2e, number=nrow(cormat[[3]]), adjust="fdr")

write.table(
    toptable_before, 
    paste0(out_file, "_limma.tsv"),
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)

write.table(
    toptable_after, 
    paste0(out_file, "_caiman_limma.tsv"),
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)

c_scale <- colorRampPalette(c('dodgerblue', 'snow2', 'violetred2'))

table <- toptable_before[toptable_before$adj.P.Val < 0.05, ]

message(nrow(table) / 2, ' false discovery edges (before correction).')
message(nrow(toptable_after[toptable_after$adj.P.Val < 0.01, ]) / 2, ' false discovery edges (after correction).')

edges <- t(matrix(unlist(c(strsplit(row.names(table), "_"))),2))
z <- cbind(edges, table$logFC)
g <- graph.data.frame(z, directed=F)
max_range <- max(as.numeric(z[,3]))
E(g)$weight <- as.numeric(z[,3]) 
E(g)$color <- head(c_scale(50)[as.numeric(cut(c(E(g)$weight, -1 * max_range), breaks=50))], -1)
V(g)$color <- "white"

lab.locs <- radian.rescale(x=1:vcount(g), direction=-1, start=0)

pdf(paste0(out_file, ".pdf"), width=6, height=5) 
par(mar=c(1,2,1,5))
plot(g, layout=layout.circle(g), vertex.size=15, vertex.label.cex=0.45, vertex.label.color="black", vertex.label.font=2, edge.width=1)
image.plot(legend.only=T, zlim=c(-1, 1), col=c_scale(50), horizontal=F, legend.shrink=0.75, legend.width=0.9, legend.mar=5.5, legend.cex=1.0, legend.lab='LogFC', axis.args=list(cex.axis=0.6, tck=-0.6))
dev.off() 

z <- cbind(edges, toptable_after[rownames(table), 'logFC'])
g <- graph.data.frame(z, directed=F)
weight <- as.numeric(z[, 3])
E(g)$weight <- weight
E(g)$color <- head(c_scale(50)[as.numeric(cut(c(E(g)$weight, max_range, -1*max_range), breaks=50))], -2)
E(g)$weight <- 1
V(g)$color <- "white"

lab.locs <- radian.rescale(x=1:vcount(g), direction=-1, start=0)

pdf(paste0(out_file, "_caiman.pdf"),  width=6, height=5) 
par(mar=c(1,2,1,5))
plot(g, layout=layout.circle(g), vertex.size=15, vertex.label.cex=0.45, vertex.label.color="black", vertex.label.font=2, edge.width=1)
image.plot(legend.only=T, zlim=c(-1, 1), col=c_scale(50), horizontal=F, legend.shrink=0.75, legend.width=0.9, legend.mar=5.5, legend.cex=1.0, legend.lab='LogFC', axis.args=list(cex.axis=0.6, tck=-0.6))
dev.off() 
