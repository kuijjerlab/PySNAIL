library(ENCODExplorer)

encode_df <- get_encode_df()

args <- commandArgs(trailingOnly=TRUE) 
outdir <- args[1]

dir.create(outdir, recursive=TRUE)
setwd(outdir)

tissues <- c(
    "heart", "liver", "forebrain", "hindbrain", "midbrain", "limb", 
    "embryonic facial prominence", "neural tube", "kidney", "lung", 
    "stomach", "intestine"
)

for(tissue in tissues){
    tissue_underscore = gsub(" ", "_", tissue)
    dir.create(tissue_underscore)
    setwd(tissue_underscore)
    query_results <- queryEncode(
        organism = "Mus musculus", 
        biosample_name = tissue,
        assay="polyA plus RNA-seq",
        file_format = "tsv",
        fixed = TRUE
    )

    filter <- query_results$output_type == "gene quantifications"
    gene_xprs <- query_results[filter]
    downloadEncode(gene_xprs)

    write.table(
        gene_xprs, paste0(tissue_underscore, "_meta.tsv"), 
        sep="\t", 
        row.names=FALSE, 
        quote=FALSE
    )
    setwd("../")
}