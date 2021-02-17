configfile: "config.yaml"

rule all:
    input: 
        expand("datasets/ENCODE/{sample}/{sample}_meta.tsv", sample=config["encode_tissues"], num=['1', '2']),
        "datasets/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
        "datasets/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 
        "datasets/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"

rule download_encode:
    output: expand("datasets/ENCODE/{sample}/{sample}_meta.tsv", sample=config["encode_tissues"], num=['1', '2'])
    shell: "Rscript scripts/download_encode.R"

rule download_gtex:
    output:
        "datasets/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
        "datasets/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        "datasets/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    shell: "bash scripts/download_gtex.sh"
