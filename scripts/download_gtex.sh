#!/usr/bin/env bash

OUTDIR="$1"

mkdir -p "$OUTDIR/GTEx/"
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz -O "$OUTDIR/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
gzip -d "$OUTDIR/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"

wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt -O "$OUTDIR/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt -O "$OUTDIR/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"