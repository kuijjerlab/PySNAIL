configfile: "config.yaml"

rule all:
    input:
        f"{config['datasets_dir']}/GTEx/filtered_samples_meta.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_meta_gender.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_xprs_qsmooth.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_xprs_caiman.tsv",
        f"{config['datasets_dir']}/ENCODE/meta.tsv",
        f"{config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv",
        f"{config['datasets_dir']}/ENCODE/xprs_caiman.tsv",
        f"{config['datasets_dir']}/ENCODE/xprs_validation.tsv",
        f"{config['out_dir']}/gtex_number_of_non_expressed_genes.html",
        f"{config['out_dir']}/gtex_xprs_distribution.html",
        f"{config['out_dir']}/encode_spikeins_expression.html",
        f"{config['out_dir']}/encode_number_of_non_expressed_genes.html",
        f"{config['out_dir']}/encode_xprs_distribution.html"

rule download_encode:
    output: expand("{dataset_dir}/ENCODE/{sample}/{sample}_meta.tsv", dataset_dir=config['datasets_dir'], sample=config['encode_tissues'])
    shell: f"Rscript --vanilla scripts/download_encode.R {config['datasets_dir']}/ENCODE/"

rule download_gtex:
    output:
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    shell: f"bash scripts/download_gtex.sh {config['datasets_dir']}/GTEx"

rule process_gtex:
    input:
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    output:
        f"{config['datasets_dir']}/GTEx/filtered_samples_meta.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_meta_gender.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_meta_tissue.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_xprs_count.tsv",
        f"{config['out_dir']}/gtex_number_of_non_expressed_genes.html",
        f"{config['out_dir']}/gtex_xprs_distribution.html"
    shell:
        f"python3 scripts/process_gtex.py"
        f" --xprs {config['datasets_dir']}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
        f" --sample {config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
        f" --subject {config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
        f" --config config.yaml"

rule process_encode:
    input:
        expand("{dataset_dir}/ENCODE/{sample}/{sample}_meta.tsv", dataset_dir=config['datasets_dir'], sample=config['encode_tissues'])
    output:
        f"{config['datasets_dir']}/ENCODE/meta.tsv",
        f"{config['datasets_dir']}/ENCODE/meta_tissue.tsv",
        f"{config['datasets_dir']}/ENCODE/xprs_count.tsv",
        f"{config['datasets_dir']}/ENCODE/xprs_validation.tsv",
        f"{config['out_dir']}/encode_spikeins_expression.html",
        f"{config['out_dir']}/encode_number_of_non_expressed_genes.html",
        f"{config['out_dir']}/encode_xprs_distribution.html"
    shell:
        f"python3 scripts/process_encode.py"
        f" --dataset {config['datasets_dir']}/ENCODE/"
        f" --config config.yaml"

rule qsmooth_normalization_gtex:
    input:
        f"{config['datasets_dir']}/GTEx/filtered_samples_xprs_count.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_meta_tissue.tsv"
    output:
        f"{config['datasets_dir']}/GTEx/filtered_samples_xprs_qsmooth.tsv",
    shell:
        f"Rscript --vanilla scripts/qsmooth_normalization.R {config['datasets_dir']}/GTEx/filtered_samples_xprs_count.tsv {config['datasets_dir']}/GTEx/filtered_samples_meta_tissue.tsv"

rule qsmooth_normalization_encode:
    input:
        f"{config['datasets_dir']}/ENCODE/xprs_count.tsv",
        f"{config['datasets_dir']}/ENCODE/meta_tissue.tsv"
    output:
        f"{config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv"
    shell:
        f"Rscript --vanilla scripts/qsmooth_normalization.R  {config['datasets_dir']}/ENCODE/xprs_count.tsv {config['datasets_dir']}/ENCODE/meta_tissue.tsv"

rule caiman_correction_gtex:
    input:
        f"{config['datasets_dir']}/GTEx/filtered_samples_xprs_qsmooth.tsv",
        f"{config['datasets_dir']}/GTEx/filtered_samples_meta_tissue.tsv"
    output:
        f"{config['datasets_dir']}/GTEx/filtered_samples_xprs_caiman.tsv",
    run:
        shell(f"python3 -m caiman {config['datasets_dir']}/GTEx/filtered_samples_xprs_qsmooth.tsv --groups {config['datasets_dir']}/GTEx/filtered_samples_meta_tissue.tsv --dist --save_model --outdir {config['datasets_dir']}/GTEx/caiman/ --verbose")
        shell(f"mv {config['datasets_dir']}/GTEx/caiman/xprs_caiman.tsv {config['datasets_dir']}/GTEx/filtered_samples_xprs_caiman.tsv")

rule caiman_correction_encode:
    input:
        f"{config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv",
        f"{config['datasets_dir']}/ENCODE/meta_tissue.tsv"
    output:
        f"{config['datasets_dir']}/ENCODE/xprs_caiman.tsv"
    run:
        shell(f"python3 -m caiman {config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv --groups {config['datasets_dir']}/ENCODE/meta_tissue.tsv --dist --save_model --outdir {config['datasets_dir']}/ENCODE/caiman/ --verbose")
        shell(f"mv {config['datasets_dir']}/ENCODE/caiman/xprs_caiman.tsv {config['datasets_dir']}/ENCODE/xprs_caiman.tsv")