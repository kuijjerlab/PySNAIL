configfile: "config.yaml"

rule all:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=[
                'meta', 'meta_tissue',
                'xprs_qsmooth', 'xprs_snail',
                'xprs_quantile', 'xprs_deseq', 'xprs_uq', 'xprs_tmm'
            ]
        ),
        expand("{out_dir}/{dataset}_spearman_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            data=['qsmooth', 'snail']
        ),
        expand("{out_dir}/{dataset}_{label}_{metrics}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            label=['comparison'],
            metrics=['auroc', 'auprc']
        ),
        expand("{out_dir}/encode_lioness_{data}.pdf",
            out_dir=config['out_dir'],
            data=['qsmooth', 'snail']
        ),
        expand("{out_dir}/{data}_comparison_{metric}.html",
            out_dir=config['out_dir'],
            data=[
                'tissue-exclusive_genes',
                'expressed_proportion_[0%,_25%)',
                'expressed_proportion_[25%,_50%)',
                'expressed_proportion_[50%,_75%)',
                'expressed_proportion_[75%,_100%)'
            ],
            metric=['auprc', 'auroc']
        )

rule download_encode:
    output: expand("{dataset_dir}/ENCODE/{sample}/{sample}_meta.tsv", dataset_dir=config['datasets_dir'], sample=config['encode_tissues'])
    shell: f"Rscript --vanilla scripts/download_encode.R {config['datasets_dir']}/ENCODE/"

rule download_gtex:
    output:
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    shell: f"bash scripts/download_gtex.sh {config['datasets_dir']}"

rule process_gtex:
    input:
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        f"{config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    output:
        expand("{dataset_dir}/GTEx/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=['meta', 'meta_tissue', 'xprs_count']
        ),
    shell:
        f"python3 scripts/process_gtex.py"
        f" --xprs {config['datasets_dir']}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
        f" --sample {config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
        f" --subject {config['datasets_dir']}/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
        f" --config config.yaml"

rule process_encode:
    input:
        expand("{dataset_dir}/ENCODE/{sample}/{sample}_meta.tsv",
            dataset_dir=config['datasets_dir'],
            sample=config['encode_tissues']
        )
    output:
        expand("{dataset_dir}/ENCODE/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=['meta', 'meta_tissue', 'xprs_validation', 'xprs_count']
        )
    shell:
        f"python3 scripts/process_encode.py"
        f" --dataset {config['datasets_dir']}/ENCODE/"
        f" --config config.yaml"

rule normalization:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['meta_tissue', 'xprs_count']
        ),
    output:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_qsmooth', 'xprs_quantile', 'xprs_deseq', 'xprs_uq', 'xprs_tmm']
        ),
    run:
        shell(
            f"Rscript --vanilla"
            f" scripts/normalize_expression.R"
            f" {config['datasets_dir']}/GTEx/xprs_count.tsv"
            f" {config['datasets_dir']}/GTEx/meta_tissue.tsv"
        ),
        shell(
            f"Rscript --vanilla"
            f" scripts/normalize_expression.R"
            f" {config['datasets_dir']}/ENCODE/xprs_count.tsv"
            f" {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
        )

rule snail_correction:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_count']
        )
    output:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_snail']
        )
    run:
        shell(
            f"pysnail "
            f" --groups {config['datasets_dir']}/GTEx/meta_tissue.tsv"
            f" --outdir {config['datasets_dir']}/GTEx/snail"
            f" {config['datasets_dir']}/GTEx/xprs_count.tsv"
        ),
        shell(
            f"mv {config['datasets_dir']}/GTEx/snail/xprs_norm.tsv"
            f" {config['datasets_dir']}/GTEx/xprs_snail.tsv"
        ),
        shell(
            f"pysnail "
            f" --groups {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
            f" --outdir {config['datasets_dir']}/ENCODE/snail"
            f" {config['datasets_dir']}/ENCODE/xprs_count.tsv"
        ),
        shell(
            f"mv {config['datasets_dir']}/ENCODE/snail/xprs_norm.tsv"
            f" {config['datasets_dir']}/ENCODE/xprs_snail.tsv"
        )

rule evaluation:
    input: 
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=[
                'meta_tissue',
                'xprs_qsmooth', 'xprs_snail',
            ]
        ),
        f"{config['datasets_dir']}/GTEx/xprs_count.tsv",
        f"{config['datasets_dir']}/ENCODE/xprs_validation.tsv"
    output:
        expand("{out_dir}/{dataset}_comparison_{metrics}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            metrics=['auroc', 'auprc']
        ),
        expand("{out_dir}/{dataset}_spearman_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            data=['qsmooth', 'snail']
        ),
        expand("{dataset_dir}/ENCODE/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=[
                'lioness_input_ref_qsmooth_snail', 
                'lioness_input_qsmooth',
                'lioness_input_snail',
            ]
        )
    run:
        shell(
            f"python3 scripts/evaluation.py"
            f" -r {config['datasets_dir']}/ENCODE/xprs_validation.tsv"
            f" -t {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
            f" -d ENCODE -x validation -y Qsmooth -z SNAIL"
            f" -a {config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv"
            f" -b {config['datasets_dir']}/ENCODE/xprs_snail.tsv"
            f" -c config.yaml"
        ),
        shell(
            f"python3 scripts/evaluation.py"
            f" -r {config['datasets_dir']}/GTEx/xprs_count.tsv"
            f" -t {config['datasets_dir']}/GTEx/meta_tissue.tsv"
            f" -d GTEx -x count -y Qsmooth -z SNAIL"
            f" -a {config['datasets_dir']}/GTEx/xprs_qsmooth.tsv"
            f" -b {config['datasets_dir']}/GTEx/xprs_snail.tsv"
            f" -c config.yaml"
        )

rule lioness:
    input:
        expand("{dataset_dir}/ENCODE/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=[
                'lioness_input_ref_qsmooth_snail', 
                'lioness_input_qsmooth',
                'lioness_input_snail',
            ]
        ),
    output:
        expand("{out_dir}/encode_lioness_{data}.pdf",
            out_dir=config['out_dir'],
            data=['qsmooth', 'snail']
        )
    shell:
        f"Rscript --vanilla scripts/lioness.R"
        f" {config['datasets_dir']}/ENCODE/lioness_input_ref_qsmooth_snail.tsv"
        f" {config['datasets_dir']}/ENCODE/lioness_input_qsmooth.tsv"
        f" {config['datasets_dir']}/ENCODE/lioness_input_snail.tsv"
        f" {config['out_dir']}/encode_lioness_qsmooth"
        f" {config['out_dir']}/encode_lioness_snail"

rule comparison:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            data=[
                'meta_tissue', 'xprs_validation', 'xprs_deseq', 'xprs_tmm', 'xprs_uq',
                'xprs_qsmooth', 'xprs_snail',
            ]
        ),
    output:
        expand("{out_dir}/{data}_comparison_{metric}.html",
            out_dir=config['out_dir'],
            data=[
                'tissue-exclusive_genes',
                'expressed_proportion_[0%,_25%)',
                'expressed_proportion_[25%,_50%)',
                'expressed_proportion_[50%,_75%)',
                'expressed_proportion_[75%,_100%)'
            ],
            metric=['auprc', 'auroc']
        )
    shell:
        "python3 scripts/comparison.py config.yaml validation deseq tmm uq qsmooth snail"
