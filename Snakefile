configfile: "config.yaml"

rule all:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=[
                'meta', 'meta_tissue',
                'xprs_qsmooth', 'xprs_snail',
                'xprs_qsmooth_symbol', 'xprs_snail_symbol',
                'xprs_quantile', 'xprs_rle', 'xprs_tmm'
            ]
        ),
        expand("{out_dir}/ENCODE/{dataset}/{metrics}.html",
            out_dir=config['out_dir'],
            dataset=['encode'],
            metrics=['auroc', 'auprc']
        ),
        expand("{out_dir}/{dataset}/spearman_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['GTEx/gtex', 'ENCODE/encode'],
            data=['qsmooth', 'snail']
        ),
        expand("{out_dir}/{dataset}/spearman_pvalue_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['GTEx/gtex', 'ENCODE/encode'],
            data=['qsmooth', 'snail']
        ),
        expand("{out_dir}/{dataset}/pvalue_distribution.png",
            out_dir=config['out_dir'],
            dataset=['GTEx/gtex', 'ENCODE/encode'],
        ),
        expand("{out_dir}/ENCODE/encode/lioness_{data}.pdf",
            out_dir=config['out_dir'],
            data=['qsmooth', 'snail']
        ),
        expand("{out_dir}/ENCODE/comparison/{data}/{metric}.html",
            out_dir=config['out_dir'],
            data=[
                'random_genes',
                'tissue_exclusive_genes',
            ],
            metric=['auprc', 'auroc']
        ),
        expand("{dataset_dir}/{dataset}/gene_prediction_kegg_{norm}.png", 
            dataset_dir=config['out_dir'],
            dataset=['GTEx/gtex'],
            norm=[
                'qsmooth', 'snail'
            ]
        ),
        expand("{dataset_dir}/{dataset}/gene_prediction_kegg_{norm}.png", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
            norm=[
                'validation', 'qsmooth', 'snail'
            ]
        ),
        expand("{dataset_dir}/{dataset}/annotation_{metric}_curve.html", 
            dataset_dir=config['out_dir'],
            dataset=['GTEx/gtex'],
            metric=['roc', 'prc']
        ),
        expand("{dataset_dir}/{dataset}/hubscore.png", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
        ),
        expand("{dataset_dir}/{dataset}/limma_analysis.png", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
        ),
        #expand("{out_dir}/tissue_distribution_{norm}.png",
        #    out_dir=config['out_dir'],
        #    norm=["SNAIL", "RLE", "Qsmooth", "COUNT", #"TMM"]
        #),
        #expand("{out_dir}/{data}_proportion_scatter.png",
        #    out_dir=config['out_dir'],
        #    data=['GTEx', 'ENCODE']
        #)

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
            data=['xprs_qsmooth', 'xprs_quantile', 'xprs_rle', 'xprs_tmm']
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
        expand("{out_dir}/ENCODE/{dataset}/{metrics}.html",
            out_dir=config['out_dir'],
            dataset=['encode'],
            metrics=['auroc', 'auprc']
        ),
        expand("{out_dir}/{dataset}/spearman_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['GTEx/gtex', 'ENCODE/encode'],
            data=['qsmooth', 'snail']
        ),
        expand("{out_dir}/{dataset}/spearman_pvalue_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['GTEx/gtex', 'ENCODE/encode'],
            data=['qsmooth', 'snail']
        ),
        expand("{out_dir}/{dataset}/pvalue_distribution.png",
            out_dir=config['out_dir'],
            dataset=['GTEx/gtex', 'ENCODE/encode'],
        ),
        expand("{dataset_dir}/ENCODE/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=[
                'tissue_exclusive_xprs_validation', 
                'tissue_exclusive_xprs_qsmooth',
                'tissue_exclusive_xprs_snail',
            ]
        )
    run:
        shell(
            f"python3 scripts/evaluation.py config.yaml GTEx"
        ),
        shell(
            f"python3 scripts/evaluation.py config.yaml ENCODE"
        ),
        #shell(
        #    f"python3 scripts/evaluation.py config.yaml #ENCODE cutoff"
        #)

rule lioness:
    input:
        expand("{dataset_dir}/ENCODE/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=[
                'tissue_exclusive_xprs_validation_symbol', 
                'tissue_exclusive_xprs_qsmooth_symbol',
                'tissue_exclusive_xprs_snail_symbol',
            ]
        ),
    output:
        expand("{out_dir}/ENCODE/encode/lioness_{data}.pdf",
            out_dir=config['out_dir'],
            data=['qsmooth', 'snail']
        ),
        expand("{dataset_dir}/{dataset}/lioness_{norm}_limma.tsv", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
            norm=['snail', 'qsmooth']
        )
    shell:
        f"Rscript --vanilla scripts/lioness.R"
        f" {config['datasets_dir']}/ENCODE/tissue_exclusive_xprs_validation_symbol.tsv"
        f" {config['datasets_dir']}/ENCODE/tissue_exclusive_xprs_qsmooth_symbol.tsv"
        f" {config['datasets_dir']}/ENCODE/tissue_exclusive_xprs_snail_symbol.tsv"
        f" {config['out_dir']}/ENCODE/encode/lioness_qsmooth"
        f" {config['out_dir']}/ENCODE/encode/lioness_snail"

rule comparison:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            data=[
                'meta_tissue', 'xprs_validation', 'xprs_rle', 'xprs_tmm', 'xprs_qsmooth'
            ]
        ),
    output:
        expand("{out_dir}/ENCODE/comparison/{data}/{metric}.html",
            out_dir=config['out_dir'],
            data=[
                'random_genes',
                'tissue_exclusive_genes',
            ],
            metric=['auprc', 'auroc']
        ),
    shell:
        "python3 scripts/comparison.py config.yaml validation rle tmm qsmooth snail"

rule distribution:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx'],
            data=[
                'meta_tissue', 'xprs_count', 'xprs_rle', 'xprs_tmm', 'xprs_qsmooth', 'xprs_snail'
            ]
        ),
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            data=['xprs_count']
        )
    output:
        expand("{out_dir}/tissue_distribution_{norm}.png",
            out_dir=config['out_dir'],
            norm=["SNAIL", "RLE", "Qsmooth", "COUNT", "TMM"]
        ),
        expand("{out_dir}/{data}_proportion_scatter.png",
            out_dir=config['out_dir'],
            data=['GTEx', 'ENCODE']
        )
    shell:
        "python3 scripts/distribution_analysis.py config.yaml"

rule convert_gene_name:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx'],
            data=[
                'xprs_count', 'xprs_qsmooth', 'xprs_snail'
            ]
        ),
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            data=['xprs_count', 'xprs_validation', 'xprs_qsmooth', 'xprs_snail']
        )
    output:
        expand("{dataset_dir}/{dataset}/{data}_symbol.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx'],
            data=[
                'xprs_count', 'xprs_qsmooth', 'xprs_snail', 'tissue_exclusive_xprs_count', 'tissue_exclusive_xprs_qsmooth', 'tissue_exclusive_xprs_snail'
            ]
        ),
        expand("{dataset_dir}/{dataset}/{data}_symbol.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            data=['xprs_validation', 'xprs_qsmooth', 'xprs_snail', 'tissue_exclusive_xprs_validation', 'tissue_exclusive_xprs_qsmooth', 'tissue_exclusive_xprs_snail']
        )
    run:
        shell("python3 scripts/convert_gene_name.py config.yaml GTEx"),
        shell("python3 scripts/convert_gene_name.py config.yaml ENCODE")

rule predict_gene_function:
    input:
        expand("{dataset_dir}/{dataset}/{data}_symbol.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            data=['xprs_validation', 'xprs_qsmooth', 'xprs_snail', 'tissue_exclusive_xprs_validation', 'tissue_exclusive_xprs_qsmooth', 'tissue_exclusive_xprs_snail']
        ),
        expand("{dataset_dir}/{dataset}/{data}_symbol.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx'],
            data=['xprs_count', 'xprs_qsmooth', 'xprs_snail', 'tissue_exclusive_xprs_count', 'tissue_exclusive_xprs_qsmooth', 'tissue_exclusive_xprs_snail']
        )
    output:
        expand("{dataset_dir}/{dataset}/gene_prediction_kegg_{norm}.png", 
            dataset_dir=config['out_dir'],
            dataset=['GTEx/gtex'],
            norm=[
                'qsmooth', 'snail'
            ]
        ),
        expand("{dataset_dir}/{dataset}/gene_prediction_kegg_{norm}.png", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
            norm=[
                'validation', 'qsmooth', 'snail'
            ]
        ),
    run:
        shell("python3 scripts/predict_gene_function.py config.yaml GTEx"),
        shell("python3 scripts/predict_gene_function.py config.yaml ENCODE")

rule annotate_false_positives:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx'],
            data=[
                'xprs_count', 'xprs_qsmooth', 'xprs_snail'
            ]
        ),
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            data=['xprs_count', 'xprs_validation', 'xprs_qsmooth', 'xprs_snail']
        )
    output:
        expand("{dataset_dir}/{dataset}/annotation_{metric}_curve.html", 
            dataset_dir=config['out_dir'],
            dataset=['GTEx/gtex'],
            metric=['roc', 'prc']
        ),
    run:
        shell("python3 scripts/annotate_false_positives.py config.yaml"),

rule hubscore_limma:
    input:
        expand("{dataset_dir}/{dataset}/tissue_exclusive_xprs_{ref}_symbol.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx'],
            ref=['count', 'qsmooth', 'snail']
        ),
        expand("{dataset_dir}/{dataset}/tissue_exclusive_xprs_{ref}_symbol.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['ENCODE'],
            ref=['validation', 'qsmooth', 'snail']
        ),
        expand("{dataset_dir}/{dataset}/lioness_{norm}_limma.tsv", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
            norm=['snail', 'qsmooth']
        )
    output:
        expand("{dataset_dir}/{dataset}/hubscore.png", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
        ),
        expand("{dataset_dir}/{dataset}/limma_analysis.png", 
            dataset_dir=config['out_dir'],
            dataset=['ENCODE/encode'],
        ),
    run:
        shell("python3 scripts/compute_lioness.py config.yaml ENCODE"),
        shell("python3 scripts/compute_lioness.py config.yaml GTEx")