configfile: "config.yaml"

rule all:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=[
                'meta', 'meta_tissue',
                'xprs_qsmooth', 'xprs_qsmooth_round',
                'xprs_quantile', 'xprs_quantile_round',
                'xprs_deseq', 'xprs_uq', 'xprs_tmm'
            ]
        ),
        expand("{out_dir}/{dataset}_spearman_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            data=['qsmooth', 'qsmooth_round', 'quantile', 'quantile_round']
        ),
        expand("{out_dir}/{dataset}_{data}_{metrics}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            data=['qsmooth', 'quantile'],
            metrics=['auroc', 'auprc']
        ),
        expand("{out_dir}/encode_lioness_{data}.pdf",
            out_dir=config['out_dir'],
            data=['qsmooth', 'qsmooth_round']
        ),
        #f"{config['out_dir']}/cpu_time_usage.html",
        #f"{config['out_dir']}/memory_usage.html",
        #f"{config['out_dir']}/elapsed_time.html"

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

rule qsmooth_normalization:
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
            data=['xprs_qsmooth', 'xprs_quantile']
        ),
    run:
        shell(
            f"Rscript --vanilla"
            f" scripts/qsmooth_normalization.R"
            f" {config['datasets_dir']}/GTEx/xprs_count.tsv"
            f" {config['datasets_dir']}/GTEx/meta_tissue.tsv"
        ),
        shell(
            f"Rscript --vanilla"
            f" scripts/qsmooth_normalization.R"
            f" {config['datasets_dir']}/ENCODE/xprs_count.tsv"
            f" {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
        )

rule deseq_normalization:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_count']
        ),
    output:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_uq', 'xprs_tmm', 'xprs_deseq']
        ),
    run:
        shell(
            f"Rscript --vanilla"
            f" scripts/other_normalization.R"
            f" {config['datasets_dir']}/GTEx/xprs_count.tsv"
        ),
        shell(
            f"Rscript --vanilla"
            f" scripts/other_normalization.R"
            f" {config['datasets_dir']}/ENCODE/xprs_count.tsv"
        )

rule round_correction:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_qsmooth', 'xprs_quantile']
        ),
    output:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_qsmooth_round', 'xprs_quantile_round']
        ),
    run:
        shell(
            f"python3 scripts/round_correction.py  {config['datasets_dir']}/GTEx/xprs_qsmooth.tsv {config['datasets_dir']}/GTEx/xprs_qsmooth_round.tsv",
        ),
        shell(
            f"python3 scripts/round_correction.py {config['datasets_dir']}/GTEx/xprs_quantile.tsv {config['datasets_dir']}/GTEx/xprs_quantile_round.tsv"
        ),
        shell(
            f"python3 scripts/round_correction.py {config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv {config['datasets_dir']}/ENCODE/xprs_qsmooth_round.tsv",
        ),
        shell(
            f"python3 scripts/round_correction.py {config['datasets_dir']}/ENCODE/xprs_quantile.tsv {config['datasets_dir']}/ENCODE/xprs_quantile_round.tsv",
        ),

"""rule caiman_correction:
    input:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_qsmooth', 'xprs_quantile']
        ),
    output:
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=['xprs_qsmooth_caiman', 'xprs_quantile_caiman']
        ),
    run:
        shell(
            f"python3 -m caiman "
            f" {config['datasets_dir']}/GTEx/xprs_qsmooth.tsv"
            f" --groups {config['datasets_dir']}/GTEx/meta_tissue.tsv"
            f" --outdir {config['datasets_dir']}/GTEx/caiman/qsmooth"
            f" --dist --save_model  --verbose"
        ),
        shell(
            f"mv {config['datasets_dir']}/GTEx/caiman/qsmooth/xprs_caiman.tsv"
            f" {config['datasets_dir']}/GTEx/xprs_qsmooth_caiman.tsv"
        ),
        shell(
            f"python3 -m caiman "
            f" {config['datasets_dir']}/GTEx/xprs_qsmooth.tsv"
            f" --groups {config['datasets_dir']}/GTEx/meta_tissue.tsv"
            f" --outdir {config['datasets_dir']}/GTEx/caiman/quantile"
            f" --dist --save_model  --verbose"
        ),
        shell(
            f"mv {config['datasets_dir']}/GTEx/caiman/quantile/xprs_caiman.tsv"
            f" {config['datasets_dir']}/GTEx/xprs_quantile_caiman.tsv"
        ),
        shell(
            f"python3 -m caiman "
            f" {config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv"
            f" --groups {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
            f" --outdir {config['datasets_dir']}/ENCODE/caiman/qsmooth"
            f" --dist --save_model  --verbose"
        ),
        shell(
            f"mv {config['datasets_dir']}/ENCODE/caiman/qsmooth/xprs_caiman.tsv"
            f" {config['datasets_dir']}/ENCODE/xprs_qsmooth_caiman.tsv"
        ),
        shell(
            f"python3 -m caiman "
            f" {config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv"
            f" --groups {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
            f" --outdir {config['datasets_dir']}/ENCODE/caiman/quantile"
            f" --dist --save_model  --verbose"
        ),
        shell(
            f"mv {config['datasets_dir']}/ENCODE/caiman/quantile/xprs_caiman.tsv"
            f" {config['datasets_dir']}/ENCODE/xprs_quantile_caiman.tsv"
        ),"""

rule evaluation:
    input: 
        expand("{dataset_dir}/{dataset}/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            dataset=['GTEx', 'ENCODE'],
            data=[
                'meta_tissue',
                'xprs_qsmooth', 'xprs_qsmooth_round',
                'xprs_quantile', 'xprs_quantile_round',
            ]
        ),
        f"{config['datasets_dir']}/GTEx/xprs_count.tsv",
        f"{config['datasets_dir']}/ENCODE/xprs_validation.tsv"
    output:
        expand("{out_dir}/{dataset}_{data}_{metrics}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            data=['qsmooth', 'quantile'],
            metrics=['auroc', 'auprc']
        ),
        expand("{out_dir}/{dataset}_spearman_heatmap_{data}.html",
            out_dir=config['out_dir'],
            dataset=['gtex', 'encode'],
            data=['qsmooth', 'qsmooth_round', 'quantile', 'quantile_round']
        ),
        expand("{dataset_dir}/ENCODE/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=[
                'lioness_input_ref_qsmooth', 
                'lioness_input_before_qsmooth',
                'lioness_input_after_qsmooth',
            ]
        ),
    run:
        shell(
            f"python3 scripts/evaluation.py"
            f" -r {config['datasets_dir']}/ENCODE/xprs_validation.tsv"
            f" -t {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
            f" -d ENCODE -x validation -y qsmooth"
            f" -b {config['datasets_dir']}/ENCODE/xprs_qsmooth.tsv"
            f" -a {config['datasets_dir']}/ENCODE/xprs_qsmooth_round.tsv"
            f" -c config.yaml"
        ),
        shell(
            f"python3 scripts/evaluation.py"
            f" -r {config['datasets_dir']}/ENCODE/xprs_validation.tsv"
            f" -t {config['datasets_dir']}/ENCODE/meta_tissue.tsv"
            f" -d ENCODE -x validation -y quantile"
            f" -b {config['datasets_dir']}/ENCODE/xprs_quantile.tsv"
            f" -a {config['datasets_dir']}/ENCODE/xprs_quantile_round.tsv"
            f" -c config.yaml"
        ),
        shell(
            f"python3 scripts/evaluation.py"
            f" -r {config['datasets_dir']}/GTEx/xprs_count.tsv"
            f" -t {config['datasets_dir']}/GTEx/meta_tissue.tsv"
            f" -d GTEx -x count -y qsmooth"
            f" -b {config['datasets_dir']}/GTEx/xprs_qsmooth.tsv"
            f" -a {config['datasets_dir']}/GTEx/xprs_qsmooth_round.tsv"
            f" -c config.yaml"
        ),
        shell(
            f"python3 scripts/evaluation.py"
            f" -r {config['datasets_dir']}/GTEx/xprs_count.tsv"
            f" -t {config['datasets_dir']}/GTEx/meta_tissue.tsv"
            f" -d GTEx -x count -y quantile"
            f" -b {config['datasets_dir']}/GTEx/xprs_quantile.tsv"
            f" -a {config['datasets_dir']}/GTEx/xprs_quantile_round.tsv"
            f" -c config.yaml"
        )

rule lioness:
    input:
        expand("{dataset_dir}/ENCODE/{data}.tsv", 
            dataset_dir=config['datasets_dir'],
            data=[
                'lioness_input_ref_qsmooth', 
                'lioness_input_before_qsmooth',
                'lioness_input_after_qsmooth',
            ]
        ),
    output:
        expand("{out_dir}/encode_lioness_{data}.pdf",
            out_dir=config['out_dir'],
            data=['qsmooth', 'qsmooth_round']
        )
    shell:
        f"Rscript --vanilla scripts/lioness.R"
        f" {config['datasets_dir']}/ENCODE/lioness_input_ref_qsmooth.tsv"
        f" {config['datasets_dir']}/ENCODE/lioness_input_before_qsmooth.tsv"
        f" {config['datasets_dir']}/ENCODE/lioness_input_after_qsmooth.tsv"
        f" {config['out_dir']}/encode_lioness_qsmooth"

"""rule memory_time_usage:
    input:
        f"{config['datasets_dir']}/ENCODE/xprs_validation.tsv"
    output:
        f"{config['out_dir']}/cpu_time_usage.html",
        f"{config['out_dir']}/memory_usage.html",
        f"{config['out_dir']}/elapsed_time.html"
    run:
        shell(f"rm -rf ./tmp/correction.log")
        shell(
            f"python3 scripts/memory_time_usage.py" 
            f" -x ./manuscript_analysis/datasets/ENCODE/xprs_validation.tsv"
            f" -c config.yaml"
        )"""