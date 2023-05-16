pysnail --groups ./manuscript_analysis_20220330/datasets/GTEx/meta_tissue.tsv --outdir ./manuscript_analysis_20220330/datasets/GTEx/snail_010 ./manuscript_analysis_20220330/datasets/GTEx/xprs_count.tsv -c 0.10
pysnail --groups ./manuscript_analysis_20220330/datasets/GTEx/meta_tissue.tsv --outdir ./manuscript_analysis_20220330/datasets/GTEx/snail_020 ./manuscript_analysis_20220330/datasets/GTEx/xprs_count.tsv -c 0.20
pysnail --groups ./manuscript_analysis_20220330/datasets/GTEx/meta_tissue.tsv --outdir ./manuscript_analysis_20220330/datasets/GTEx/snail_025 ./manuscript_analysis_20220330/datasets/GTEx/xprs_count.tsv -c 0.25

pysnail --groups ./manuscript_analysis_20220330/datasets/ENCODE/meta_tissue.tsv --outdir ./manuscript_analysis_20220330/datasets/ENCODE/snail_010 ./manuscript_analysis_20220330/datasets/ENCODE/xprs_count.tsv -c 0.10
pysnail --groups ./manuscript_analysis_20220330/datasets/ENCODE/meta_tissue.tsv --outdir ./manuscript_analysis_20220330/datasets/ENCODE/snail_020 ./manuscript_analysis_20220330/datasets/ENCODE/xprs_count.tsv -c 0.20
pysnail --groups ./manuscript_analysis_20220330/datasets/ENCODE/meta_tissue.tsv --outdir ./manuscript_analysis_20220330/datasets/ENCODE/snail_025 ./manuscript_analysis_20220330/datasets/ENCODE/xprs_count.tsv -c 0.25

mkdir -p manuscript_analysis_20220330_010/results
mkdir -p manuscript_analysis_20220330_010/datasets/ENCODE
mkdir -p manuscript_analysis_20220330_010/ENCODE
python3 scripts/evaluation.py -r ./manuscript_analysis_20220330/datasets/ENCODE/xprs_validation.tsv -t ./manuscript_analysis_20220330/datasets/ENCODE/meta_tissue.tsv -d ENCODE -x validation -y Qsmooth -z SNAIL_010 -a ./manuscript_analysis_20220330/datasets/ENCODE/xprs_qsmooth.tsv -b ./manuscript_analysis_20220330/datasets/ENCODE/snail_010/xprs_norm.tsv -c config.yaml

mkdir -p manuscript_analysis_20220330_020/results
mkdir -p manuscript_analysis_20220330_020/datasets/ENCODE
mkdir -p manuscript_analysis_20220330_020/ENCODE
python3 scripts/evaluation.py -r ./manuscript_analysis_20220330/datasets/ENCODE/xprs_validation.tsv -t ./manuscript_analysis_20220330/datasets/ENCODE/meta_tissue.tsv -d ENCODE -x validation -y Qsmooth -z SNAIL_020 -a ./manuscript_analysis_20220330/datasets/ENCODE/xprs_qsmooth.tsv -b ./manuscript_analysis_20220330/datasets/ENCODE/snail_020/xprs_norm.tsv -c config.yaml

mkdir -p manuscript_analysis_20220330_025/results
mkdir -p manuscript_analysis_20220330_025/datasets/ENCODE
mkdir -p manuscript_analysis_20220330_025/ENCODE
python3 scripts/evaluation.py -r ./manuscript_analysis_20220330/datasets/ENCODE/xprs_validation.tsv -t ./manuscript_analysis_20220330/datasets/ENCODE/meta_tissue.tsv -d ENCODE -x validation -y Qsmooth -z SNAIL_025 -a ./manuscript_analysis_20220330/datasets/ENCODE/xprs_qsmooth.tsv -b ./manuscript_analysis_20220330/datasets/ENCODE/snail_025/xprs_norm.tsv -c config.yaml
