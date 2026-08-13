# pY1H vs ChIP-seq analysis

Compares pY1H interactions with ChIP-seq data from GTRD. Tests if the overlap is higher than random.

## Scripts

scripts/s0_download_gtrd.sh
Downloads GTRD data.

scripts/s1_generate_all_tfs_chipseq.py
Generates the list of TFs with ChIP-seq data in GTRD.

scripts/s2_extract_chipseq_evidence.py
Looks up ChIP-seq peaks for each TF-cytokine pair.

scripts/s3_randomization_and_noevidence.R
Runs the randomization test and separates the pairs with no evidence.

## How to run

Needs Linux or WSL.

```
conda env create -f environment.yml
conda activate tf_het_analysis_2026
cd analysis_2026
python scripts/s1_generate_all_tfs_chipseq.py
python scripts/s2_extract_chipseq_evidence.py
Rscript scripts/s3_randomization_and_noevidence.R
```

For a faster s2, copy the GTRD peak files to a local disk and point to them:

```
export GTRD_MACS2_DIR_OVERRIDE=/local/path/MACS2
```
