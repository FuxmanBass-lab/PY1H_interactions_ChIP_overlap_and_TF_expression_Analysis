# pY1H vs ChIP-seq analysis

Compares pY1H interactions with ChIP-seq data from GTRD. Tests if the overlap is higher than random.

## Scripts for pY1H vs ChIP-seq analysis
ChIP_Seq_Analysis/scripts/s0_download_gtrd.sh

Downloads GTRD data.

ChIP_Seq_Analysis/scripts/s1_generate_all_tfs_chipseq.py
Generates the list of TFs with ChIP-seq data in GTRD.

ChIP_Seq_Analysis/scripts/s2_extract_chipseq_evidence.py
Looks up ChIP-seq peaks for each TF-cytokine pair.

ChIP_Seq_Analysis/scripts/s3_randomization_and_noevidence.R
Runs the randomization test and separates the pairs with no evidence.

## How to run pY1H vs ChIP-seq analysis

Needs Linux or WSL.

```
conda env create -f environment.yml
conda activate tf_het_analysis_2026
cd analysis_2026
python ChIP_Seq_Analysis/scripts/s1_generate_all_tfs_chipseq.py
python ChIP_Seq_Analysis/scripts/s2_extract_chipseq_evidence.py
Rscript ChIP_Seq_Analysis/scripts/s3_randomization_and_noevidence.R
```

For a faster s2, copy the GTRD peak files to a local disk and point to them:

```
export GTRD_MACS2_DIR_OVERRIDE=/local/path/MACS2
```

## Scripts for TF expression among tissues and tissue-cell types analysis: tissue cell type specificity score analysis


Tabula sapien tissue specific score analysis/tissue.specificity.score.R.

## How to run TF expression among tissues and tissue-cell types analysis: tissue cell type specificity score analysis

Use R to run "Tabula sapien tissue specific score analysis/tissue.specificity.score.R".




## Scripts for TF expression among tissues and tissue-cell types analysis: jaccard index and simpson index of transcription factors


Tabula sapien jaccard and simpson index index/Jaccard.Index.and.Simpson.Index.Calculation.Rmd

## How to run TF expression among tissues and tissue-cell types analysis: jaccard index and simpson index of transcription factors

Use R studio to render "Tabula sapien jaccard and simpson index index/Jaccard.Index.and.Simpson.Index.Calculation.Rmd".

A rendered html file is in "Tabula sapien jaccard and simpson index index/Jaccard.Index.and.Simpson.Index.Calculation.html".



