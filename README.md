# pY1H reproducibility

Analysis code and data for the pY1H Nat Comm 2026 paper. This repo currently covers three analyses: pY1H vs. ChIP-seq overlap, TF/cytokine tissue-cell-type expression specificity, and single-cell RNA-seq processing.

## pY1H PDIs vs. ChIP-seq (GTRD) overlap analysis

Compares protein-DNA interactions (PDIs) detected by pY1H with independent ChIP-seq evidence from GTRD, and tests whether the overlap is higher than expected by chance with a degree-preserving network randomization.

### s0 - Download GTRD ChIP-seq data
Downloads the GTRD MACS2 peak files and metadata table used by the other steps.

**Input Files**
- None - pulls directly from the GTRD server.

**Scripts**
- `ChIP_Seq_Analysis/scripts/s0_download_gtrd.sh`: Downloads GTRD MACS2 ChIP-seq peak files (hg38, `.bb`) and the metadata table via `wget`.

**Output Files**
- `ChIP_Seq_Analysis/inputs/ChIP-seq.metadata.txt`: GTRD ChIP-seq experiment metadata.
- GTRD MACS2 peak files (`.bb`, hg38): not stored in this repo. Point `GTRD_MACS2_DIR_OVERRIDE` at wherever you downloaded them for s2.

### s1 - Build the TF ChIP-seq availability list
Resolves every TF Uniprot ID in the GTRD metadata to a gene symbol, so later steps know which TFs have ChIP-seq data.

**Input Files**
- `ChIP_Seq_Analysis/inputs/ChIP-seq.metadata.txt`: GTRD ChIP-seq experiment metadata (from s0).
- `ChIP_Seq_Analysis/manual_lookups/no_found_biomart.txt`: manual Uniprot ID → gene symbol lookups for TFs BioMart can't resolve.

**Scripts**
- `ChIP_Seq_Analysis/scripts/s1_generate_all_tfs_chipseq.py`: Builds the list of TFs with GTRD ChIP-seq data (gene symbol + Uniprot ID) via Ensembl BioMart, filled in with manual lookups where BioMart fails.

**Output Files**
- `ChIP_Seq_Analysis/outputs/all_tfs_chipseq.txt`: TF gene symbol ↔ Uniprot ID table for TFs with GTRD data.
- `ChIP_Seq_Analysis/outputs/unresolved_uniprot_ids.txt`: Uniprot IDs neither BioMart nor the manual lookup table could resolve.

### s2 - Extract ChIP-seq evidence for pY1H PDIs
For every TF in the pY1H list, checks its GTRD peaks against every cytokine promoter, building the full TF x cytokine reference network plus a peak count for each reported PDI.

**Input Files**
- `ChIP_Seq_Analysis/inputs/py1h PDIs for ChIP-seq comparison.xlsx`: pY1H TF-cytokine PDIs to check against ChIP-seq evidence.
- `ChIP_Seq_Analysis/inputs/ChIP-seq.metadata.txt`: GTRD ChIP-seq experiment metadata (from s0).
- `ChIP_Seq_Analysis/outputs/all_tfs_chipseq.txt`: TFs with GTRD data (from s1).
- GTRD MACS2 peak files (from s0); set `GTRD_MACS2_DIR_OVERRIDE` to their local path.

**Scripts**
- `ChIP_Seq_Analysis/scripts/s2_extract_chipseq_evidence.py`: Looks up ChIP-seq peak evidence for every TF x cytokine-promoter pair using `pybbi` to read the GTRD `.bb` peak files.

**Output Files**
- `ChIP_Seq_Analysis/outputs/py1h_macs2_peaks_full.xlsx`: Every TF x cytokine ChIP-seq peak match (the full reference network).
- `ChIP_Seq_Analysis/outputs/py1h_chipseq_evidence.xlsx`: Number of ChIP-seq peaks (n_evidence) supporting each pY1H PDI.

### s3 - Randomization test and no-evidence set
Tests whether the true pY1H/ChIP-seq overlap exceeds chance, and separates out PDIs with no supporting ChIP-seq peak.

**Input Files**
- `ChIP_Seq_Analysis/outputs/py1h_chipseq_evidence.xlsx` (from s2).
- `ChIP_Seq_Analysis/outputs/py1h_macs2_peaks_full.xlsx` (from s2).

**Scripts**
- `ChIP_Seq_Analysis/scripts/s3_randomization_and_noevidence.R`: Runs a degree-preserving edge-rewiring randomization test (100,000 iterations) of the pY1H PDI network against the ChIP-seq reference network, and separates out PDIs with no ChIP-seq evidence.

**Output Files**
- `ChIP_Seq_Analysis/outputs/py1h_randomization_results.csv`: Overlap counts for the true network and all 100,000 randomized networks.
- `ChIP_Seq_Analysis/outputs/py1h_randomization_histogram.png`: Histogram of randomized overlaps vs. the true overlap.
- `ChIP_Seq_Analysis/outputs/py1h_noevidence_TF_available.xlsx`: PDIs where the TF has GTRD data but no ChIP-seq peak was found for that specific TF-cytokine pair.

### How to run
Needs Linux or WSL.
```
conda env create -f ChIP_Seq_Analysis/environment.yml
conda activate tf_het_analysis_2026
bash ChIP_Seq_Analysis/scripts/s0_download_gtrd.sh
python ChIP_Seq_Analysis/scripts/s1_generate_all_tfs_chipseq.py
python ChIP_Seq_Analysis/scripts/s2_extract_chipseq_evidence.py
Rscript ChIP_Seq_Analysis/scripts/s3_randomization_and_noevidence.R
```
For a faster s2, copy the GTRD peak files to a local disk and point to them:
```
export GTRD_MACS2_DIR_OVERRIDE=/local/path/MACS2
```

### Software / package versions
Defined in `ChIP_Seq_Analysis/environment.yml` (conda env `tf_het_analysis_2026`):
- Python 3.10, with `pandas`, `openpyxl`, and (via pip) `pybbi`, `biomart`.
- R (`r-base`, conda-forge/bioconda) with `dplyr`, `readxl`, `writexl`, `igraph`, `reshape2`, `ggplot2`, `readr`.
- `s0_download_gtrd.sh` additionally needs `wget`.

## TF and cytokine tissue/cell-type expression specificity (Tabula Sapiens)

Computes a per-gene tissue/cell-type expression specificity score (an entropy-style score over CPM) for TFs and cytokines from Tabula Sapiens bulk CPM tables.

### Input Files
- `Tabula sapien tissue specific score analysis/Transcription.factors.Tissue_Celltype.CPM.fixed.tsv` and `...fixed.log2.transformed.tsv`: TF CPM expression per tissue/cell type.
- `Cytokine.Tissue_Celltype.CPM.fixed.tsv`: cytokine CPM expression per tissue/cell type, read by the second half of the script. Not included in this repo - place it alongside the TF table before running.

### Scripts
- `Tabula sapien tissue specific score analysis/tissue.specificity.score.R`: Computes the tissue/cell-type specificity score and the number of tissues/cell types each gene is expressed in, run once for TFs and once for cytokines.

### Output Files
- `Tabula sapien tissue specific score analysis/Tabula.sapiens.TF.tissue.celltype.specificity.score.tsv`
- `Tabula sapien tissue specific score analysis/Tabula.sapiens.Cytokine.tissue.celltype.specificity.score.tsv`

### Software / package versions
- Base R only (`read.table`/`write.table`), no additional packages required. `R_system_info.html` in this folder records R 4.4.2 Patched (2025-02-15 r87725) as the R version used.

## Single-cell RNA-seq analysis (Tabula Sapiens)

Processes Tabula Sapiens 10x single-cell RNA-seq counts per sample: QC filtering, then per-tissue cluster-marker heatmaps and cell-type composition plots.

### Input Files
- `Counts/<sample>/`: 10x Genomics count matrices, one folder per sample (read with `Read10X`). Not included in this repo.
- `metadata.tsv`: per-sample metadata (`donor`, `organ_tissue`, `gender`). Not included in this repo.
- Raw input data is available from Tabula Sapiens: https://doi.org/10.6084/m9.figshare.27921984.v1

### Scripts
- `single-cell-rna-seq-analysis/single.cell.processing.script.R`: Loads and merges all samples into one Seurat object; QC-filters cells (`nFeature_RNA` 500-7500, `nCount_RNA` < 10000, `percent.mt` < 25); then, per pre-clustered `.rds` file in `Results/`, normalizes, finds cluster markers, and plots a marker-gene heatmap and cell-type composition bar chart.

### Output Files
All in `single-cell-rna-seq-analysis/Results/`:
- `QC.pdf`: nFeature_RNA / nCount_RNA / percent.mt violin plots per sample.
- `Unfiltered.rds` / `Filtered.rds`: merged Seurat object before/after QC filtering.
- `<Tissue>.Celltype.Heatmap.pdf`: top marker-gene heatmap per tissue.
- `<Tissue>.Celltype.composition.pdf`: cell-type composition bar chart per tissue.

### Software / package versions
- R with the Seurat package (`Read10X`, `CreateSeuratObject`, `PercentageFeatureSet`, `VlnPlot`, `NormalizeData`, `FindAllMarkers`, `ScaleData`, `DoHeatmap`) and `dplyr`, `ggplot2`.
