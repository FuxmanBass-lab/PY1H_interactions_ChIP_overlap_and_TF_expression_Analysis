"""
Looks up ChIP-seq peak evidence for each TF-cytokine pair in the pY1H list.

For every TF in the list, checks its GTRD ChIP-seq peaks against every
cytokine promoter region (all TFs x all cytokines, not just the pairs
reported as PDIs - the full set is needed as the reference network for
s3's randomization test).

Requirements
------------
    pip install pandas openpyxl pybbi

Inputs
------
    inputs/py1h PDIs for ChIP-seq comparison.xlsx
    inputs/ChIP-seq.metadata.txt                   (from GTRD upload, s0)
    outputs/all_tfs_chipseq.txt                    (from s1)
    GTRD MACS2 peak files (set GTRD_MACS2_DIR_OVERRIDE for a local copy)

Outputs
-------
    outputs/py1h_macs2_peaks_full.xlsx    every peak match, TF x cytokine
    outputs/py1h_chipseq_evidence.xlsx    n_evidence per PDI
"""

import os
import re
import sys

import pandas as pd

try:
    import bbi
except ImportError:
    sys.exit(
        "Missing dependency 'bbi'. Install with: pip install pybbi\n"
        "(bbi is required to read the GTRD MACS2 .bb peak files.)"
    )

HERE = os.path.dirname(os.path.abspath(__file__))  # analysis_2026/scripts/
ANALYSIS_DIR = os.path.dirname(HERE)  # analysis_2026/
ROOT = os.path.dirname(ANALYSIS_DIR)  # TF_Heterodimers_Clean/
LAB_ROOT = os.path.dirname(ROOT)  # BostonLab/

INPUT_XLSX = os.path.join(ANALYSIS_DIR, "inputs", "py1h PDIs for ChIP-seq comparison.xlsx")
TF_CHIP_LIST = os.path.join(ANALYSIS_DIR, "outputs", "all_tfs_chipseq.txt")
METADATA_TXT = os.path.join(ANALYSIS_DIR, "inputs", "ChIP-seq.metadata.txt")
# Override with GTRD_MACS2_DIR_OVERRIDE if you've staged a local copy of
# the needed files on a faster filesystem (e.g. native ext4 instead of an
# NTFS mount via WSL) - repeated random-access bbi.fetch_intervals() calls
# over ~3,500 files were dramatically slower over an NTFS mount in testing.
GTRD_MACS2_DIR = os.environ.get("GTRD_MACS2_DIR_OVERRIDE") or os.path.join(
    LAB_ROOT,
    "Anna_heterodimer_sequencing",
    "gtrd",
    "gtrd2",
    "egrid",
    "bigBeds",
    "hg38",
    "ChIP-seq",
    "Peaks",
    "MACS2",
)
OUTPUT_DIR = os.path.join(ANALYSIS_DIR, "outputs")
FULL_PEAKS_XLSX = os.path.join(OUTPUT_DIR, "py1h_macs2_peaks_full.xlsx")
PDI_EVIDENCE_XLSX = os.path.join(OUTPUT_DIR, "py1h_chipseq_evidence.xlsx")

SEARCH_MARGIN = 200_000  # window around each cytokine bait to search for peaks

COORD_RE = re.compile(r"^(chr[\w]+):(\d+)-(\d+)$", re.IGNORECASE)


def parse_coords(coord_str):
    m = COORD_RE.match(str(coord_str).strip())
    if not m:
        raise ValueError(f"Unrecognized coordinate format: {coord_str!r}")
    chrom, start, end = m.group(1).lower(), int(m.group(2)), int(m.group(3))
    if end < start:
        start, end = end, start
    return chrom, start, end


def load_tf_chip_availability():
    if not os.path.exists(TF_CHIP_LIST):
        sys.exit(
            f"Missing {TF_CHIP_LIST}.\n"
            "Run s1_generate_all_tfs_chipseq.py first to produce it."
        )
    htfs = pd.read_csv(TF_CHIP_LIST, sep="\t")
    # add TFAP2A / AP2A, Uniprot P05549 (missing from the GTRD-derived list)
    htfs = pd.concat(
        [htfs, pd.DataFrame([{"ID": "AP2A", "Uniprot": "P05549"}])],
        ignore_index=True,
    )
    return set(htfs["ID"].unique())


def index_gtrd_files():
    if not os.path.isdir(GTRD_MACS2_DIR):
        sys.exit(f"GTRD MACS2 peak directory not found:\n  {GTRD_MACS2_DIR}")
    files_by_tf = {}
    for fname in os.listdir(GTRD_MACS2_DIR):
        # filename pattern: PEAKS<id>_<TF>_<uniprot>_MACS2_<n>.bb
        parts = fname.split("_")
        if len(parts) < 2:
            continue
        tf = parts[1]
        files_by_tf.setdefault(tf, []).append(fname)
    return files_by_tf


def build_full_cross_product_peaks(pdis, tfs_in_chip, files_by_tf, metadata):
    """Every TF in the input list x every cytokine bait in the input
    list, searched against that TF's GTRD peak files."""

    cytokines = pdis["cytokine"].unique().tolist()
    all_tfs = [t for t in pdis["tf"].unique().tolist() if t in tfs_in_chip]

    coords_by_cytokine = {}
    for cytokine in cytokines:
        coord_values = pdis.loc[pdis["cytokine"] == cytokine,
                                 "Genomic Coordinates (GRCh38/hg38)"].unique()
        if len(coord_values) != 1:
            print(f"WARNING: {cytokine} has {len(coord_values)} distinct "
                  f"coordinate values, using the first: {coord_values}")
        coords_by_cytokine[cytokine] = parse_coords(coord_values[0])

    peak_rows = []

    # TF is the outer loop so each TF's GTRD files are opened once,
    # queried against all cytokines, then closed before the next TF.
    for i, tf in enumerate(all_tfs, 1):
        tf_handles = []
        for fname in files_by_tf.get(tf, []):
            try:
                bbi_file = bbi.open(os.path.join(GTRD_MACS2_DIR, fname))
            except OSError as e:
                print(f"  WARNING: skipping unreadable file {fname}: {e}")
                continue
            peak_id = fname.split("_")[0]
            tf_handles.append((bbi_file, peak_id))

        for cytokine in cytokines:
            chrom, start, end = coords_by_cytokine[cytokine]

            for bbi_file, peak_id in tf_handles:
                chromsizes = dict(bbi_file.chromsizes)
                if chrom not in chromsizes:
                    continue

                peaks = bbi_file.fetch_intervals(
                    chrom, max(0, start - SEARCH_MARGIN), end + SEARCH_MARGIN
                )
                if peaks.shape[0] == 0:
                    continue

                peaks = peaks[(peaks["abs_summit"] >= start) & (peaks["abs_summit"] <= end)]
                if peaks.shape[0] == 0:
                    continue

                tmp_info = metadata[metadata["input"] == peak_id]

                for _, peak in peaks.iterrows():
                    peak_rows.append({
                        "tf": tf,
                        "cytokine": cytokine,
                        "name": f"{tf}-{cytokine}",
                        "peak_id": peak_id,
                        "abs_summit": peak["abs_summit"],
                        "cell_line": tmp_info.iloc[0, 8] if tmp_info.shape[0] else None,
                        "antibody": tmp_info.iloc[0, 1] if tmp_info.shape[0] else None,
                        "treatment": tmp_info.iloc[0, 3] if tmp_info.shape[0] else None,
                    })

        for bbi_file, _ in tf_handles:
            bbi_file.close()

        print(f"  processed TF {tf} ({i}/{len(all_tfs)}, {len(tf_handles)} files) "
              f"across {len(cytokines)} cytokines")

    return pd.DataFrame(peak_rows)


def write_pdi_evidence_sheet(pdis, tfs_in_chip, full_peaks):
    filtered = pdis[pdis["tf"].isin(tfs_in_chip)].copy()
    if full_peaks.empty:
        filtered["n_evidence"] = 0
    else:
        counts = full_peaks.groupby(["tf", "cytokine"]).size().rename("n_evidence")
        filtered = filtered.merge(counts, how="left", left_on=["tf", "cytokine"],
                                   right_index=True)
        filtered["n_evidence"] = filtered["n_evidence"].fillna(0).astype(int)

    name_col = "TF-cytokine" if "TF-cytokine" in filtered.columns else None
    out = pd.DataFrame({
        "cytokine": filtered["cytokine"],
        "tf": filtered["tf"],
        "n_evidence": filtered["n_evidence"],
        "name": filtered[name_col] if name_col else filtered["tf"] + "-" + filtered["cytokine"],
    })
    out.to_excel(PDI_EVIDENCE_XLSX, sheet_name="ChIP-seq results", index=False)

    n_with_evidence = (out["n_evidence"] > 0).sum()
    print(f"\nWrote {PDI_EVIDENCE_XLSX}")
    print(f"PDIs with TF ChIP-seq data available: {len(out)}")
    print(f"PDIs with >=1 ChIP-seq peak: {n_with_evidence}")
    print(f"PDIs with TF ChIP-seq data but 0 peaks (no-evidence set): "
          f"{len(out) - n_with_evidence}")
    return out


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    pdis = pd.read_excel(INPUT_XLSX, sheet_name="Sheet1")
    tfs_in_chip = load_tf_chip_availability()
    files_by_tf = index_gtrd_files()
    metadata = pd.read_csv(METADATA_TXT, sep="\t")
    metadata = metadata[metadata["specie"] == "Homo sapiens"]

    print(f"PDIs total: {len(pdis)}")
    print(f"Unique TFs in list with GTRD data: "
          f"{len([t for t in pdis['tf'].unique() if t in tfs_in_chip])} / {pdis['tf'].nunique()}")
    print(f"Unique cytokines in list: {pdis['cytokine'].nunique()}")

    full_peaks = build_full_cross_product_peaks(pdis, tfs_in_chip, files_by_tf, metadata)
    full_peaks.to_excel(FULL_PEAKS_XLSX, sheet_name="peaks", index=False)
    print(f"\nWrote {FULL_PEAKS_XLSX} ({len(full_peaks)} peak rows, "
          f"{full_peaks[['tf','cytokine']].drop_duplicates().shape[0] if len(full_peaks) else 0} unique TF-cytokine edges)")

    write_pdi_evidence_sheet(pdis, tfs_in_chip, full_peaks)


if __name__ == "__main__":
    main()
