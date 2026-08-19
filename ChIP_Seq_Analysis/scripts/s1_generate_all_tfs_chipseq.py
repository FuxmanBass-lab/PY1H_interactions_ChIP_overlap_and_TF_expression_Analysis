"""
Generates all_tfs_chipseq.txt: the list of human TFs that have GTRD
ChIP-seq data, with their Uniprot IDs.

Reads ChIP-seq.metadata.txt, gets the unique TF Uniprot IDs, and
queries Ensembl BioMart to get the gene symbol for each one. IDs
BioMart can't resolve are filled in from manual_lookups/no_found_biomart.txt.
Anything still unresolved goes to outputs/unresolved_uniprot_ids.txt.

Requirements
------------
    pip install pandas biomart

Inputs
------
    inputs/ChIP-seq.metadata.txt
    manual_lookups/no_found_biomart.txt

Outputs
-------
    outputs/all_tfs_chipseq.txt          ID<TAB>Uniprot
    outputs/unresolved_uniprot_ids.txt   only written if non-empty
"""

import math
import os
import sys

import pandas as pd

try:
    import biomart
except ImportError:
    sys.exit("Missing dependency 'biomart'. Install with: pip install biomart")

HERE = os.path.dirname(os.path.abspath(__file__))  # analysis_2026/scripts/
ANALYSIS_DIR = os.path.dirname(HERE)  # analysis_2026/

METADATA_TXT = os.path.join(ANALYSIS_DIR, "inputs", "ChIP-seq.metadata.txt")
MANUAL_LOOKUP_TXT = os.path.join(ANALYSIS_DIR, "manual_lookups", "no_found_biomart.txt")
OUTPUT_DIR = os.path.join(ANALYSIS_DIR, "outputs")
OUTPUT_TXT = os.path.join(OUTPUT_DIR, "all_tfs_chipseq.txt")
UNRESOLVED_TXT = os.path.join(OUTPUT_DIR, "unresolved_uniprot_ids.txt")

BATCH_SIZE = 20


def load_metadata_uniprot_ids():
    df = pd.read_csv(METADATA_TXT, sep="\t")
    df = df[df["specie"] == "Homo sapiens"]
    ids = df["tf_uniprot_id"].dropna().unique().tolist()
    ids = [i for i in ids if i != "NULL"]
    return sorted(set(ids))


def query_biomart(uniprot_ids):
    server = biomart.BiomartServer("http://useast.ensembl.org/biomart")
    dataset = server.datasets["hsapiens_gene_ensembl"]

    resolved = {}  # uniprot -> gene symbol
    n_batches = math.ceil(len(uniprot_ids) / BATCH_SIZE)
    for i in range(n_batches):
        batch = uniprot_ids[i * BATCH_SIZE:(i + 1) * BATCH_SIZE]
        response = dataset.search({
            "filters": {"uniprot_gn_id": batch},
            "attributes": ["hgnc_symbol", "uniprot_gn_id"],
        })
        for line in response.iter_lines():
            symbol, uniprot = line.decode("utf-8").split("\t")
            if symbol:
                resolved[uniprot] = symbol
        print(f"BioMart batch {i + 1}/{n_batches} done")
    return resolved


def load_manual_lookup():
    """no_found_biomart.txt is Uniprot<TAB>ID (reversed vs. BioMart's
    output order) - flipped here so everything downstream is ID, Uniprot."""
    manual = {}
    if not os.path.exists(MANUAL_LOOKUP_TXT):
        print(f"WARNING: {MANUAL_LOOKUP_TXT} not found, skipping manual lookups")
        return manual
    with open(MANUAL_LOOKUP_TXT) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            uniprot, symbol = line.split("\t")
            manual[uniprot] = symbol
    return manual


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    uniprot_ids = load_metadata_uniprot_ids()
    print(f"Unique human TF Uniprot IDs in ChIP-seq.metadata.txt: {len(uniprot_ids)}")

    resolved = query_biomart(uniprot_ids)
    print(f"Resolved by BioMart: {len(resolved)} / {len(uniprot_ids)}")

    manual = load_manual_lookup()

    still_missing = [u for u in uniprot_ids if u not in resolved and u not in manual]

    rows = [{"ID": symbol, "Uniprot": uniprot} for uniprot, symbol in resolved.items()]
    rows += [{"ID": symbol, "Uniprot": uniprot} for uniprot, symbol in manual.items()
             if uniprot in uniprot_ids]

    out_df = pd.DataFrame(rows).drop_duplicates()
    out_df.to_csv(OUTPUT_TXT, sep="\t", index=False)
    print(f"\nWrote {OUTPUT_TXT} ({len(out_df)} TFs)")

    if still_missing:
        with open(UNRESOLVED_TXT, "w") as f:
            f.write("\n".join(still_missing) + "\n")
        print(f"{len(still_missing)} Uniprot IDs unresolved - see {UNRESOLVED_TXT}")
    else:
        print("All TFs resolved.")


if __name__ == "__main__":
    main()
