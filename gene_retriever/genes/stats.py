#!/usr/bin/env python3

import pandas as pd

def update_gene_stats(found_genes, genus, species, accession, fasta_paths, old_gene_stats):
    """
    Update gene statistics DataFrame.

    Columns (final 2_gene_stats_found.tsv):
    col1: GENE_NAME
    col2: TOTAL_MATCHES          (step1 + step2)
    col3: STEP-1_NUM_MATCHES     (from input, unchanged)
    col4: STEP-1_ORGS            (from input, unchanged)
    col5: STEP-1_PATHS           (from input, unchanged)
    col6: STEP-2_NUM_MATCHES     (matches found in THIS run only)
    col7: STEP-2_ORGS            (orgs for THIS run only, comma-separated)
    col8: STEP-2_PATHS           (paths for THIS run only, comma-separated)
    """
    if not isinstance(found_genes, (list, tuple)):
        found_genes = []
    if not isinstance(fasta_paths, (list, tuple)):
        fasta_paths = []

    new_rows = []

    for gene, fasta_path in zip(found_genes, fasta_paths):
        if gene in old_gene_stats["GENE_NAME"].values:
            # Existing gene: bump total and step-2 counts
            idx = old_gene_stats.index[old_gene_stats["GENE_NAME"] == gene][0]

            # TOTAL_MATCHES = previous TOTAL_MATCHES + 1 (robust cast)
            prev_total = old_gene_stats.loc[idx, "TOTAL_MATCHES"]
            try:
                prev_total_int = int(prev_total)
            except (ValueError, TypeError):
                prev_total_int = 0
            old_gene_stats.loc[idx, "TOTAL_MATCHES"] = prev_total_int + 1

            # STEP-2_NUM_MATCHES = previous STEP-2_NUM_MATCHES + 1
            prev_step2 = old_gene_stats.loc[idx, "STEP-2_NUM_MATCHES"]
            try:
                prev_step2_int = int(prev_step2)
            except (ValueError, TypeError):
                prev_step2_int = 0
            old_gene_stats.loc[idx, "STEP-2_NUM_MATCHES"] = prev_step2_int + 1

            # append orgs and paths for this run
            prev_orgs = str(old_gene_stats.loc[idx, "STEP-2_ORGS"])
            prev_paths = str(old_gene_stats.loc[idx, "STEP-2_PATHS"])

            if prev_orgs.strip():
                new_orgs = prev_orgs + f",{genus} {species}"
            else:
                new_orgs = f"{genus} {species}"

            if prev_paths.strip():
                new_paths = prev_paths + f",{fasta_path}"
            else:
                new_paths = fasta_path

            old_gene_stats.loc[idx, "STEP-2_ORGS"] = new_orgs
            old_gene_stats.loc[idx, "STEP-2_PATHS"] = new_paths

        else:
            # New gene: only discovered in this run
            new_rows.append({
                "GENE_NAME": gene,
                "TOTAL_MATCHES": 1,             # only step2 matches so far
                "STEP-1_NUM_MATCHES": 0,
                "STEP-1_ORGS": "",
                "STEP-1_PATHS": "",
                "STEP-2_NUM_MATCHES": 1,
                "STEP-2_ORGS": f"{genus} {species}",
                "STEP-2_PATHS": fasta_path,
            })

    if new_rows:
        old_gene_stats = pd.concat([old_gene_stats, pd.DataFrame(new_rows)], ignore_index=True)

    return old_gene_stats
