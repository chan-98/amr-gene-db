#!/usr/bin/env python3
import os
import sys
import traceback
import pandas as pd
from genes.stats import update_gene_stats
from filesystem.paths import find_existing_fasta_for_gene, save_fasta_step2
from ncbi.entrez_client import fetch_gene_from_species

def process_input_row(row, output_basedir):
    genus = str(row.get("GENUS", "") or "")
    species = str(row.get("SPECIES", "") or "")
    accession = str(row.get("ACCESSION_STATUS", "") or "")
    
    if not accession or accession.upper() == "ORGANISM_NOT_FOUND":
        print(f"Invalid accession: {accession}", flush=True)
        return row.get("INPUT_GENES", ""), row.get("GENES_FOUND", ""), row.get("GENES_NOT_FOUND", ""), None, []
    
    genes_not_found_raw = row.get("GENES_NOT_FOUND", "")
    if isinstance(genes_not_found_raw, float):
        genes_not_found_raw = ""
    genes_to_search = [g.strip() for g in str(genes_not_found_raw).split(",") if g.strip()]
    print(f"Searching for the genes: {', '.join(genes_to_search)}", flush=True)
    
    found, not_found, fasta_paths = [], [], []
    
    for gene_symbol in genes_to_search:
        if not gene_symbol:
            continue
        
        print(f"\nGene undergoing pre-check: {gene_symbol}", flush=True)
        existing_paths = find_existing_fasta_for_gene(gene_symbol, accession, genus, species, output_basedir)
        
        if existing_paths:
            print(f"    ⤴ Skipping NCBI query for {gene_symbol} — existing file(s) found: {existing_paths}", flush=True)
            for pth in existing_paths:
                fasta_paths.append(pth)
            found.append(gene_symbol)
            continue
        
        cds_records = fetch_gene_from_species(gene_symbol, genus, species, 10)
        if cds_records:
            fasta_path = save_fasta_step2(cds_records, gene_symbol, genus, species, accession, output_basedir)
            print(f"\n✅ === Saved {gene_symbol} for {genus} {species} →→→ {fasta_path} ===", flush=True)
            found.append(gene_symbol)
            fasta_paths.append(fasta_path)
        else:
            print(f"    ✗ NOT FOUND: {gene_symbol}", flush=True)
            not_found.append(gene_symbol)
    
    already_found_raw = row.get("GENES_FOUND", "")
    if isinstance(already_found_raw, float):
        already_found_raw = ""
    already_found = [g.strip() for g in str(already_found_raw).split(",") if g.strip() and g.strip().upper() != "NIL"]
    all_found = list(dict.fromkeys(already_found + found))
    
    print(f"\nSUMMARY:", flush=True)
    print(f"  Genes to search: {len(genes_to_search)}", flush=True)
    print(f"  Found: {len(found)}", flush=True)
    print(f"  Not found: {len(not_found)}", flush=True)
    
    return row.get("INPUT_GENES", ""), ",".join(all_found), ",".join(not_found), fasta_paths, found

def run_step2(in_run_stats, in_gene_stats, gene_sequence_dir):
    print("="*80)
    print("STEP 2: Querying NCBI with gene+genus+species")
    print("="*80)
    
    out_gene_stats = "2_gene_stats_found.tsv"
    out_run_stats = "2_run_stats.tsv"
    
    in_df = pd.read_csv(in_run_stats, sep="\t")
    gene_stats_in = pd.read_csv(in_gene_stats, sep="\t")
    
    if "STEP-2_NUM_MATCHES" not in gene_stats_in.columns:
        gene_stats_in["STEP-2_NUM_MATCHES"] = 0
    else:
        gene_stats_in["STEP-2_NUM_MATCHES"] = 0
    if "STEP-2_ORGS" not in gene_stats_in.columns:
        gene_stats_in["STEP-2_ORGS"] = ""
    else:
        gene_stats_in["STEP-2_ORGS"] = ""
    if "STEP-2_PATHS" not in gene_stats_in.columns:
        gene_stats_in["STEP-2_PATHS"] = ""
    else:
        gene_stats_in["STEP-2_PATHS"] = ""
    
    run_rows = []
    
    for _, row in in_df.iterrows():
        genus, species, acc = row.get("GENUS", ""), row.get("SPECIES", ""), row.get("ACCESSION_STATUS", "")
        print(f"\n\n======================================== Processing {genus} {species} ({acc}) ========================================\n", flush=True)
        
        try:
            input_genes, found_str, not_found_str, cds_paths, new_found = process_input_row(row, gene_sequence_dir)
            run_rows.append([genus, species, acc, input_genes, found_str, not_found_str])
            if cds_paths and new_found:
                gene_stats_in = update_gene_stats(new_found, genus, species, acc, cds_paths, gene_stats_in)

        except Exception as e:
            print(f"===Error processing {genus} {species}: {e}", flush=True)
            traceback.print_exc()
            run_rows.append([
                genus, species, acc,
                row.get("INPUT_GENES", ""),
                row.get("GENES_FOUND", ""),
                row.get("GENES_NOT_FOUND", "")
            ])

        finally:
            run_df = pd.DataFrame(run_rows, columns=[
                "GENUS", "SPECIES", "ACCESSION_STATUS", "INPUT_GENES", "GENES_FOUND", "GENES_NOT_FOUND"
            ])
            run_df.to_csv(out_run_stats, sep="\t", index=False)
            gene_stats_in.sort_values(by="GENE_NAME", inplace=True)
            gene_stats_in.to_csv(out_gene_stats, sep="\t", index=False)
            print("UPDATED 2_run_stats.tsv AND 2_gene_stats_found.tsv", flush=True)

    print("\nStep 2 completed ✅", flush=True)
    return True
