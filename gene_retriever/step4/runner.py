#!/usr/bin/env python3
import pandas as pd
import os
import sys
from collections import defaultdict
from genes.normalisation import clean_gene_term, generate_gene_patterns, normalise_hyphens
from filesystem.paths import save_fasta_step4, sanitize_filename
from ncbi.entrez_client import query_ncbi, efetch_ncbi, get_assembly_accession_from_nuccore

def parse_gene_list_from_file(genes_not_found_file):
    """Read genes from GENES_NOT_FOUND.txt"""
    genes = []
    if os.path.exists(genes_not_found_file):
        with open(genes_not_found_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.append(gene)
    return genes

def safe_int(x, default=0):
    try:
        return int(x)
    except (ValueError, TypeError):
        return default

def determine_taxa(organism_str):
    """Determine taxa from organism string"""
    organism_lower = organism_str.lower()
    if any(term in organism_lower for term in ['virus', 'phage', 'viridae']):
        return "Viruses"
    elif any(term in organism_lower for term in ['archaeon', 'archaea', 'archaeal']):
        return "Archaea"
    elif any(term in organism_lower for term in ['fungus', 'fungi', 'yeast', 'mycota']):
        return "Fungi"
    else:
        return "Bacteria"

def fetch_gene_any_organism(gene_symbol_raw, taxa_filter, max_results=1000, batch_size=50):
    """
    Fetch gene sequences from NCBI for any organism with taxa filter.
    Returns: cds_recs (list), seq_meta (list of dicts)
    """
    patterns = generate_gene_patterns(gene_symbol_raw)
    uid_list = []
    
    for pat in patterns:
        query = f"{pat}[Gene] AND {taxa_filter}[filter]"
        ids = query_ncbi(query, "nucleotide", max_results)
        if ids:
            print(f"Found {len(ids)} hits for pattern '{pat}' with filter '{taxa_filter}', stopping further pattern checks.", flush=True)
            uid_list = ids
            break
        else:
            print(f"No hits for pattern '{pat}' with filter '{taxa_filter}'.", flush=True)
    
    if not uid_list:
        print(f"No search results found for {gene_symbol_raw} with filter '{taxa_filter}'.", flush=True)
        return [], []
    
    print(f"Using {len(uid_list)} UIDs for gene '{gene_symbol_raw}' with filter '{taxa_filter}'.", flush=True)
    
    cds_recs = []
    seq_meta = []
    
    # Process in batches
    total_uids = len(uid_list)
    for batch_start in range(0, total_uids, batch_size):
        batch_end = min(batch_start + batch_size, total_uids)
        batch_ids = uid_list[batch_start:batch_end]
        
        print(f"\n--- Processing batch {batch_start//batch_size + 1}: UIDs {batch_start+1}-{batch_end} of {total_uids} ---", flush=True)
        
        for i, idx in enumerate(batch_ids, 1):
            print(f"Fetching result {batch_start + i}/{total_uids} (ID: {idx})", flush=True)
            
            gb_record = efetch_ncbi(idx, "gb", "nucleotide", 1)
            if not gb_record:
                print(f"No GenBank record for UID {idx}", flush=True)
                continue
            
            fasta_record = efetch_ncbi(idx, "fasta", "nucleotide", 1)
            if not fasta_record:
                print("No FASTA sequence found for this entry.", flush=True)
                continue
            
            organism = gb_record.annotations.get("organism", gb_record.description)
            
            # Extract strain/isolate info
            strain = None
            for feat in gb_record.features:
                if feat.type == "source":
                    q = feat.qualifiers
                    strain = q.get("strain", [None])[0] or q.get("isolate", [None])[0] or q.get("serovar", [None])[0]
                    break
            
            organism_full = f"{organism} {strain}" if strain else organism
            
            # Parse genus and species
            parts = organism.split()
            genus = parts[0] if len(parts) > 0 else "Unknown"
            species = parts[1] if len(parts) > 1 else "sp"
            
            nuccore_id = gb_record.id
            assembly_acc = get_assembly_accession_from_nuccore(nuccore_id)
            
            taxa = determine_taxa(organism)
            
            per_uid_seq_count = 0
            
            for rec in gb_record.features:
                if rec.type != "CDS":
                    continue
                
                gene_names = [normalise_hyphens(g).lower() for g in rec.qualifiers.get("gene", [])]
                if not gene_names:
                    continue
                
                if not any(g in patterns for g in gene_names):
                    continue
                
                try:
                    cds_seq = rec.extract(fasta_record.seq)
                    print(f"=== FOUND {gene_symbol_raw} at {rec.location} of {gb_record.description} -- (Nuccore: {nuccore_id}, Assembly: {assembly_acc}, Organism: {organism_full}) ===", flush=True)
                    per_uid_seq_count += 1
                    
                    header = (
                        f">{gb_record.id} [gene={rec.qualifiers.get('gene', [gene_symbol_raw])[0]}] "
                        f"[locus_tag={rec.qualifiers.get('locus_tag', [''])[0]}] "
                        f"[protein={rec.qualifiers.get('product', [''])[0]}] "
                        f"[protein_id={rec.qualifiers.get('protein_id', [''])[0]}] "
                        f"[location={rec.location}] [organism={organism_full}] "
                        f"[description={gb_record.description}]"
                    )
                    cds_rec = f"{header}\n{cds_seq}"
                    cds_recs.append(cds_rec)
                    
                    file_gene_name = normalise_hyphens(f"{rec.qualifiers.get('gene', [gene_symbol_raw])[0]}").strip()
                    
                    seq_meta.append({
                        "organism": organism_full,
                        "genus": genus,
                        "species": species,
                        "assembly_accession": assembly_acc,
                        "nuccore_accession": nuccore_id,
                        "fna_gene_name": file_gene_name,
                        "taxa": taxa
                    })
                except Exception as e:
                    print(f"Error extracting CDS: {e}", flush=True)
            
            if per_uid_seq_count == 0:
                print(f"!!! CDS sequence for {gene_symbol_raw} not found in UID {idx}.", flush=True)
    
    return cds_recs, seq_meta

def run_step4(in_gene_stats, genes_not_found_file, gene_sequence_dir):
    print("="*80)
    print("STEP 4: Querying NCBI with gene only (any organism)")
    print("="*80)
    
    out_gene_stats = "4_gene_stats_found.tsv"
    
    # Load gene stats
    gene_stats_in = pd.read_csv(in_gene_stats, sep="\t")
    
    # Normalize existing fields
    sep_cols = [
        "STEP-1_ORGS", "STEP-1_PATHS", "STEP-2_ORGS", "STEP-2_PATHS",
        "STEP-3_GENUS", "STEP-3_SOURCE_ORGS", "STEP-3_ACCESSIONS", "STEP-3_PATHS",
        "STEP-4_TAXA", "STEP-4_SOURCE_ORGS", "STEP-4_ACCESSIONS", "STEP-4_PATHS"
    ]
    
    for col in sep_cols:
        if col in gene_stats_in.columns:
            gene_stats_in[col] = (
                gene_stats_in[col].astype(str)
                .apply(lambda x: "" if x.lower() == "nan" else x.replace(",", "; "))
            )
    
    # Ensure STEP-4 columns exist
    if "STEP-4_NUM_MATCHES" not in gene_stats_in.columns:
        gene_stats_in["STEP-4_NUM_MATCHES"] = 0
    else:
        gene_stats_in["STEP-4_NUM_MATCHES"] = gene_stats_in["STEP-4_NUM_MATCHES"].fillna(0)
    
    for col in ["STEP-4_TAXA", "STEP-4_SOURCE_ORGS", "STEP-4_ACCESSIONS", "STEP-4_PATHS"]:
        if col not in gene_stats_in.columns:
            gene_stats_in[col] = ""
    
    # Build canonical gene index
    gene_index_by_canonical = {}
    if "GENE_NAME" in gene_stats_in.columns:
        for idx, gname in gene_stats_in["GENE_NAME"].items():
            g_can = clean_gene_term(gname)
            if g_can not in gene_index_by_canonical:
                gene_index_by_canonical[g_can] = idx
    
    # Load genes not found
    genes_to_search = parse_gene_list_from_file(genes_not_found_file)
    print(f"Loaded {len(genes_to_search)} genes from {genes_not_found_file}", flush=True)
    
    # Gene-centric stats for STEP-4
    gene_step4_info = defaultdict(lambda: {
        "num_matches": 0,
        "taxa_set": set(),
        "source_orgs": [],
        "assembly_accessions": [],
        "nuccore_accessions": [],
        "paths": set(),
    })
    
    taxa_filters = ["bacteria", "archaea", "fungi", "viruses"]
    
    for gene_raw in genes_to_search:
        gene_raw = gene_raw.strip()
        if not gene_raw:
            continue
        
        gene_can = clean_gene_term(gene_raw)
        gene_safe = sanitize_filename(gene_can)
        
        print(f"\n\n{'='*80}\nProcessing gene: {gene_raw} (canonical: {gene_can})\n{'='*80}\n", flush=True)
        
        all_cds_records = []
        all_seq_meta = []
        
        # Query all 4 taxa filters
        for taxa_filter in taxa_filters:
            print(f"\n--- Querying {taxa_filter} for gene {gene_raw} ---", flush=True)
            cds_records, seq_meta = fetch_gene_any_organism(gene_raw, taxa_filter, max_results=1000, batch_size=50)
            
            if cds_records:
                print(f"✅ Found {len(cds_records)} sequences for {gene_raw} in {taxa_filter}", flush=True)
                all_cds_records.extend(cds_records)
                all_seq_meta.extend(seq_meta)
            else:
                print(f"✗ No sequences found for {gene_raw} in {taxa_filter}", flush=True)
        
        if not all_cds_records:
            print(f"\n✗ NO SEQUENCES FOUND for {gene_raw} across all taxa", flush=True)
            continue
        
        # Group sequences by organism and save
        organism_groups = defaultdict(list)
        organism_meta = {}
        
        for cds_rec, meta in zip(all_cds_records, all_seq_meta):
            key = (meta["assembly_accession"], meta["genus"], meta["species"])
            organism_groups[key].append(cds_rec)
            if key not in organism_meta:
                organism_meta[key] = meta
        
        print(f"\n✅ Total sequences found for {gene_raw}: {len(all_cds_records)} from {len(organism_groups)} organisms", flush=True)
        
        saved_paths = []
        
        for (assembly_acc, genus, species), records in organism_groups.items():
            meta = organism_meta[(assembly_acc, genus, species)]
            fasta_path = save_fasta_step4(records, gene_safe, assembly_acc, genus, species, gene_sequence_dir)
            
            if fasta_path:
                print(f"Saved {len(records)} sequences to {fasta_path}", flush=True)
                saved_paths.append(fasta_path)
                
                # Update gene_step4_info
                ginfo = gene_step4_info[gene_can]
                ginfo["num_matches"] += len(records)
                ginfo["taxa_set"].add(meta["taxa"])
                ginfo["source_orgs"].append(meta["organism"])
                ginfo["assembly_accessions"].append(meta["assembly_accession"])
                ginfo["nuccore_accessions"].append(meta["nuccore_accession"])
                ginfo["paths"].add(fasta_path)
        
        print(f"\n✅ Saved gene {gene_raw} to {len(saved_paths)} organism directories", flush=True)
    
    # Update gene stats
    for gene_can, info in gene_step4_info.items():
        num_new = info["num_matches"]
        if num_new <= 0:
            continue
        
        taxa_str = ";".join(sorted(info["taxa_set"])) if info["taxa_set"] else ""
        source_orgs_str = ";".join(info["source_orgs"]) if info["source_orgs"] else ""
        accessions_str = ";".join(info["assembly_accessions"]) if info["assembly_accessions"] else ""
        paths_str = ";".join(sorted(info["paths"])) if info["paths"] else ""
        
        idx = gene_index_by_canonical.get(gene_can)
        
        if idx is not None:
            # Update existing row
            prev_total = safe_int(gene_stats_in.loc[idx, "TOTAL_MATCHES"])
            gene_stats_in.loc[idx, "TOTAL_MATCHES"] = prev_total + num_new
            
            prev_step4 = safe_int(gene_stats_in.loc[idx, "STEP-4_NUM_MATCHES"])
            gene_stats_in.loc[idx, "STEP-4_NUM_MATCHES"] = prev_step4 + num_new
            
            def append_field(col_name, add_str):
                if not add_str:
                    return
                prev = str(gene_stats_in.loc[idx, col_name]) if col_name in gene_stats_in.columns else ""
                prev = prev if prev != "nan" else ""
                prev = prev.strip()
                if prev:
                    gene_stats_in.loc[idx, col_name] = prev + "; " + add_str
                else:
                    gene_stats_in.loc[idx, col_name] = add_str
            
            append_field("STEP-4_TAXA", taxa_str)
            append_field("STEP-4_SOURCE_ORGS", source_orgs_str)
            append_field("STEP-4_ACCESSIONS", accessions_str)
            append_field("STEP-4_PATHS", paths_str)
        else:
            # New gene discovered only in step 4
            new_row = {
                "GENE_NAME": gene_can,
                "TOTAL_MATCHES": num_new,
                "STEP-1_NUM_MATCHES": 0,
                "STEP-1_ORGS": "",
                "STEP-1_PATHS": "",
                "STEP-2_NUM_MATCHES": 0,
                "STEP-2_ORGS": "",
                "STEP-2_PATHS": "",
                "STEP-3_NUM_MATCHES": 0,
                "STEP-3_GENUS": "",
                "STEP-3_SOURCE_ORGS": "",
                "STEP-3_ACCESSIONS": "",
                "STEP-3_PATHS": "",
                "STEP-4_NUM_MATCHES": num_new,
                "STEP-4_TAXA": taxa_str,
                "STEP-4_SOURCE_ORGS": source_orgs_str,
                "STEP-4_ACCESSIONS": accessions_str,
                "STEP-4_PATHS": paths_str,
            }
            
            for col in gene_stats_in.columns:
                if col not in new_row:
                    new_row[col] = "" if gene_stats_in[col].dtype == object else 0
            
            gene_stats_in = pd.concat([gene_stats_in, pd.DataFrame([new_row])], ignore_index=True)
            gene_index_by_canonical[gene_can] = gene_stats_in.index[-1]
    
    # Sort and save
    if "GENE_NAME" in gene_stats_in.columns:
        gene_stats_in.sort_values(by="GENE_NAME", inplace=True)
    
    gene_stats_in.to_csv(out_gene_stats, sep="\t", index=False)
    print(f"\nWrote {out_gene_stats}", flush=True)
    
    print("\nStep 4 completed ✅", flush=True)
    return True
