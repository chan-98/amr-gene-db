import pandas as pd
import os
import sys
from collections import defaultdict
from genes.normalisation import clean_gene_term
from filesystem.paths import save_fasta_step3
from ncbi.entrez_client import fetch_gene_for_genus

def parse_gene_list(val):
    if val is None:
        return []
    try:
        if pd.isna(val):
            return []
    except Exception:
        pass
    s = str(val).strip()
    if not s:
        return []
    genes = [g.strip() for g in s.split(",")]
    return [g for g in genes if g and g.upper() != "NIL"]

def safe_int(x, default=0):
    try:
        return int(x)
    except (ValueError, TypeError):
        return default

def run_step3(in_run_stats, in_gene_stats, gene_sequence_dir):
    print("="*80)
    print("STEP 3: Querying NCBI with gene+genus only")
    print("="*80)
    
    out_run_stats = "3_run_stats.tsv"
    out_gene_stats = "3_gene_stats_found.tsv"
    out_genus_stats = "3_genus_stats_found.tsv"
    
    in_df = pd.read_csv(in_run_stats, sep="\t")
    gene_stats_in = pd.read_csv(in_gene_stats, sep="\t")
    
    sep_cols = ["STEP-1_ORGS", "STEP-1_PATHS", "STEP-2_ORGS", "STEP-2_PATHS", 
                "STEP-3_GENUS", "STEP-3_SOURCE_ORGS", "STEP-3_ACCESSIONS", "STEP-3_PATHS"]
    
    for col in sep_cols:
        if col in gene_stats_in.columns:
            gene_stats_in[col] = (
                gene_stats_in[col].astype(str)
                .apply(lambda x: "" if x.lower() == "nan" else x.replace(",", "; "))
            )
    
    if "STEP-3_NUM_MATCHES" not in gene_stats_in.columns:
        gene_stats_in["STEP-3_NUM_MATCHES"] = 0
    else:
        gene_stats_in["STEP-3_NUM_MATCHES"] = gene_stats_in["STEP-3_NUM_MATCHES"].fillna(0)
    
    for col in ["STEP-3_GENUS", "STEP-3_SOURCE_ORGS", "STEP-3_ACCESSIONS", "STEP-3_PATHS"]:
        if col not in gene_stats_in.columns:
            gene_stats_in[col] = ""
    
    if "STEP-2_NUM_MATCHES" not in gene_stats_in.columns:
        gene_stats_in["STEP-2_NUM_MATCHES"] = 0
    if "STEP-2_ORGS" not in gene_stats_in.columns:
        gene_stats_in["STEP-2_ORGS"] = ""
    if "STEP-2_PATHS" not in gene_stats_in.columns:
        gene_stats_in["STEP-2_PATHS"] = ""
    
    gene_index_by_canonical = {}
    if "GENE_NAME" in gene_stats_in.columns:
        for idx, gname in gene_stats_in["GENE_NAME"].items():
            g_can = clean_gene_term(gname)
            if g_can not in gene_index_by_canonical:
                gene_index_by_canonical[g_can] = idx
    
    genus_initial_found = defaultdict(set)
    genus_gene_input_info = defaultdict(lambda: {"species": set(), "accessions": set()})
    
    for _, row in in_df.iterrows():
        genus = str(row.get("GENUS", "") or "")
        species = str(row.get("SPECIES", "") or "")
        accession = str(row.get("ACCESSION_STATUS", "") or "")
        
        gf_list = parse_gene_list(row.get("GENES_FOUND", ""))
        for g in gf_list:
            g_can = clean_gene_term(g)
            genus_initial_found[genus].add(g_can)
        
        gnf_list = parse_gene_list(row.get("GENES_NOT_FOUND", ""))
        for g in gnf_list:
            g_can = clean_gene_term(g)
            genus_gene_input_info[(genus, g_can)]["species"].add(f"{genus} {species}")
            genus_gene_input_info[(genus, g_can)]["accessions"].add(accession)
    
    genus_found_genes = {g: set(s) for g, s in genus_initial_found.items()}
    genus_nohit_genes = defaultdict(set)
    query_cache = {}
    
    gene_step3_info = defaultdict(lambda: {
        "num_matches": 0, "genus_set": set(), "source_orgs": [],
        "assembly_accessions": [], "nuccore_accessions": [], "paths": set(),
    })
    
    genus_gene_result_info = defaultdict(lambda: {
        "source_orgs": [], "assembly_accessions": [],
        "nuccore_accessions": [], "fna_path": None,
    })
    
    genera = sorted(set(str(g) for g in in_df["GENUS"].dropna()))
    first_genus_for_runstats = True
    
    for genus in genera:
        print(f"\n\n##############################  PROCESSING GENUS: {genus}  ##############################\n\n", flush=True)
        
        genus_mask = in_df["GENUS"].astype(str) == genus
        genus_df = in_df[genus_mask]
        
        for _, row in genus_df.iterrows():
            species = str(row.get("SPECIES", "") or "")
            acc = str(row.get("ACCESSION_STATUS", "") or "")
            
            print(f"\n======================================== Species {genus} {species} ({acc}) ========================================\n", flush=True)
            
            genes_not_found_row = parse_gene_list(row.get("GENES_NOT_FOUND", ""))
            genes_to_search_count = len(genes_not_found_row)
            
            found_this_species = 0
            not_found_this_species = 0
            
            for gene_raw in genes_not_found_row:
                gene_raw = gene_raw.strip()
                if not gene_raw:
                    continue
                
                gene_can = clean_gene_term(gene_raw)
                
                if gene_can in genus_found_genes.get(genus, set()):
                    print(f"Skipping {gene_raw} (canonical {gene_can}) for {genus}: already known FOUND in this genus.", flush=True)
                    found_this_species += 1
                    continue
                
                if gene_can in genus_nohit_genes.get(genus, set()):
                    print(f"Skipping {gene_raw} (canonical {gene_can}) for {genus}: already known NO HIT in this genus.", flush=True)
                    not_found_this_species += 1
                    continue
                
                if (genus, gene_can) in query_cache:
                    result = query_cache[(genus, gene_can)]
                    if result["hit"]:
                        genus_found_genes.setdefault(genus, set()).add(gene_can)
                        found_this_species += 1
                    else:
                        genus_nohit_genes[genus].add(gene_can)
                        not_found_this_species += 1
                    print(f"Reusing cached result for {gene_raw} (canonical {gene_can}) in genus {genus}: hit={result['hit']}", flush=True)
                    continue
                
                print(f"\n--- Querying NCBI for gene {gene_raw} (canonical {gene_can}) in genus {genus} (gene+genus only) ---", flush=True)
                cds_records, seq_meta = fetch_gene_for_genus(gene_raw, genus, max_results=10)
                
                if seq_meta and isinstance(seq_meta, list) and seq_meta[0].get("fna_gene_name"):
                    file_gene_name = seq_meta[0]["fna_gene_name"]
                else:
                    file_gene_name = gene_can
                
                if cds_records:
                    fasta_path = save_fasta_step3(cds_records, file_gene_name, genus, gene_sequence_dir)
                    print(f"\n✅ Saved {len(cds_records)} sequences for {gene_raw} (canonical {gene_can}) in genus {genus} → {fasta_path}", flush=True)
                    
                    genus_found_genes.setdefault(genus, set()).add(gene_can)
                    found_this_species += 1
                    query_cache[(genus, gene_can)] = {"hit": True, "path": fasta_path, "meta": seq_meta}
                    
                    ginfo = gene_step3_info[gene_can]
                    ginfo["num_matches"] += len(seq_meta)
                    ginfo["genus_set"].add(genus)
                    ginfo["paths"].add(fasta_path)
                    for sm in seq_meta:
                        ginfo["source_orgs"].append(sm["organism"])
                        ginfo["assembly_accessions"].append(sm["assembly_accession"])
                        ginfo["nuccore_accessions"].append(sm["nuccore_accession"])
                    
                    rinfo = genus_gene_result_info[(genus, gene_can)]
                    rinfo["fna_path"] = fasta_path
                    for sm in seq_meta:
                        rinfo["source_orgs"].append(sm["organism"])
                        rinfo["assembly_accessions"].append(sm["assembly_accession"])
                        rinfo["nuccore_accessions"].append(sm["nuccore_accession"])
                else:
                    print(f"    ✗ NO HITS for {gene_raw} (canonical {gene_can}) in genus {genus}", flush=True)
                    genus_nohit_genes[genus].add(gene_can)
                    not_found_this_species += 1
                    query_cache[(genus, gene_can)] = {"hit": False, "path": None, "meta": []}
            
            print("\nSUMMARY:", flush=True)
            print(f"  Genes to search (this species): {genes_to_search_count}", flush=True)
            print(f"  Found (for this genus): {found_this_species}", flush=True)
            print(f"  Not found (still missing in this genus): {not_found_this_species}\n", flush=True)
        
        genus_run_rows = []
        for _, row in genus_df.iterrows():
            species = str(row.get("SPECIES", "") or "")
            acc = str(row.get("ACCESSION_STATUS", "") or "")
            input_genes_str = str(row.get("INPUT_GENES", "") or "")
            
            found_set_for_genus = genus_found_genes.get(genus, set())
            found_out = ",".join(sorted(found_set_for_genus)) if found_set_for_genus else ""
            
            orig_not_found_list = parse_gene_list(row.get("GENES_NOT_FOUND", ""))
            not_found_out_list = [g for g in orig_not_found_list if clean_gene_term(g) not in found_set_for_genus]
            not_found_out = ",".join(not_found_out_list) if not_found_out_list else ""
            
            genus_run_rows.append([genus, species, acc, input_genes_str, found_out, not_found_out])
        
        genus_run_df = pd.DataFrame(genus_run_rows, columns=[
            "GENUS", "SPECIES", "ACCESSION_STATUS", "INPUT_GENES", "GENES_FOUND", "GENES_NOT_FOUND"
        ])
        
        mode = "w" if first_genus_for_runstats else "a"
        header = first_genus_for_runstats
        genus_run_df.to_csv(out_run_stats, sep="\t", index=False, mode=mode, header=header)
        first_genus_for_runstats = False
        print(f"Appended run stats for genus {genus} to {out_run_stats}", flush=True)
        
        genus_rows = []
        for (g, gene_can), rinfo in genus_gene_result_info.items():
            if g != genus:
                continue
            fna_path = rinfo.get("fna_path")
            if not fna_path:
                continue
            
            input_info = genus_gene_input_info.get((genus, gene_can), {"species": set(), "accessions": set()})
            input_orgs_str = ";".join(sorted(input_info["species"])) if input_info["species"] else ""
            input_acc_str = ";".join(sorted(input_info["accessions"])) if input_info["accessions"] else ""
            
            source_orgs = rinfo.get("source_orgs", [])
            asm_accs = rinfo.get("assembly_accessions", [])
            nuc_accs = rinfo.get("nuccore_accessions", [])
            
            source_orgs_str = ";".join(source_orgs) if source_orgs else ""
            source_asm_str = ";".join(asm_accs) if asm_accs else ""
            source_nuc_str = ";".join(nuc_accs) if nuc_accs else ""
            
            genus_rows.append([gene_can, input_orgs_str, input_acc_str, source_orgs_str, source_asm_str, source_nuc_str, fna_path])
        
        if genus_rows:
            genus_dir = os.path.join(gene_sequence_dir, genus)
            os.makedirs(genus_dir, exist_ok=True)
            per_genus_stats_path = os.path.join(genus_dir, "stats.tsv")
            per_genus_df = pd.DataFrame(genus_rows, columns=[
                "GENE", "INPUT_ORGS", "INPUT_ACCESSIONS", "SOURCE_ORGS",
                "SOURCE_ASSEMBLY_ACCESSIONS", "SOURCE_NUCCORE_ACCESSIONS", "FNA_PATH"
            ])
            per_genus_df.to_csv(per_genus_stats_path, sep="\t", index=False)
            print(f"Wrote per-genus stats for {genus} → {per_genus_stats_path}", flush=True)
        else:
            print(f"No sequences found for genus {genus}, skipping per-genus stats.tsv", flush=True)
    
    for gene_can, info in gene_step3_info.items():
        num_new = info["num_matches"]
        if num_new <= 0:
            continue
        
        genus_str = ";".join(sorted(info["genus_set"])) if info["genus_set"] else ""
        source_orgs_str = ";".join(info["source_orgs"]) if info["source_orgs"] else ""
        accessions_str = ";".join(info["assembly_accessions"]) if info["assembly_accessions"] else ""
        paths_str = ";".join(sorted(info["paths"])) if info["paths"] else ""
        
        idx = gene_index_by_canonical.get(gene_can)
        
        if idx is not None:
            prev_total = safe_int(gene_stats_in.loc[idx, "TOTAL_MATCHES"])
            gene_stats_in.loc[idx, "TOTAL_MATCHES"] = prev_total + num_new
            
            prev_step3 = safe_int(gene_stats_in.loc[idx, "STEP-3_NUM_MATCHES"])
            gene_stats_in.loc[idx, "STEP-3_NUM_MATCHES"] = prev_step3 + num_new
            
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
            
            append_field("STEP-3_GENUS", genus_str)
            append_field("STEP-3_SOURCE_ORGS", source_orgs_str)
            append_field("STEP-3_ACCESSIONS", accessions_str)
            append_field("STEP-3_PATHS", paths_str)
        else:
            new_row = {
                "GENE_NAME": gene_can, "TOTAL_MATCHES": num_new,
                "STEP-1_NUM_MATCHES": 0, "STEP-1_ORGS": "", "STEP-1_PATHS": "",
                "STEP-2_NUM_MATCHES": 0, "STEP-2_ORGS": "", "STEP-2_PATHS": "",
                "STEP-3_NUM_MATCHES": num_new, "STEP-3_GENUS": genus_str,
                "STEP-3_SOURCE_ORGS": source_orgs_str, "STEP-3_ACCESSIONS": accessions_str,
                "STEP-3_PATHS": paths_str,
            }
            for col in gene_stats_in.columns:
                if col not in new_row:
                    new_row[col] = "" if gene_stats_in[col].dtype == object else 0
            gene_stats_in = pd.concat([gene_stats_in, pd.DataFrame([new_row])], ignore_index=True)
            gene_index_by_canonical[gene_can] = gene_stats_in.index[-1]
    
    if "GENE_NAME" in gene_stats_in.columns:
        gene_stats_in.sort_values(by="GENE_NAME", inplace=True)
    
    gene_stats_in.to_csv(out_gene_stats, sep="\t", index=False)
    print(f"\nWrote {out_gene_stats}", flush=True)
    
    genus_rows_global = []
    for (genus, gene_can), rinfo in genus_gene_result_info.items():
        fna_path = rinfo.get("fna_path")
        if not fna_path:
            continue
        
        input_info = genus_gene_input_info.get((genus, gene_can), {"species": set(), "accessions": set()})
        input_orgs_str = ";".join(sorted(input_info["species"])) if input_info["species"] else ""
        input_acc_str = ";".join(sorted(input_info["accessions"])) if input_info["accessions"] else ""

        source_orgs = rinfo.get("source_orgs", [])
        asm_accs = rinfo.get("assembly_accessions", [])
        nuc_accs = rinfo.get("nuccore_accessions", [])
    
        source_orgs_str = ";".join(source_orgs) if source_orgs else ""
        source_asm_str = ";".join(asm_accs) if asm_accs else ""
        source_nuc_str = ";".join(nuc_accs) if nuc_accs else ""
    
        genus_rows_global.append([genus, gene_can, input_orgs_str, input_acc_str, source_orgs_str, source_asm_str, source_nuc_str, fna_path])

    genus_df_global = pd.DataFrame(genus_rows_global, columns=[
        "GENUS", "GENE", "INPUT_ORGS", "INPUT_ACCESSIONS",
        "SOURCE_ORGS", "SOURCE_ASSEMBLY_ACCESSIONS", "SOURCE_NUCCORE_ACCESSIONS", "FNA_PATH"
    ])
    genus_df_global.to_csv(out_genus_stats, sep="\t", index=False)
    print(f"Wrote {out_genus_stats}", flush=True)

    print("\nStep 3 completed ✅", flush=True)
    return True
