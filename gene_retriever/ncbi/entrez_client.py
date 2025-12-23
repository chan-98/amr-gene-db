#!/usr/bin/env python3
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
import time
import traceback
from genes.normalisation import normalise_hyphens, generate_gene_patterns

Entrez.email = "chandini@acrannolife.com"
LAST_CALL = 0.0

def rate_limited_sleep(min_interval=0.4):
    global LAST_CALL
    now = time.time()
    if LAST_CALL == 0.0:
        LAST_CALL = now
        return
    elapsed = now - LAST_CALL
    if elapsed < min_interval:
        time.sleep(min_interval - elapsed)
    LAST_CALL = time.time()

def query_ncbi(query, db="nucleotide", max_results=10):
    print(f"----------------------------------------------SEARCHING FOR: {query}-----------------------------------------------", flush=True)
    for attempt in range(3):
        try:
            rate_limited_sleep()
            search_handle = Entrez.esearch(db=db, term=query, retmax=max_results, sort="relevance")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            return search_results.get("IdList", [])
        except HTTPError as e:
            print(f"HTTPError: {e}. Retrying ({attempt+1}/3)...", flush=True)
            time.sleep(3 * (attempt + 1))
        except RuntimeError as e:
            print(f"NCBI search backend error for '{query}': {e}. Skipping.", flush=True)
            return []
        except Exception as e:
            print(f"Unexpected error during NCBI search for '{query}': {e}", flush=True)
            time.sleep(3)
    return []

def efetch_ncbi(result_id, return_type, ncbi_db="nucleotide", max_results=1):
    print(f"Fetching {return_type} for ID {result_id} ...", flush=True)
    for attempt in range(3):
        try:
            rate_limited_sleep()
            fetch_handle = Entrez.efetch(db=ncbi_db, id=result_id, rettype=return_type, retmax=max_results)
            record_format = "genbank" if return_type == "gb" else "fasta"
            record = SeqIO.read(fetch_handle, record_format)
            fetch_handle.close()
            return record
        except HTTPError as e:
            print(f"HTTPError during efetch for ID {result_id}: {e}. Retrying ({attempt+1}/3)...", flush=True)
            time.sleep(3 * (attempt + 1))
        except Exception as e:
            print(f"Error fetching {return_type} for ID {result_id}: {e}", flush=True)
            return None
    return None

def get_assembly_accession_from_nuccore(nuccore_id: str) -> str:
    try:
        rate_limited_sleep()
        link_handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuccore_id, rettype="xml")
        link_record = Entrez.read(link_handle)
        link_handle.close()
        linksets = link_record[0].get("LinkSetDb", [])
        if not linksets:
            return "NA"
        links = linksets[0].get("Link", [])
        if not links:
            return "NA"
        asm_uid = links[0]["Id"]
        rate_limited_sleep()
        sum_handle = Entrez.esummary(db="assembly", id=asm_uid, report="full")
        sum_record = Entrez.read(sum_handle)
        sum_handle.close()
        docsum = sum_record["DocumentSummarySet"]["DocumentSummary"][0]
        asm_acc = docsum.get("AssemblyAccession", None)
        if asm_acc:
            return asm_acc
        return "NA"
    except Exception as e:
        print(f"Warning: could not resolve assembly accession for {nuccore_id}: {e}", flush=True)
        return "NA"

def fetch_gene_from_species(gene_symbol, genus, species, max_results=20):
    patterns = generate_gene_patterns(gene_symbol)
    id_list = []
    for p in patterns:
        query = f"{p}[Gene] AND ({genus} {species}[Organism])"
        id_list = query_ncbi(query, "nucleotide", max_results)
        if id_list:
            print(f"Found {len(id_list)} sequences.", flush=True)
            break
    if not id_list:
        print("No search results found for the query.", flush=True)
        return []
    cds_recs = []
    for i, idx in enumerate(id_list, 1):
        print(f"Fetching result {i}/{len(id_list)} (ID: {idx})", flush=True)
        gb_record = efetch_ncbi(idx, "gb", "nucleotide", 10)
        if not gb_record:
            print(f"No genbank record for UID {idx}", flush=True)
            continue
        print(f"Found genbank record for UID {idx}: {gb_record.description}", flush=True)
        fasta_record = efetch_ncbi(idx, "fasta", "nucleotide", 10)
        if not fasta_record:
            print("No sequences found for this entry.", flush=True)
            continue
        seq_count = 0
        for rec in gb_record.features:
            if rec.type == "CDS" and any(gene_symbol == g for g in rec.qualifiers.get("gene", [])):
                try:
                    cds_seq = rec.extract(fasta_record.seq)
                    print(f"=== FOUND {gene_symbol} at {rec.location} of {gb_record.description} -- (ID: {gb_record.id}; Protein ID: {rec.qualifiers.get('protein_id', [''])[0]}) ===", flush=True)
                    seq_count += 1
                    header = (
                        f">{gb_record.id} [gene={gene_symbol}] "
                        f"[locus_tag={rec.qualifiers.get('locus_tag', [''])[0]}] "
                        f"[protein={rec.qualifiers.get('product', [''])[0]}] "
                        f"[protein_id={rec.qualifiers.get('protein_id', [''])[0]}] "
                        f"[location={rec.location}] [description={gb_record.description}]"
                    )
                    cds_rec = f"{header}\n{cds_seq}"
                    cds_recs.append(cds_rec)
                except Exception as e:
                    print(f"Error extracting CDS: {e}", flush=True)
        if seq_count == 0:
            print(f"!!! CDS sequence for {gene_symbol} not found for UID {idx}.", flush=True)
    return cds_recs

def fetch_gene_for_genus(gene_symbol_raw, genus, max_results=10):
    patterns = generate_gene_patterns(gene_symbol_raw)
    uid_list = []
    for pat in patterns:
        query = f"{pat}[Gene] AND {genus}[Organism]"
        ids = query_ncbi(query, "nucleotide", max_results)
        if ids:
            print(f"Found {len(ids)} hits for pattern '{pat}', stopping further pattern checks.", flush=True)
            uid_list = ids
            break
        else:
            print(f"No hits for pattern '{pat}' in genus {genus}.", flush=True)
    if not uid_list:
        print(f"No search results found for {gene_symbol_raw} in genus {genus}.", flush=True)
        return [], []
    print(f"Using {len(uid_list)} UIDs for further efetch for gene '{gene_symbol_raw}' in genus '{genus}'.", flush=True)
    cds_recs = []
    seq_meta = []
    for i, idx in enumerate(uid_list, 1):
        print(f"Fetching result {i}/{len(uid_list)} (ID: {idx})", flush=True)
        gb_record = efetch_ncbi(idx, "gb", "nucleotide", 1)
        if not gb_record:
            print(f"No GenBank record for UID {idx}", flush=True)
            continue
        fasta_record = efetch_ncbi(idx, "fasta", "nucleotide", 1)
        if not fasta_record:
            print("No FASTA sequence found for this entry.", flush=True)
            continue
        organism = gb_record.annotations.get("organism", gb_record.description)
        strain = None
        for feat in gb_record.features:
            if feat.type == "source":
                q = feat.qualifiers
                strain = q.get("strain", [None])[0] or q.get("isolate", [None])[0] or q.get("serovar", [None])[0]
                break
        organism_full = f"{organism} {strain}" if strain else organism
        nuccore_id = gb_record.id
        assembly_acc = get_assembly_accession_from_nuccore(nuccore_id)
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
                print(f"=== FOUND {gene_symbol_raw} at {rec.location} of {gb_record.description} -- (Nuccore: {nuccore_id}, Assembly: {assembly_acc}, Gene Symbol in Record: {rec.qualifiers.get('gene', [gene_symbol_raw])[0]}) ===", flush=True)
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
                    "assembly_accession": assembly_acc,
                    "nuccore_accession": nuccore_id,
                    "fna_gene_name": file_gene_name
                })
            except Exception as e:
                print(f"Error extracting CDS: {e}", flush=True)
        if per_uid_seq_count == 0:
            print(f"!!! CDS sequence for {gene_symbol_raw} not found in UID {idx}.", flush=True)
    return cds_recs, seq_meta
