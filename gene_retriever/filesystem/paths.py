import os
import re
from genes.normalisation import (
    canonical_gene_name,
    generate_gene_patterns,
    strip_wrapping_quotes,
    normalise_hyphens,
    normalise_apostrophes
)

def sanitize_filename(name: str) -> str:
    if not isinstance(name, str):
        name = str(name)
    name = name.replace('"', '')
    name = re.sub(r'[\/\\\s:"\*\?<>\|]+', '_', name)
    name = re.sub(r'_+', '_', name)
    name = name.strip('_')
    return name

# Step 2: Species-based paths
def build_output_dir_for_accession(accession, genus, species, outdir):
    return os.path.join(outdir, f"{accession}_{genus}_{species}", "2_QUERY_GENE_GENUS_SPECIES")

def generate_filename_candidates(gene_symbol):
    candidates = []
    canonical = sanitize_filename(canonical_gene_name(gene_symbol))
    if canonical:
        candidates.append(canonical)
    normalised = sanitize_filename(
        normalise_apostrophes(normalise_hyphens(strip_wrapping_quotes(gene_symbol)))
    )
    if normalised and normalised not in candidates:
        candidates.append(normalised)
    for pattern in generate_gene_patterns(gene_symbol):
        patt_safe = sanitize_filename(pattern)
        if patt_safe and patt_safe not in candidates:
            candidates.append(patt_safe)
    return candidates

def find_existing_fasta_for_gene(gene_symbol, accession, genus, species, outdir):
    found_paths = []
    target_dir = build_output_dir_for_accession(accession, genus, species, outdir)
    if not os.path.isdir(target_dir):
        print(f"Directory does not exist for this organism: {target_dir}. Gene needs to be retrieved. Continuing...", flush=True)
        return found_paths
    filenames = generate_filename_candidates(gene_symbol)
    print(f"Checking {len(filenames)} filename variants for {gene_symbol}: {filenames}")
    for fname in filenames:
        candidate = os.path.join(target_dir, f"{fname}.fna")
        if os.path.isfile(candidate):
            found_paths.append(os.path.abspath(candidate))
        else:
            print(f"Did not find any existing files stored for {fname}. Gene needs to be retrieved. Continuing...", flush=True)
    return found_paths

def save_fasta_step2(records, gene, genus, species, accession, outdir):
    if not records:
        return None
    dir_path = build_output_dir_for_accession(accession, genus, species, outdir)
    os.makedirs(dir_path, exist_ok=True)
    safe_name = sanitize_filename(gene)
    cds_path = os.path.join(dir_path, f"{safe_name}.fna")
    with open(cds_path, "w") as f:
        for rec in records:
            f.write(f"{rec}\n")
    return os.path.abspath(cds_path)

# Step 3: Genus-based paths
def save_fasta_step3(records, gene_canonical, genus, outdir):
    if not records:
        return None
    dir_path = os.path.join(outdir, genus)
    os.makedirs(dir_path, exist_ok=True)
    cds_path = os.path.join(dir_path, f"{gene_canonical}.fna")
    with open(cds_path, "w") as f:
        for rec in records:
            f.write(f"{rec}\n")
    return os.path.abspath(cds_path)

# Step 4: Organism-specific paths
def save_fasta_step4(records, gene_safe, accession, genus, species, outdir):
    if not records:
        return None
    dir_path = os.path.join(outdir, f"{accession}_{genus}_{species}", "4_QUERY_GENE")
    os.makedirs(dir_path, exist_ok=True)
    cds_path = os.path.join(dir_path, f"{gene_safe}.fna")
    with open(cds_path, "w") as f:
        for rec in records:
            f.write(f"{rec}\n")
    return os.path.abspath(cds_path)
