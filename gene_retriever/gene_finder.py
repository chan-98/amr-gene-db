import sys
import os
from pathlib import Path
from step1.cds_parser import (
    search_gene_in_cds,
    sanitize_filename,
    load_gene_stats,
    save_gene_stats
)

def parse_genes_list(genes_raw):
    """Parse comma-separated genes list."""
    genes = [g.strip() for g in genes_raw.split(',')]
    return [g for g in genes if g]

def main():
    if len(sys.argv) != 8:
        print("Usage: gene_finder.py <genus> <species> <accession> <cds_file> <genes_raw> <run_stats_file> <gene_seq_basedir>")
        sys.exit(1)
    
    genus = sys.argv[1]
    species = sys.argv[2]
    accession = sys.argv[3]
    cds_file = sys.argv[4]
    genes_raw = sys.argv[5]
    run_stats_file = sys.argv[6]
    gene_seq_basedir = sys.argv[-1]
    
    gene_stats_file = os.path.join(os.path.dirname(run_stats_file) or '.', '1_gene_stats_found.tsv')
    
    print(f"Searching for genes in {genus} {species} (Accession: {accession})")
    
    genes_list = parse_genes_list(genes_raw)
    print(f"Genes to search: {genes_list}")
    
    genes_found = []
    genes_not_found = []
    
    output_dir = None
    current_dir = os.getcwd()
    
    gene_stats = load_gene_stats(gene_stats_file)
    organism_name = f"{genus} {species}"
    
    for gene in genes_list:
        print(f"  Searching for gene: {gene}")
        
        result = search_gene_in_cds(gene, cds_file)
        
        if result:
            header, sequence, matched_gene_name = result
            print(f"    ✓ FOUND: {gene} (matched as: {matched_gene_name})")
            genes_found.append(matched_gene_name)
            
            if output_dir is None:
                output_dir = f"{gene_seq_basedir}/{accession}_{genus}_{species}/1_EXACT_MATCH_FROM_CDS"
                Path(output_dir).mkdir(parents=True, exist_ok=True)
                print(f"    Created output directory: {output_dir}")
            
            safe_name = sanitize_filename(matched_gene_name)
            output_file = f"{output_dir}/{safe_name}.fna"
            with open(output_file, 'w') as f:
                f.write(f"{header}\n{sequence}\n")
            
            abs_path = os.path.join(current_dir, output_file)
            print(f"    Saved to: {output_file}")
            
            if matched_gene_name in gene_stats:
                gene_stats[matched_gene_name]['count'] += 1
                gene_stats[matched_gene_name]['organisms'].append(organism_name)
                gene_stats[matched_gene_name]['paths'].append(abs_path)
            else:
                gene_stats[matched_gene_name] = {
                    'count': 1,
                    'organisms': [organism_name],
                    'paths': [abs_path]
                }
        else:
            print(f"    ✗ NOT FOUND: {gene}")
            genes_not_found.append(gene)
    
    genes_found_str = ','.join(genes_found) if genes_found else ""
    genes_not_found_str = ','.join(genes_not_found) if genes_not_found else ""
    
    print(f"\nSummary:")
    print(f"  Total genes: {len(genes_list)}")
    print(f"  Found: {len(genes_found)}")
    print(f"  Not found: {len(genes_not_found)}")
    
    with open(run_stats_file, 'a') as f:
        f.write(f"{genus}\t{species}\t{accession}\t{genes_raw}\t{genes_found_str}\t{genes_not_found_str}\n")
    
    save_gene_stats(gene_stats_file, gene_stats)
    
    print(f"Results appended to: {run_stats_file}")
    print(f"Gene statistics updated in: {gene_stats_file}")

if __name__ == "__main__":
    main()
