#!/usr/bin/env python3
import re
import os
from genes.normalisation import (
    normalise_hyphens,
    normalise_apostrophes,
    strip_wrapping_quotes,
    clean_gene_term,
    generate_gene_patterns
)

def sanitize_filename(name: str) -> str:
    name = strip_wrapping_quotes(name)
    name = normalise_hyphens(name)
    name = normalise_apostrophes(name)
    name = re.sub(r'[\/\\\s:"\*?<>\|]+', '_', name)
    name = re.sub(r'_+', '_', name).strip('_')
    return name

def parse_fasta_headers(cds_file):
    gene_records = {}
    current_header = None
    current_sequence = []
    current_gene_lower = None
    current_gene_original = None
    
    with open(cds_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_gene_lower and current_header:
                    seq_single_line = ''.join(current_sequence).replace('\n', '').replace('\r', '')
                    gene_records[current_gene_lower] = (current_header, seq_single_line, current_gene_original)
                current_header = line.rstrip('\n\r')
                current_sequence = []
                current_gene_lower = None
                current_gene_original = None
                gene_match = re.search(r'\[gene=([^\]]+)\]', line)
                if gene_match:
                    current_gene_original = strip_wrapping_quotes(gene_match.group(1))
                    key = normalise_hyphens(current_gene_original)
                    key = normalise_apostrophes(key).lower()
                    current_gene_lower = key
            else:
                if current_header:
                    current_sequence.append(line.rstrip('\n\r'))
        if current_gene_lower and current_header:
            seq_single_line = ''.join(current_sequence).replace('\n', '').replace('\r', '')
            gene_records[current_gene_lower] = (current_header, seq_single_line, current_gene_original)
    return gene_records

def search_gene_in_cds(query_gene, cds_file):
    cleaned_gene = clean_gene_term(query_gene)
    patterns = generate_gene_patterns(cleaned_gene)
    gene_records = parse_fasta_headers(cds_file)
    for pattern in patterns:
        if pattern in gene_records:
            header, sequence, original_gene_name = gene_records[pattern]
            return (header, sequence, original_gene_name)
    return None

def load_gene_stats(gene_stats_file):
    gene_stats = {}
    if os.path.exists(gene_stats_file):
        with open(gene_stats_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 5:
                    gene = parts[0]
                    count = int(parts[2])
                    organisms = parts[3].split(',') if parts[3] else []
                    paths = parts[4].split(',') if parts[4] else []
                    gene_stats[gene] = {'count': count, 'organisms': organisms, 'paths': paths}
    return gene_stats

def save_gene_stats(gene_stats_file, gene_stats):
    with open(gene_stats_file, 'w') as f:
        f.write("GENE_NAME\tTOTAL_MATCHES\tSTEP-1_NUM_MATCHES\tSTEP-1_ORGS\tSTEP-1_PATHS\n")
        for gene in sorted(gene_stats.keys()):
            count = gene_stats[gene]['count']
            organisms = ','.join(gene_stats[gene]['organisms'])
            paths = ','.join(gene_stats[gene]['paths'])
            f.write(f"{gene}\t{count}\t{count}\t{organisms}\t{paths}\n")
