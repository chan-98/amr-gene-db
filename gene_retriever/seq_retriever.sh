#!/bin/bash

# Usage: ./sequence_retrieval.sh input.tsv /path/REFSEQ_DATASET /path/GENE_SEQUENCES

input_file=$1
refseq_basedir=$2
geneseq_basedir=$3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

shopt -s nullglob
# Create necessary directories
mkdir -p "${refseq_basedir}" "${geneseq_basedir}"

# Initialize run stats file
run_stats_file="1_run_stats.tsv"
echo -e "GENUS\tSPECIES\tACCESSION_STATUS\tINPUT_GENES\tGENES_FOUND\tGENES_NOT_FOUND" > "$run_stats_file"

# Read the input file line by line (skip header)
tail -n +2 "$input_file" | while IFS=$'\t' read -r genus species genes_raw; do
    # Clean up carriage returns and whitespace
    genus=$(echo "$genus" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    species=$(echo "$species" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    genes_raw=$(echo "$genes_raw" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    
    echo "========================================"
    echo "Processing: $genus $species"
    echo "Raw genes input: $genes_raw"
    
    # Check if genes field is empty, nil, or NA
    if [[ -z "$genes_raw" ]] || [[ "$genes_raw" == "nil" ]] || [[ "$genes_raw" == "NIL" ]] || [[ "$genes_raw" == "Nil" ]] || [[ "$genes_raw" == "NA" ]] || [[ "$genes_raw" == "na" ]]; then
        genes_status="NIL"
        echo "No AMR genes specified for this organism."
    else
        genes_status="$genes_raw"
    fi
    
    # Check if dataset already exists
    dataset_found=""
    acc=""
    
    # Search for existing dataset directory matching genus and species
    for dir in "${refseq_basedir}"/*_"${genus}"_"${species}"; do
        if [[ -d "$dir" ]]; then
            dataset_found="$dir"
            # Extract accession from directory name (e.g., GCF_029866465.1 from GCF_029866465.1_Genus_species)
            dirname=$(basename "$dir")
            acc="${dirname%%_${genus}_${species}}"
            echo "Found existing dataset: $dir (Accession: $acc)"
            break
        fi
    done
    
    # If dataset not found, download it
    if [[ -z "$dataset_found" ]]; then
        echo "Dataset not found. Fetching accession for $genus $species..."
        acc=$(datasets summary genome taxon "$genus $species" --reference --assembly-source RefSeq --limit 1 2>/dev/null | jq -r '.reports[0].accession // "ORGANISM_NOT_FOUND"')
        
        if [[ "$acc" == "ORGANISM_NOT_FOUND" ]] || [[ -z "$acc" ]]; then
            echo "ERROR: Organism not found in NCBI: $genus $species"
            echo -e "$genus\t$species\tORGANISM_NOT_FOUND\t$genes_status\t\t" >> "$run_stats_file"
            echo
            continue
        fi
        
        echo "Accession found: $acc"
        echo "Downloading reference genome..."
        
        datasets download genome accession "$acc" --include gff3,gtf,cds,genome --filename "${refseq_basedir}/${acc}_ncbi_dataset.zip"
        
        if [[ $? -eq 0 ]]; then
            echo "Extracting dataset..."
            unzip -q "${refseq_basedir}/${acc}_ncbi_dataset.zip"
            mv "ncbi_dataset/data/${acc}" "${refseq_basedir}/${acc}_${genus}_${species}"
            rm -rf ncbi_dataset md5sum.txt README.md "${refseq_basedir}/${acc}_ncbi_dataset.zip"
            dataset_found="${refseq_basedir}/${acc}_${genus}_${species}"
            echo "Dataset downloaded and extracted successfully. Location: ${dataset_found}"
        else
            echo "ERROR: Failed to download dataset for $genus $species"
            echo -e "$genus\t$species\tDOWNLOAD_FAILED\t$genes_status\t\t" >> "$run_stats_file"
            echo
            continue
        fi
        
        sleep 5
    fi
    
    # Now process genes if they exist
    if [[ "$genes_status" != "NIL" ]]; then
        # Check if CDS file exists
        cds_file="${dataset_found}/cds_from_genomic.fna"
        
        if [[ ! -f "$cds_file" ]]; then
            echo "ERROR: CDS file not found: $cds_file"
            echo -e "$genus\t$species\t$acc\t$genes_status\t\tCDS_FILE_NOT_FOUND" >> "$run_stats_file"
            echo
            continue
        fi
        
        echo "CDS file found: $cds_file"
        echo "Searching for genes using Python script..."
        
        # Call Python script to search for genes
        python3 "${SCRIPT_DIR}/gene_finder.py" "$genus" "$species" "$acc" "$cds_file" "$genes_raw" "$run_stats_file" "$geneseq_basedir"
        
    else
        # No genes to search for
        echo "No genes to search for this organism."
        echo -e "$genus\t$species\t$acc\tNIL\t\t" >> "$run_stats_file"
    fi
    
    echo
done

echo "========================================"
echo "Processing complete!"
echo "Results saved in: $run_stats_file"
