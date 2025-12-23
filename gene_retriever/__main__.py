#!/usr/bin/env python3
"""
Main pipeline runner for AMR gene sequence retrieval.

Usage:
    python3 main.py --step 1 --input input.tsv --refseq /path/REFSEQ --genedir /path/GENE_SEQUENCES
    python3 main.py --step 2 --genedir /path/GENE_SEQUENCES
    python3 main.py --step 3 --genedir /path/GENE_SEQUENCES
    python3 main.py --step 4 --genedir /path/GENE_SEQUENCES --genes-not-found GENES_NOT_FOUND.txt
    python3 main.py --step all --input input.tsv --refseq /path/REFSEQ --genedir /path/GENE_SEQUENCES
"""

import sys
import argparse
import os
import shutil

def main():
    parser = argparse.ArgumentParser(
        description='AMR Gene Sequence Retrieval Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Run Step 1:
    python3 gene_finder --step 1 --input input.tsv --refseq ./REFSEQ_DATASET --genedir ./GENE_SEQUENCES
  
  Run Step 2:
    python3 gene_finder --step 2 --genedir ./GENE_SEQUENCES
  
  Run Step 3:
    python3 gene_finder --step 3 --genedir ./GENE_SEQUENCES
  
  Run Step 4:
    python3 gene_finder --step 4 --genedir ./GENE_SEQUENCES --genes-not-found GENES_NOT_FOUND.txt
  
  Run all steps:
    python3 gene_finder --step all --input input.tsv --refseq ./REFSEQ_DATASET --genedir ./GENE_SEQUENCES
        """
    )
    
    parser.add_argument('--step', required=True, choices=['1', '2', '3', '4', 'all'],
                        help='Which step to run (1, 2, 3, 4, or all)')
    parser.add_argument('--input', help='Input TSV file (required for step 1)')
    parser.add_argument('--refseq', help='RefSeq dataset directory (required for step 1)')
    parser.add_argument('--genedir', required=True, help='Gene sequences output directory')
    parser.add_argument('--run-stats', help='*_run_stats.tsv file (required for steps 2 and 3)')
    parser.add_argument('--gene-stats', help='*_gene_stats_found.tsv file (required for steps 2, 3 and 4)')
    parser.add_argument('--genes-not-found', help='GENES_NOT_FOUND.txt file (required for step 4)')
    
    args = parser.parse_args()

    ############# Setting up runner args for each step

    if not args.genedir:
            parser.error("--genedir is required for step 1")
            sys.exit(1)
    geneseq_basedir = args.genedir
    ############ Create gene directory if it doesn't exist
    os.makedirs(args.genedir, exist_ok=True)
    
    if args.step in ['1', 'all']:
        if not ( args.input and os.path.exists(args.input) ):
            parser.error("--input is required for step 1")
            sys.exit(1)
        input_amr_genes = args.input
        if not args.refseq:
            parser.error("--refseq is required for step 1")
            sys.exit(1)
        refseq_basedir = args.refseq
        os.makedirs(args.refseq, exist_ok=True)        
    elif args.step in ['2', '3']:
        if not (args.run_stats and os.path.exists(args.run_stats)):
            if os.path.exists(str(f"{int(args.step)-1}_run_stats.tsv")):
                in_run_stats = str(f"{int(args.step)-1}_run_stats.tsv")
            else:
                parser.error(f"--run-stats is required for steps 2 & 3. Run step {int(args.step)-1} first.")
                sys.exit(1)
        else:
            in_run_stats = args.run_stats
        if not (args.gene_stats and os.path.exists(args.gene_stats)):
            if os.path.exists(str(f"{int(args.step)-1}_gene_stats_found.tsv")):
                in_gene_stats = str(f"{int(args.step)-1}_gene_stats_found.tsv")
            else:
                parser.error(f"--gene-stats is required for steps 2 & 3. Run step {int(args.step)-1} first.")
                sys.exit(1)
        else:
            in_gene_stats = args.gene_stats
    elif args.step == '4':
        if not (args.gene_stats and os.path.exists(args.gene_stats)):
            if os.path.exists(str(f"{int(args.step)-1}_gene_stats_found.tsv")):
                in_gene_stats = str(f"{int(args.step)-1}_gene_stats_found.tsv")
            else:
                parser.error(f"--gene-stats is required for step 4. Run step {int(args.step)-1} first.")
        else:
            in_gene_stats = args.gene_stats
        if os.path.exists(args.run_stats):
            src = os.path.abspath(args.run_stats)
            dst = os.path.abspath("3_run_stats.tsv")
            if src != dst:
                in_run_stats = shutil.move(src, dst)
            else:
                in_run_stats = "3_run_stats.tsv"
        if not ( args.genes_not_found and os.path.exists(args.genes_not_found) ):
            if not ( in_gene_stats and in_run_stats ):
                parser.error(f"--genes-not-found is required for step 4. Run step {int(args.step)-1} first. Or provide --gene-stats and --run-stats to generate GENES_NOT_FOUND.txt.")
                sys.exit(1)
            else:                
                cmd = 'cut -f6 3_run_stats.tsv | tr "," "\\n" | grep -v "^GENES_NOT_FOUND$" | grep -v "^$" | python3 utils/normalize_genes.py | awk \'FNR==NR{a[$1]=1; next} !($1 in a)\' 3_gene_stats_found.tsv - > GENES_NOT_FOUND.txt'
                os.system(cmd)
                if not os.path.exists("GENES_NOT_FOUND.txt"):
                    print("ERROR: Failed to generate GENES_NOT_FOUND.txt", flush=True)
                    sys.exit(1)
                genes_not_found = "GENES_NOT_FOUND.txt"
        else:
             genes_not_found = args.genes_not_found

    
    success = True
    
    ############ Run requested step(s)
    if args.step == '1' or args.step == 'all':
        from step1.runner import run_step1
        success = run_step1(input_amr_genes, refseq_basedir, geneseq_basedir)
        if not success:
            print("\n❌ Step 1 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
    
    if args.step == '2' or args.step == 'all':
        from step2.runner import run_step2
        success = run_step2(in_run_stats, in_gene_stats, geneseq_basedir)
        if not success:
            print("\n❌ Step 2 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
    
    if args.step == '3' or args.step == 'all':
        from step3.runner import run_step3
        success = run_step2(in_run_stats, in_gene_stats, geneseq_basedir)
        if not success:
            print("\n❌ Step 3 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
    
    if args.step == '4' or args.step == 'all':
        from step4.runner import run_step4        
        success = run_step4(in_gene_stats, genes_not_found, geneseq_basedir)
        if not success:
            print("\n❌ Step 4 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
    
    print("\n" + "="*80)
    if args.step == 'all':
        print("✅ ALL STEPS COMPLETED SUCCESSFULLY")
    else:
        print(f"✅ STEP {args.step} COMPLETED SUCCESSFULLY")
    print("="*80)

if __name__ == "__main__":
    main()
