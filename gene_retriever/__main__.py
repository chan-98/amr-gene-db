#!/usr/bin/env python3
"""
Main pipeline runner for AMR gene sequence retrieval.

Usage:
    python3 gene_finder --step 1 --input input.tsv --refseq /path/REFSEQ --genedir /path/GENE_SEQUENCES
    python3 gene_finder --step 2,3 --genedir /path/GENE_SEQUENCES
    python3 gene_finder --step 1,2,3 --input input.tsv --refseq /path/REFSEQ --genedir /path/GENE_SEQUENCES
    python3 gene_finder --step all --input input.tsv --refseq /path/REFSEQ --genedir /path/GENE_SEQUENCES
"""

import sys
import argparse
import os
import shutil

def parse_steps(step_arg):
    """
    Parse step argument and validate it's a continuous range.
    Returns list of steps to run.
    """
    if step_arg == 'all':
        return [1, 2, 3, 4]
    
    # Parse comma-separated steps
    try:
        steps = [int(s.strip()) for s in step_arg.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid step format: {step_arg}")
    
    # Validate steps are 1-4
    if not all(1 <= s <= 4 for s in steps):
        raise argparse.ArgumentTypeError("Steps must be between 1 and 4")
    
    # Check for duplicates
    if len(steps) != len(set(steps)):
        raise argparse.ArgumentTypeError("Duplicate steps found")
    
    # Sort steps
    steps.sort()
    
    # Validate continuous range
    if len(steps) > 1:
        for i in range(len(steps) - 1):
            if steps[i+1] - steps[i] != 1:
                raise argparse.ArgumentTypeError(
                    f"Steps must be continuous (e.g., 1,2,3 or 2,3 or 3,4). "
                    f"Invalid: {step_arg}"
                )
    
    return steps

def main():
    parser = argparse.ArgumentParser(
        description='AMR Gene Sequence Retrieval Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Run single step:
    python3 gene_finder --step 1 --input input.tsv --refseq ./REFSEQ_DATASET --genedir ./GENE_SEQUENCES
  
  Run continuous steps:
    python3 gene_finder --step 1,2,3 --input input.tsv --refseq ./REFSEQ_DATASET --genedir ./GENE_SEQUENCES
    python3 gene_finder --step 2,3 --genedir ./GENE_SEQUENCES
    python3 gene_finder --step 3,4 --genedir ./GENE_SEQUENCES
  
  Run all steps:
    python3 gene_finder --step all --input input.tsv --refseq ./REFSEQ_DATASET --genedir ./GENE_SEQUENCES
  
  Invalid (non-continuous):
    python3 gene_finder --step 1,3 --genedir ./GENE_SEQUENCES  ❌ (skips step 2)
    python3 gene_finder --step 2,4 --genedir ./GENE_SEQUENCES  ❌ (skips step 3)
        """
    )
    
    parser.add_argument('--step', required=True,
                        help='Which step(s) to run: single (1, 2, 3, 4), continuous range (1,2 or 2,3,4), or all')
    parser.add_argument('--input', help='Input TSV file (required for step 1)')
    parser.add_argument('--refseq', help='RefSeq dataset directory (required for step 1)')
    parser.add_argument('--genedir', required=True, help='Gene sequences output directory')
    parser.add_argument('--run-stats', help='*_run_stats.tsv file (optional for steps 2 and 3)')
    parser.add_argument('--gene-stats', help='*_gene_stats_found.tsv file (optional for steps 2, 3 and 4)')
    parser.add_argument('--genes-not-found', help='GENES_NOT_FOUND.txt file (optional for step 4)')
    
    args = parser.parse_args()

    # Parse and validate steps
    try:
        steps_to_run = parse_steps(args.step)
    except argparse.ArgumentTypeError as e:
        parser.error(str(e))
    
    print(f"\n{'='*80}")
    print(f"Running steps: {', '.join(map(str, steps_to_run))}")
    print(f"{'='*80}\n")

    ############# Setting up runner args for each step

    if not args.genedir:
        parser.error("--genedir is required")
        sys.exit(1)
    geneseq_basedir = args.genedir
    
    ############ Create gene directory if it doesn't exist
    os.makedirs(args.genedir, exist_ok=True)
    
    # Variables to track across steps
    input_amr_genes = None
    refseq_basedir = None
    in_run_stats = None
    in_gene_stats = None
    genes_not_found = None
    success = True
    
    # Validate step 1 requirements
    if 1 in steps_to_run:
        if not (args.input and os.path.exists(args.input)):
            parser.error("--input is required for step 1")
            sys.exit(1)
        input_amr_genes = args.input
        if not args.refseq:
            parser.error("--refseq is required for step 1")
            sys.exit(1)
        refseq_basedir = args.refseq
        os.makedirs(args.refseq, exist_ok=True)
    
        # Run step 1
        from step1.runner import run_step1
        print(f"\n{'='*80}")
        print("EXECUTING STEP 1")
        print(f"{'='*80}")
        print(f"""
        Arguments:
        \t--input:\t{input_amr_genes}
        \t--refseq:\t{refseq_basedir}
        \t--genedir:\t{geneseq_basedir}
        """)
        success = run_step1(input_amr_genes, refseq_basedir, geneseq_basedir)
        if not success:
            print("\n❌ Step 1 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
        # Update variables for next step
        in_run_stats = "1_run_stats.tsv"
        in_gene_stats = "1_gene_stats_found.tsv"
    
    # Validate step 2 requirements
    if 2 in steps_to_run:
        # Try provided arg, then look for file from previous step
        if args.run_stats and os.path.exists(args.run_stats):
            in_run_stats = args.run_stats
        elif os.path.exists("1_run_stats.tsv"):
            in_run_stats = "1_run_stats.tsv"
        else:
            parser.error("--run-stats is required for step 2, or run step 1 first.")
            sys.exit(1)
        
        if args.gene_stats and os.path.exists(args.gene_stats):
            in_gene_stats = args.gene_stats
        elif os.path.exists("1_gene_stats_found.tsv"):
            in_gene_stats = "1_gene_stats_found.tsv"
        else:
            parser.error("--gene-stats is required for step 2, or run step 1 first.")
            sys.exit(1)
        
        # Run step 2
        from step2.runner import run_step2
        print(f"\n{'='*80}")
        print("EXECUTING STEP 2")
        print(f"{'='*80}")
        print(f"""
        Arguments:
        \t--run-stats:\t{in_run_stats}
        \t--gene-stats:\t{in_gene_stats}
        \t--genedir:\t{geneseq_basedir}
        """)
        success = run_step2(in_run_stats, in_gene_stats, geneseq_basedir)
        if not success:
            print("\n❌ Step 2 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
        # Update variables for next step
        in_run_stats = "2_run_stats.tsv"
        in_gene_stats = "2_gene_stats_found.tsv"
    
    # Validate step 3 requirements
    if 3 in steps_to_run:
        # Try provided arg, then look for file from previous step
        if args.run_stats and os.path.exists(args.run_stats):
            in_run_stats = args.run_stats
        elif os.path.exists("2_run_stats.tsv"):
            in_run_stats = "2_run_stats.tsv"
        else:
            parser.error("--run-stats is required for step 3, or run step 2 first.")
            sys.exit(1)
        
        if args.gene_stats and os.path.exists(args.gene_stats):
            in_gene_stats = args.gene_stats
        elif os.path.exists("2_gene_stats_found.tsv"):
            in_gene_stats = "2_gene_stats_found.tsv"
        else:
            parser.error("--gene-stats is required for step 3, or run step 2 first.")
            sys.exit(1)
        
        # Run step 3
        from step3.runner import run_step3
        print(f"\n{'='*80}")
        print("EXECUTING STEP 3")
        print(f"{'='*80}")
        print(f"""
        Arguments:
        \t--run-stats:\t{in_run_stats}
        \t--gene-stats:\t{in_gene_stats}
        \t--genedir:\t{geneseq_basedir}
        """)
        success = run_step3(in_run_stats, in_gene_stats, geneseq_basedir)
        if not success:
            print("\n❌ Step 3 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
        # Update variables for next step
        in_run_stats = "3_run_stats.tsv"
        in_gene_stats = "3_gene_stats_found.tsv"
    
    # Validate step 4 requirements
    if 4 in steps_to_run:
        if args.gene_stats and os.path.exists(args.gene_stats):
            in_gene_stats = args.gene_stats
        elif os.path.exists("3_gene_stats_found.tsv"):
            in_gene_stats = "3_gene_stats_found.tsv"
        else:
            parser.error("--gene-stats is required for step 4, or run step 3 first.")
            sys.exit(1)
        
        # Handle run stats for step 4 (needed for GENES_NOT_FOUND generation)
        if args.run_stats and os.path.exists(args.run_stats):
            src = os.path.abspath(args.run_stats)
            dst = os.path.abspath("3_run_stats.tsv")
            if src != dst:
                shutil.copy(src, dst)
            in_run_stats = "3_run_stats.tsv"
        elif os.path.exists("3_run_stats.tsv"):
            in_run_stats = "3_run_stats.tsv"
        
        # Handle GENES_NOT_FOUND.txt
        if args.genes_not_found and os.path.exists(args.genes_not_found):
            genes_not_found = args.genes_not_found
        elif in_gene_stats and in_run_stats:
            # Auto-generate GENES_NOT_FOUND.txt
            print("\nGenerating GENES_NOT_FOUND.txt...", flush=True)
            cmd = 'cut -f6 3_run_stats.tsv | tr "," "\\n" | grep -v "^GENES_NOT_FOUND$" | grep -v "^$" | python3 utils/normalize_genes.py | awk \'FNR==NR{a[$1]=1; next} !($1 in a)\' 3_gene_stats_found.tsv - > GENES_NOT_FOUND.txt'
            result = os.system(cmd)
            if result != 0 or not os.path.exists("GENES_NOT_FOUND.txt"):
                print("ERROR: Failed to generate GENES_NOT_FOUND.txt", flush=True)
                sys.exit(1)
            genes_not_found = "GENES_NOT_FOUND.txt"
            print("✓ Generated GENES_NOT_FOUND.txt", flush=True)
        else:
            parser.error("--genes-not-found is required for step 4, or provide --gene-stats and --run-stats to auto-generate it.")
            sys.exit(1)
        
        # Run step 4
        from step4.runner import run_step4
        print(f"\n{'='*80}")
        print("EXECUTING STEP 4")
        print(f"{'='*80}")
        print(f"""
        Arguments:
        \t--gene-stats:\t{in_gene_stats}
        \t--genes-not-found:\t{genes_not_found}
        \t--genedir:\t{geneseq_basedir}
        """)
        success = run_step4(in_gene_stats, genes_not_found, geneseq_basedir)
        if not success:
            print("\n❌ Step 4 failed. Stopping pipeline.", flush=True)
            sys.exit(1)
    
    # Final success message
    print("\n" + "="*80)
    if len(steps_to_run) == 1:
        print(f"✅ STEP {steps_to_run[0]} COMPLETED SUCCESSFULLY")
    elif args.step == 'all':
        print("✅ ALL STEPS COMPLETED SUCCESSFULLY")
    else:
        steps_str = ','.join(map(str, steps_to_run))
        print(f"✅ STEPS {steps_str} COMPLETED SUCCESSFULLY")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()
