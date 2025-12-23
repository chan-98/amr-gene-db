#!/usr/bin/env python3
import sys
import os
import subprocess
from pathlib import Path

def run_step1(input_file, refseq_basedir, geneseq_basedir):
    """Run Step 1: CDS extraction using bash script"""
    print("="*80)
    print("STEP 1: Extracting genes from CDS files")
    print("="*80)
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up one level to the package root
    package_root = os.path.dirname(script_dir)    
    # Look for seq_retriever.sh in the package root
    bash_script = os.path.join(package_root, "seq_retriever.sh")

    if not os.path.exists(bash_script):
        print(f"ERROR: {bash_script} not found in current directory")
        return False
    
    # Make script executable
    os.chmod(bash_script, 0o755)
    
    # Run the bash script
    cmd = [bash_script, input_file, refseq_basedir, geneseq_basedir]
    print(f"Running: {' '.join(cmd)}", flush=True)
    
    try:
        result = subprocess.run(cmd, check=True)
        print("\nStep 1 completed successfully ✅", flush=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"\nStep 1 failed with error: {e}", flush=True)
        return False
