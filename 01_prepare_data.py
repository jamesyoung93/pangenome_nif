#!/usr/bin/env python3
"""
Step 1: Prepare data and label genomes as diazotroph/non-diazotroph
Filter for complete genomes and download protein sequences
"""

import pandas as pd
import numpy as np
import ast
import os
import sys
from pathlib import Path
import subprocess
import time

def parse_organism_name(name_str):
    """Parse organism name from dictionary string"""
    try:
        parsed = ast.literal_eval(name_str) if pd.notna(name_str) and name_str.startswith('{') else {}
        return parsed.get('organism_name', '')
    except:
        return ''

def extract_genus(organism_name):
    """Extract genus from organism name"""
    return organism_name.split()[0] if organism_name else ''

def main():
    print("="*80)
    print("STEP 1: DATA PREPARATION AND LABELING")
    print("="*80)
    
    # Load the CSV
    csv_path = sys.argv[1] if len(sys.argv) > 1 else 'nif_hdk_hits_enriched_with_quality_checkm.csv'
    print(f"\nLoading data from: {csv_path}")
    
    df = pd.read_csv(csv_path, on_bad_lines='skip', engine='python')
    print(f"Total genomes loaded: {len(df)}")
    
    # Parse organism names
    print("\nParsing organism names...")
    df['organism_full'] = df['organism_name'].apply(parse_organism_name)
    df['genus'] = df['organism_full'].apply(extract_genus)
    
    # Convert e-values to float
    print("Converting e-values...")
    for col in ['nifH_best_evalue', 'nifD_best_evalue', 'nifK_best_evalue']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Filter for ONLY Complete Genome (NOT Chromosome)
    print("\nFiltering for Complete Genome ONLY...")
    complete = df[df['assembly_level'] == 'Complete Genome'].copy()
    print(f"Complete genomes: {len(complete)}")
    
    # Label as diazotroph if ALL nifH, nifD, nifK are below 1e-50
    print("\nLabeling diazotrophs (all nifH, D, K < 1e-50)...")
    complete['is_diazotroph'] = (
        (complete['nifH_best_evalue'] < 1e-50) & 
        (complete['nifD_best_evalue'] < 1e-50) & 
        (complete['nifK_best_evalue'] < 1e-50)
    )
    
    print(f"  Diazotrophs: {complete['is_diazotroph'].sum()}")
    print(f"  Non-diazotrophs: {(~complete['is_diazotroph']).sum()}")
    print(f"  Unique genera: {complete['genus'].nunique()}")
    
    # Save labeled dataset
    output_file = 'complete_genomes_labeled.csv'
    complete.to_csv(output_file, index=False)
    print(f"\nSaved labeled dataset to: {output_file}")
    
    # Create list of accessions for download
    accessions_file = 'genome_accessions.txt'
    complete['assembly_accession'].to_csv(accessions_file, index=False, header=False)
    print(f"Saved accession list to: {accessions_file}")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"\nTotal complete genomes: {len(complete)}")
    print(f"Diazotrophs: {complete['is_diazotroph'].sum()} ({complete['is_diazotroph'].sum()/len(complete)*100:.1f}%)")
    print(f"Non-diazotrophs: {(~complete['is_diazotroph']).sum()} ({(~complete['is_diazotroph']).sum()/len(complete)*100:.1f}%)")
    
    print("\nTop 10 genera:")
    print(complete['genus'].value_counts().head(10))
    
    print("\nDiazotroph distribution by top genera (n >= 5):")
    genus_stats = complete.groupby('genus').agg({
        'is_diazotroph': ['sum', 'count']
    })
    genus_stats.columns = ['diazotroph_count', 'total_count']
    genus_stats['percent_diazotroph'] = (genus_stats['diazotroph_count'] / genus_stats['total_count'] * 100).round(1)
    genus_stats = genus_stats[genus_stats['total_count'] >= 5].sort_values('total_count', ascending=False)
    print(genus_stats.head(15))
    
    print("\n" + "="*80)
    print("STEP 1 COMPLETE")
    print("="*80)

if __name__ == '__main__':
    main()
