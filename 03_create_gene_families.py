#!/usr/bin/env python3
"""
Step 3: Create gene families using MMseqs2 clustering
- Cluster all proteins at 80% identity
- Create presence/absence matrix
- Filter for families with at least 40 genomes
"""
import os
min_genomes_per_family = int(os.environ.get("MIN_GENOMES_PER_FAMILY", "40"))
print(f"\nUsing minimum {min_genomes_per_family} genomes per gene family")
import pandas as pd
import numpy as np
import subprocess
import os
from pathlib import Path
from collections import defaultdict
import sys

def concatenate_protein_sequences(protein_dir, output_file):
    """Concatenate all protein fasta files with genome IDs"""
    print("Concatenating protein sequences...")
    
    genome_protein_counts = {}
    
    with open(output_file, 'w') as outf:
        for faa_file in Path(protein_dir).glob('*.faa'):
            genome_id = faa_file.stem
            count = 0
            
            with open(faa_file, 'r') as inf:
                for line in inf:
                    if line.startswith('>'):
                        # Modify header to include genome ID
                        protein_id = line[1:].split()[0]
                        outf.write(f'>{genome_id}|{protein_id}\n')
                        count += 1
                    else:
                        outf.write(line)
            
            genome_protein_counts[genome_id] = count
    
    total_proteins = sum(genome_protein_counts.values())
    print(f"  Total proteins: {total_proteins:,}")
    print(f"  From {len(genome_protein_counts)} genomes")
    
    return genome_protein_counts

def run_mmseqs2_clustering(input_fasta, output_prefix, identity=0.8, coverage=0.8, threads=8):
    """Run MMseqs2 clustering"""
    print(f"\nRunning MMseqs2 clustering (identity={identity}, coverage={coverage})...")
    
    # Create MMseqs2 database
    print("  Creating database...")
    db_file = f'{output_prefix}_DB'
    cmd = ['mmseqs', 'createdb', input_fasta, db_file]
    subprocess.run(cmd, check=True, capture_output=True)
    
    # Cluster
    print("  Clustering (this may take 30-60 minutes)...")
    cluster_db = f'{output_prefix}_cluster'
    tmp_dir = f'{output_prefix}_tmp'
    os.makedirs(tmp_dir, exist_ok=True)
    
    cmd = [
        'mmseqs', 'cluster',
        db_file,
        cluster_db,
        tmp_dir,
        '--min-seq-id', str(identity),
        '-c', str(coverage),
        '--threads', str(threads),
        '--cluster-mode', '0',  # Set-cover clustering
        '--cov-mode', '0'  # Coverage of target
    ]
    subprocess.run(cmd, check=True)
    
    # Create TSV output
    print("  Creating TSV output...")
    tsv_file = f'{output_prefix}_clusters.tsv'
    cmd = [
        'mmseqs', 'createtsv',
        db_file,
        db_file,
        cluster_db,
        tsv_file,
        '--threads', str(threads)
    ]
    subprocess.run(cmd, check=True)
    
    print(f"  Clustering complete. Results in: {tsv_file}")
    
    # Clean up temporary files
    subprocess.run(['rm', '-rf', tmp_dir], capture_output=True)
    
    return tsv_file

def parse_clusters(tsv_file):
    """Parse MMseqs2 cluster file and create gene family matrix"""
    print("\nParsing clusters...")
    
    clusters = defaultdict(set)
    
    with open(tsv_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                representative = parts[0]
                member = parts[1]
                
                # Extract genome ID
                genome_id = member.split('|')[0]
                clusters[representative].add(genome_id)
    
    print(f"  Total gene families: {len(clusters)}")
    
    # Get cluster sizes
    cluster_sizes = {rep: len(genomes) for rep, genomes in clusters.items()}
    
    return clusters, cluster_sizes

def create_presence_absence_matrix(clusters, complete_genomes_df, min_genomes=40):
    """Create presence/absence matrix for gene families"""
    print(f"\nCreating presence/absence matrix (min {min_genomes} genomes per family)...")
    
    # Get list of all genomes
    all_genomes = set(complete_genomes_df['assembly_accession'])
    
    # Filter clusters by size
    filtered_clusters = {
        rep: genomes 
        for rep, genomes in clusters.items() 
        if len(genomes) >= min_genomes
    }
    
    print(f"  Gene families after filtering: {len(filtered_clusters)}")
    
    # Create matrix
    genome_list = sorted(all_genomes)
    cluster_list = sorted(filtered_clusters.keys())
    
    matrix = np.zeros((len(genome_list), len(cluster_list)), dtype=np.int8)
    
    for j, cluster_rep in enumerate(cluster_list):
        genomes_in_cluster = filtered_clusters[cluster_rep]
        for i, genome in enumerate(genome_list):
            if genome in genomes_in_cluster:
                matrix[i, j] = 1
    
    # Create DataFrame
    df_matrix = pd.DataFrame(
        matrix,
        index=genome_list,
        columns=[f'GF_{i:05d}' for i in range(len(cluster_list))]
    )
    
    print(f"  Matrix shape: {df_matrix.shape}")
    print(f"  Sparsity: {(1 - matrix.mean())*100:.1f}%")
    
    return df_matrix, filtered_clusters

def main():
    print("="*80)
    print("STEP 3: GENE FAMILY CLUSTERING")
    print("="*80)
    
    # Parameters
    protein_dir = 'downloads/proteins'
    identity = 0.8
    coverage = 0.8
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    
    # Check inputs
    if not os.path.exists(protein_dir):
        print(f"Error: {protein_dir} not found. Run 02_download_proteins.py first.")
        sys.exit(1)
    
    if not os.path.exists('complete_genomes_labeled.csv'):
        print("Error: complete_genomes_labeled.csv not found. Run 01_prepare_data.py first.")
        sys.exit(1)
    
    # Check how many protein files we have
    protein_files = list(Path(protein_dir).glob('*.faa'))
    print(f"\nFound {len(protein_files)} protein sequence files")
    
    if len(protein_files) == 0:
        print("\nERROR: No protein sequences found!")
        print("Download proteins first using 02_download_proteins.py")
        print("See MANUAL_DOWNLOAD.txt if automated download failed")
        sys.exit(1)
    
    # Load genome metadata
    print("\nLoading genome metadata...")
    complete_genomes = pd.read_csv('complete_genomes_labeled.csv')
    print(f"  Genomes in metadata: {len(complete_genomes)}")
    
    # Filter metadata to only genomes with proteins
    available_accessions = [f.stem for f in protein_files]
    complete_genomes = complete_genomes[
        complete_genomes['assembly_accession'].isin(available_accessions)
    ].copy()
    
    print(f"  Genomes with proteins: {len(complete_genomes)}")
    
    if len(complete_genomes) < 50:
        print(f"\nWARNING: Only {len(complete_genomes)} genomes available.")
        print("This may not be enough for robust analysis.")
        print("Consider downloading more proteins.")
    
    # Adjust min_genomes based on available data
    # Use 10% of genomes as minimum, but at least 10 and at most 40
    min_genomes_per_family = max(10, min(40, len(complete_genomes) // 10))
    print(f"\nUsing minimum {min_genomes_per_family} genomes per gene family")
    
    # Step 1: Concatenate all proteins
    all_proteins_file = 'all_proteins.faa'
    if not os.path.exists(all_proteins_file):
        genome_protein_counts = concatenate_protein_sequences(protein_dir, all_proteins_file)
        if sum(genome_protein_counts.values()) == 0:
            print("\nERROR: No proteins in concatenated file!")
            sys.exit(1)
    else:
        print(f"Using existing concatenated file: {all_proteins_file}")
    
    # Step 2: Run MMseqs2 clustering
    cluster_tsv = 'gene_families_clusters.tsv'
    if not os.path.exists(cluster_tsv):
        cluster_tsv = run_mmseqs2_clustering(
            all_proteins_file, 
            'gene_families',
            identity=identity,
            coverage=coverage,
            threads=threads
        )
    else:
        print(f"Using existing cluster file: {cluster_tsv}")
    
    # Step 3: Parse clusters
    clusters, cluster_sizes = parse_clusters(cluster_tsv)
    
    # Step 4: Create presence/absence matrix
    pa_matrix, filtered_clusters = create_presence_absence_matrix(
        clusters,
        complete_genomes,
        min_genomes=min_genomes_per_family
    )
    
    # Check if we have enough features
    if pa_matrix.shape[1] < 10:
        print(f"\nWARNING: Only {pa_matrix.shape[1]} gene families after filtering!")
        print("Consider:")
        print("  1. Lowering identity threshold (currently {identity})")
        print("  2. Reducing min_genomes requirement")
        print("  3. Downloading more genomes")
    
    # Save matrix
    output_file = 'gene_family_matrix.csv'
    pa_matrix.to_csv(output_file)
    print(f"\nSaved gene family matrix to: {output_file}")
    
    # Save cluster information
    cluster_info = pd.DataFrame([
        {'gene_family': f'GF_{i:05d}', 'representative': rep, 'num_genomes': len(genomes)}
        for i, (rep, genomes) in enumerate(sorted(filtered_clusters.items()))
    ])
    cluster_info.to_csv('gene_family_info.csv', index=False)
    print(f"Saved gene family info to: gene_family_info.csv")
    
    # Also save the filtered genome list
    complete_genomes.to_csv('complete_genomes_with_proteins.csv', index=False)
    print(f"Saved filtered genome metadata to: complete_genomes_with_proteins.csv")
    
    # Print summary
    print("\n" + "="*80)
    print("GENE FAMILY SUMMARY")
    print("="*80)
    print(f"Genomes analyzed: {len(complete_genomes)}")
    print(f"Total gene families (>={min_genomes_per_family} genomes): {len(filtered_clusters)}")
    print(f"Average family size: {cluster_info['num_genomes'].mean():.1f}")
    print(f"Median family size: {cluster_info['num_genomes'].median():.0f}")
    print(f"Max family size: {cluster_info['num_genomes'].max()}")
    
    print("\n" + "="*80)
    print("STEP 3 COMPLETE")
    print("="*80)

if __name__ == '__main__':
    main()
