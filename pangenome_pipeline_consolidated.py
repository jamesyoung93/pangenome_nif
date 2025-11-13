#!/usr/bin/env python3
"""
Consolidated Pangenome Analysis Pipeline
Combines steps 1-8 into a single executable script

Steps:
1. Prepare data and label genomes
2. Download protein sequences
3. Create gene families using MMseqs2
4. Classification with genus-level cross-validation
5. Analyze feature directionality
6. Merge Swiss-Prot annotations
7. Expand gene family information
8. Build narrative table

Usage:
    python pangenome_pipeline_consolidated.py --input INPUT_CSV [options]
"""

import argparse
import os
import sys
import time
import warnings
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
import subprocess
import gzip
import csv
import re
import shutil
import zipfile
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

# Data manipulation
import pandas as pd
import numpy as np
import ast

# Statistics
from scipy import stats

# Machine Learning
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, roc_curve
)
from sklearn.preprocessing import StandardScaler

# Plotting
import matplotlib
matplotlib.use('Agg')  # HPC-safe
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# Optional XGBoost
try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("Warning: XGBoost not available, will skip XGBoost models")

# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================

class PipelineConfig:
    """Configuration parameters for the pipeline"""
    def __init__(self):
        # Data paths
        self.input_csv = "nif_hdk_hits_enriched_with_quality_checkm.csv"
        self.protein_dir = "downloads/proteins"

        # Step 3: Gene family parameters
        self.mmseqs_identity = 0.8
        self.mmseqs_coverage = 0.8
        self.min_genomes_per_family = 40
        self.threads = 8

        # Step 4: Classification parameters
        self.purity_threshold = 0.9
        self.cv_folds = 5
        self.cv_seed = 42

        # Step 8: Top features
        self.top_n_features = 100

        # NCBI download parameters
        self.email = "james.young@jacks.sdstate.edu"
        self.download_jobs = 8
        self.download_sleep = 0.25
        self.download_retries = 3

config = PipelineConfig()

# ============================================================================
# STEP 1: PREPARE DATA AND LABEL GENOMES
# ============================================================================

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

def step1_prepare_data():
    """Step 1: Prepare data and label genomes as diazotroph/non-diazotroph"""
    print("\n" + "="*80)
    print("STEP 1: DATA PREPARATION AND LABELING")
    print("="*80)

    # Load the CSV
    print(f"\nLoading data from: {config.input_csv}")

    if not os.path.exists(config.input_csv):
        print(f"ERROR: Input file not found: {config.input_csv}")
        return False

    df = pd.read_csv(config.input_csv, on_bad_lines='skip', engine='python')
    print(f"Total genomes loaded: {len(df)}")

    # Parse organism names
    print("\nParsing organism names...")
    df['organism_full'] = df['organism_name'].apply(parse_organism_name)
    df['genus'] = df['organism_full'].apply(extract_genus)

    # Convert e-values to float
    print("Converting e-values...")
    for col in ['nifH_best_evalue', 'nifD_best_evalue', 'nifK_best_evalue']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Filter for ONLY Complete Genome
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

    print("\n" + "="*80)
    print("STEP 1 COMPLETE")
    print("="*80)

    return True

# ============================================================================
# STEP 2: DOWNLOAD PROTEINS (simplified for HPC)
# ============================================================================

def step2_download_proteins():
    """Step 2: Download protein sequences from NCBI"""
    print("\n" + "="*80)
    print("STEP 2: DOWNLOAD PROTEIN SEQUENCES")
    print("="*80)

    print("\nNOTE: This step downloads protein sequences from NCBI.")
    print("This may take several hours depending on the number of genomes.")
    print(f"Using {config.download_jobs} parallel downloads\n")

    # Check if accessions file exists
    if not os.path.exists('genome_accessions.txt'):
        print("ERROR: genome_accessions.txt not found. Run Step 1 first.")
        return False

    # Create output directory
    protein_dir = Path(config.protein_dir)
    protein_dir.mkdir(parents=True, exist_ok=True)

    # Check how many we already have
    existing = list(protein_dir.glob('*.faa'))
    print(f"Found {len(existing)} existing protein files")

    # For HPC environments, we'll use the existing 02_download_proteins.py
    # as it has robust retry logic and error handling
    if os.path.exists('02_download_proteins.py'):
        print("\nUsing existing download script with robust error handling...")
        cmd = [
            'python3', '02_download_proteins.py',
            '--csv', 'complete_genomes_labeled.csv',
            '--complete-only',
            '--outdir', 'downloads',
            '--jobs', str(config.download_jobs),
            '--sleep', str(config.download_sleep),
            '--retries', str(config.download_retries),
            '--email', config.email
        ]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print("\nWARNING: Download encountered errors. Check download_log.tsv")
            print("Pipeline will continue with available proteins.")
    else:
        print("\nWARNING: 02_download_proteins.py not found.")
        print("Please ensure protein sequences are available in:", protein_dir)

    # Check results
    existing_after = list(protein_dir.glob('*.faa'))
    print(f"\nProtein files available: {len(existing_after)}")

    if len(existing_after) < 50:
        print("\nWARNING: Very few protein files available.")
        print("Analysis may not be robust with limited data.")
        response = input("Continue anyway? (yes/no): ")
        if response.lower() != 'yes':
            return False

    print("\n" + "="*80)
    print("STEP 2 COMPLETE")
    print("="*80)

    return True

# ============================================================================
# STEP 3: CREATE GENE FAMILIES
# ============================================================================

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
        '--cluster-mode', '0',
        '--cov-mode', '0'
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
    """Parse MMseqs2 cluster file"""
    print("\nParsing clusters...")

    clusters = defaultdict(set)

    with open(tsv_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                representative = parts[0]
                member = parts[1]
                genome_id = member.split('|')[0]
                clusters[representative].add(genome_id)

    print(f"  Total gene families: {len(clusters)}")

    cluster_sizes = {rep: len(genomes) for rep, genomes in clusters.items()}

    return clusters, cluster_sizes

def create_presence_absence_matrix(clusters, complete_genomes_df, min_genomes=40):
    """Create presence/absence matrix for gene families"""
    print(f"\nCreating presence/absence matrix (min {min_genomes} genomes per family)...")

    all_genomes = set(complete_genomes_df['assembly_accession'])

    filtered_clusters = {
        rep: genomes
        for rep, genomes in clusters.items()
        if len(genomes) >= min_genomes
    }

    print(f"  Gene families after filtering: {len(filtered_clusters)}")

    genome_list = sorted(all_genomes)
    cluster_list = sorted(filtered_clusters.keys())

    matrix = np.zeros((len(genome_list), len(cluster_list)), dtype=np.int8)

    for j, cluster_rep in enumerate(cluster_list):
        genomes_in_cluster = filtered_clusters[cluster_rep]
        for i, genome in enumerate(genome_list):
            if genome in genomes_in_cluster:
                matrix[i, j] = 1

    df_matrix = pd.DataFrame(
        matrix,
        index=genome_list,
        columns=[f'GF_{i:05d}' for i in range(len(cluster_list))]
    )

    print(f"  Matrix shape: {df_matrix.shape}")
    print(f"  Sparsity: {(1 - matrix.mean())*100:.1f}%")

    return df_matrix, filtered_clusters

def step3_create_gene_families():
    """Step 3: Create gene families using MMseqs2"""
    print("\n" + "="*80)
    print("STEP 3: GENE FAMILY CLUSTERING")
    print("="*80)

    # Check inputs
    if not os.path.exists(config.protein_dir):
        print(f"ERROR: {config.protein_dir} not found. Run Step 2 first.")
        return False

    if not os.path.exists('complete_genomes_labeled.csv'):
        print("ERROR: complete_genomes_labeled.csv not found. Run Step 1 first.")
        return False

    protein_files = list(Path(config.protein_dir).glob('*.faa'))
    print(f"\nFound {len(protein_files)} protein sequence files")

    if len(protein_files) == 0:
        print("\nERROR: No protein sequences found!")
        return False

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

    # Adjust min_genomes based on available data
    min_genomes_per_family = max(10, min(config.min_genomes_per_family, len(complete_genomes) // 10))
    print(f"\nUsing minimum {min_genomes_per_family} genomes per gene family")

    # Step 1: Concatenate all proteins
    all_proteins_file = 'all_proteins.faa'
    if not os.path.exists(all_proteins_file):
        genome_protein_counts = concatenate_protein_sequences(config.protein_dir, all_proteins_file)
        if sum(genome_protein_counts.values()) == 0:
            print("\nERROR: No proteins in concatenated file!")
            return False
    else:
        print(f"Using existing concatenated file: {all_proteins_file}")

    # Step 2: Run MMseqs2 clustering
    cluster_tsv = 'gene_families_clusters.tsv'
    if not os.path.exists(cluster_tsv):
        cluster_tsv = run_mmseqs2_clustering(
            all_proteins_file,
            'gene_families',
            identity=config.mmseqs_identity,
            coverage=config.mmseqs_coverage,
            threads=config.threads
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

    # Save filtered genome list
    complete_genomes.to_csv('complete_genomes_with_proteins.csv', index=False)
    print(f"Saved filtered genome metadata to: complete_genomes_with_proteins.csv")

    print("\n" + "="*80)
    print("STEP 3 COMPLETE")
    print("="*80)

    return True

# ============================================================================
# STEP 4: CLASSIFICATION
# ============================================================================

def _safe_select(df_or_ser, idx):
    """Select rows by either integer positions or index labels"""
    if len(idx) == 0:
        return df_or_ser.iloc[[]]
    first = idx[0]
    if isinstance(first, (int, np.integer)):
        return df_or_ser.iloc[idx]
    else:
        return df_or_ser.loc[idx]

def filter_pure_gene_families(X, y, threshold=0.9):
    """Remove gene families that are highly pure"""
    print(f"\nFiltering gene families (removing >={threshold*100:.0f}% or <={100-threshold*100:.0f}% diazotroph)...")

    family_stats = []
    for col in X.columns:
        present = X[col] == 1
        total_count = int(present.sum())
        if total_count == 0:
            continue
        diaz_count = int(y[present].sum())
        diaz_pct = diaz_count / total_count
        family_stats.append({
            'gene_family': col,
            'total_genomes': total_count,
            'diazotroph_count': diaz_count,
            'diazotroph_pct': diaz_pct,
            'is_pure': (diaz_pct >= threshold) or (diaz_pct <= (1 - threshold))
        })

    stats_df = pd.DataFrame(family_stats).sort_values('gene_family')
    pure_families = stats_df.loc[stats_df['is_pure'], 'gene_family'].tolist()
    informative_families = stats_df.loc[~stats_df['is_pure'], 'gene_family'].tolist()

    print(f"  Total gene families: {len(stats_df)}")
    print(f"  Pure families (removed): {len(pure_families)}")
    print(f"  Informative families (kept): {len(informative_families)}")

    stats_df.to_csv('gene_family_purity_stats.csv', index=False)
    print("  Saved purity statistics to: gene_family_purity_stats.csv")

    X_filtered = X[informative_families]

    return X_filtered, stats_df

def genus_cross_validation_split(metadata_df, n_folds=5, seed=42):
    """Create genus-level cross-validation folds"""
    rng = np.random.default_rng(seed)

    meta = metadata_df.copy()
    unique_genera = meta['genus'].dropna().unique().tolist()
    rng.shuffle(unique_genera)

    folds_genera = [[] for _ in range(n_folds)]
    for i, g in enumerate(unique_genera):
        folds_genera[i % n_folds].append(g)

    folds = []
    for fold_genus in folds_genera:
        idx_labels = meta.index[meta['genus'].isin(fold_genus)].unique().tolist()
        folds.append(idx_labels)

    return folds

def evaluate_model(y_true, y_pred, y_prob):
    """Calculate evaluation metrics"""
    metrics = {
        'accuracy': accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, zero_division=0),
        'recall': recall_score(y_true, y_pred, zero_division=0),
        'f1': f1_score(y_true, y_pred, zero_division=0),
        'roc_auc': roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else 0.0
    }
    return metrics

def run_genus_cv(X, y, metadata_df, model, model_name, n_folds=5, seed=42):
    """Run genus-level cross-validation"""
    print(f"\n{'='*60}")
    print(f"Running {model_name} with Genus-level CV")
    print(f"{'='*60}")

    folds = genus_cross_validation_split(metadata_df, n_folds=n_folds, seed=seed)

    fold_metrics = []
    all_y_true = []
    all_y_prob = []
    feature_importances = []

    all_labels = X.index.tolist()

    for fold_idx, test_labels in enumerate(folds, start=1):
        test_labels = [lab for lab in test_labels if lab in X.index]
        test_set = set(test_labels)
        train_labels = [lab for lab in all_labels if lab not in test_set]

        X_train = _safe_select(X, train_labels)
        X_test  = _safe_select(X, test_labels)
        y_train = y.loc[X_train.index]
        y_test  = y.loc[X_test.index]

        print(f"\nFold {fold_idx}/{n_folds}")
        print(f"  Train: {len(X_train)} genomes, Test: {len(X_test)} genomes")

        if model_name in ['Logistic Regression', 'SVM']:
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
        else:
            X_train_scaled = X_train
            X_test_scaled = X_test

        model.fit(X_train_scaled, y_train)
        y_pred = model.predict(X_test_scaled)

        if hasattr(model, 'predict_proba'):
            y_prob = model.predict_proba(X_test_scaled)[:, 1]
        else:
            raw = model.decision_function(X_test_scaled)
            y_prob = (raw - raw.min()) / (raw.max() - raw.min() + 1e-12)

        metrics = evaluate_model(y_test, y_pred, y_prob)
        fold_metrics.append(metrics)

        print(f"  Accuracy: {metrics['accuracy']:.4f}")
        print(f"  F1: {metrics['f1']:.4f}")
        print(f"  ROC-AUC: {metrics['roc_auc']:.4f}")

        all_y_true.extend(y_test.tolist())
        all_y_prob.extend(y_prob.tolist())

        if hasattr(model, 'feature_importances_'):
            feature_importances.append(model.feature_importances_)
        elif hasattr(model, 'coef_'):
            feature_importances.append(np.abs(model.coef_[0]))

    metrics_df = pd.DataFrame(fold_metrics)
    mean_metrics = metrics_df.mean()
    std_metrics = metrics_df.std()

    print(f"\n{'='*60}")
    print(f"{model_name} - AGGREGATE RESULTS")
    print(f"{'='*60}")
    for metric in ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']:
        print(f"{metric.upper():12s}: {mean_metrics[metric]:.4f} +/- {std_metrics[metric]:.4f}")

    if feature_importances:
        avg_importance = np.mean(feature_importances, axis=0)
        feature_importance_df = pd.DataFrame({
            'gene_family': X.columns,
            'importance': avg_importance
        }).sort_values('importance', ascending=False)
    else:
        feature_importance_df = None

    return {
        'fold_metrics': metrics_df,
        'mean_metrics': mean_metrics,
        'std_metrics': std_metrics,
        'all_y_true': np.array(all_y_true),
        'all_y_prob': np.array(all_y_prob),
        'feature_importance': feature_importance_df
    }

def plot_roc_curves(results_dict, output_file='roc_curves.png'):
    """Plot ROC curves for all models"""
    plt.figure(figsize=(10, 8))

    for model_name, results in results_dict.items():
        fpr, tpr, _ = roc_curve(results['all_y_true'], results['all_y_prob'])
        auc = results['mean_metrics']['roc_auc']
        auc_std = results['std_metrics']['roc_auc']
        plt.plot(fpr, tpr, label=f"{model_name} (AUC = {auc:.3f} +/- {auc_std:.3f})", linewidth=2)

    plt.plot([0, 1], [0, 1], 'k--', label='Random Classifier', linewidth=1)
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC Curves - Genus-Level Cross-Validation', fontsize=14, fontweight='bold')
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nSaved ROC curves to: {output_file}")

def analyze_feature_importance(results_dict, top_n=20):
    """Analyze and save feature importance across models"""
    all_importances = []

    for model_name, results in results_dict.items():
        fi = results['feature_importance']
        if fi is not None:
            tmp = fi.copy()
            tmp['model'] = model_name
            all_importances.append(tmp)

    if not all_importances:
        print("No feature importances available")
        return

    combined = pd.concat(all_importances, ignore_index=True)

    avg_importance = (combined.groupby('gene_family')['importance']
                      .mean().sort_values(ascending=False))
    top_features = avg_importance.head(top_n).index

    pivot = (combined[combined['gene_family'].isin(top_features)]
             .pivot(index='gene_family', columns='model', values='importance')
             .fillna(0.0)
             .loc[top_features])

    pivot.to_csv('feature_importance_top.csv')
    print(f"\nSaved top {top_n} feature importances to: feature_importance_top.csv")

    fig, ax = plt.subplots(figsize=(12, 8))
    pivot.plot(kind='barh', ax=ax)
    ax.set_xlabel('Feature Importance', fontsize=12)
    ax.set_ylabel('Gene Family', fontsize=12)
    ax.set_title(f'Top {top_n} Most Important Gene Families', fontsize=14, fontweight='bold')
    ax.legend(title='Model', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('feature_importance_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved feature importance plot to: feature_importance_plot.png")

def step4_classification():
    """Step 4: Classification with genus-level cross-validation"""
    print("\n" + "="*80)
    print("STEP 4: CLASSIFICATION WITH GENUS-LEVEL CROSS-VALIDATION")
    print("="*80)

    # Load data
    print("\nLoading data...")
    gene_family_matrix = pd.read_csv('gene_family_matrix.csv', index_col=0)

    if os.path.exists('complete_genomes_with_proteins.csv'):
        metadata = pd.read_csv('complete_genomes_with_proteins.csv')
    else:
        metadata = pd.read_csv('complete_genomes_labeled.csv')

    # Deduplicate
    if metadata.duplicated('assembly_accession').any():
        dup_count = int(metadata.duplicated('assembly_accession').sum())
        print(f"  De-duplicating metadata (dropping {dup_count} duplicates)")
        metadata = metadata.drop_duplicates('assembly_accession', keep='first')

    # Align data
    common = sorted(set(metadata['assembly_accession']) & set(gene_family_matrix.index))
    print(f"  Common genomes: {len(common)}")

    metadata = (metadata[metadata['assembly_accession'].isin(common)]
                .copy()
                .set_index('assembly_accession')
                .loc[common])
    metadata = metadata[~metadata.index.duplicated(keep='first')]
    X_full = gene_family_matrix.loc[common]
    X_full = X_full[~X_full.index.duplicated(keep='first')]

    final_common = X_full.index.intersection(metadata.index)
    X_full = X_full.loc[final_common]
    metadata = metadata.loc[final_common]

    print(f"  Gene families: {X_full.shape[1]}")
    print(f"  Genomes: {X_full.shape[0]}")
    print(f"  Diazotrophs: {int(metadata['is_diazotroph'].sum())} ({metadata['is_diazotroph'].mean()*100:.1f}%)")

    # Filter pure families
    X, purity_stats = filter_pure_gene_families(
        X_full,
        metadata['is_diazotroph'],
        threshold=config.purity_threshold
    )

    y = metadata['is_diazotroph'].astype(int)

    # Filter to families in multiple genera
    def filter_min_genera(X, meta, k=3):
        keep = []
        for gf in X.columns:
            g = meta.loc[X[gf] == 1, "genus"].dropna().nunique()
            if g >= k:
                keep.append(gf)
        return X[keep]

    X = filter_min_genera(X, metadata, k=3)
    print(f"\nFinal feature matrix: {X.shape}")

    # Determine folds
    n_genera = metadata['genus'].nunique()
    n_folds = min(config.cv_folds, max(3, n_genera // 3))
    print(f"\nUsing {n_folds}-fold cross-validation ({n_genera} genera)")

    # Define models
    models = {
        'Random Forest': RandomForestClassifier(
            n_estimators=100, max_depth=10, min_samples_leaf=5,
            random_state=42, n_jobs=-1
        ),
        'Gradient Boosting': GradientBoostingClassifier(
            n_estimators=100, max_depth=5, learning_rate=0.1,
            random_state=42
        ),
        'Logistic Regression': LogisticRegression(
            max_iter=1000, random_state=42
        ),
    }
    if HAS_XGB:
        models['XGBoost'] = xgb.XGBClassifier(
            n_estimators=200, max_depth=5, learning_rate=0.1,
            subsample=0.8, colsample_bytree=0.8,
            random_state=42, n_jobs=-1, use_label_encoder=False,
            eval_metric='logloss'
        )

    # Run CV
    results_dict = {}
    for model_name, model in models.items():
        results = run_genus_cv(X, y, metadata, model, model_name, n_folds=n_folds, seed=config.cv_seed)
        results_dict[model_name] = results

    # Summary
    print("\n" + "="*80)
    print("SUMMARY OF ALL MODELS")
    print("="*80)

    summary_rows = []
    for model_name, results in results_dict.items():
        row = {'Model': model_name}
        for metric in ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']:
            mean = results['mean_metrics'][metric]
            std = results['std_metrics'][metric]
            row[metric] = f"{mean:.4f} +/- {std:.4f}"
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))
    summary_df.to_csv('classification_summary.csv', index=False)
    print("\nSaved summary to: classification_summary.csv")

    # Plots
    plot_roc_curves(results_dict)
    analyze_feature_importance(results_dict, top_n=30)

    # Save per-model results
    for model_name, results in results_dict.items():
        fold_df = results['fold_metrics'].copy()
        fold_df['model'] = model_name
        fold_df.to_csv(f"fold_metrics_{model_name.replace(' ', '_')}.csv", index=False)

        if results['feature_importance'] is not None:
            results['feature_importance'].to_csv(
                f"feature_importance_{model_name.replace(' ', '_')}.csv", index=False
            )

    print("\n" + "="*80)
    print("STEP 4 COMPLETE")
    print("="*80)

    return True

# ============================================================================
# STEP 5: FEATURE DIRECTIONALITY
# ============================================================================

def calculate_feature_directionality(X, y, feature_importance_df):
    """Calculate feature directionality"""
    print("\nCalculating feature directionality...")
    results = []

    y = y.astype(int).loc[X.index]

    for feature in X.columns:
        present = (X[feature] == 1)
        n_present = int(present.sum())
        n_absent  = int((~present).sum())

        diaz_rate_present = y[present].mean() if n_present > 0 else 0.0
        diaz_rate_absent  = y[~present].mean() if n_absent  > 0 else 0.0

        effect_size = float(diaz_rate_present - diaz_rate_absent)
        direction = 'Positive' if effect_size > 0 else ('Negative' if effect_size < 0 else 'Neutral')

        try:
            table = pd.crosstab(present, y)
            chi2, p_value = (0.0, 1.0)
            if table.shape == (2, 2):
                chi2, p_value, _, _ = stats.chi2_contingency(table)
        except Exception:
            chi2, p_value = 0.0, 1.0

        results.append({
            'gene_family': feature,
            'direction': direction,
            'effect_size': effect_size,
            'diaz_rate_present': float(diaz_rate_present),
            'diaz_rate_absent': float(diaz_rate_absent),
            'n_genomes_with_feature': n_present,
            'chi2': float(chi2),
            'p_value': float(p_value)
        })

    direction_df = pd.DataFrame(results)

    if feature_importance_df is not None and len(feature_importance_df):
        avg_imp = (feature_importance_df
                   .groupby('gene_family', as_index=False)['importance']
                   .mean())
        direction_df = direction_df.merge(avg_imp, on='gene_family', how='left')
        direction_df['importance'] = direction_df['importance'].fillna(0.0)
    else:
        direction_df['importance'] = 0.0

    direction_df = direction_df.sort_values(['importance', 'effect_size'], ascending=[False, False])
    return direction_df

def plot_feature_effects(direction_df, top_n=30, output_file='feature_effects.png'):
    """Plot feature effects"""
    top_n = min(top_n, len(direction_df))
    if top_n == 0:
        return

    top = direction_df.head(top_n).copy()
    top['signed_importance'] = top['importance'] * np.sign(top['effect_size'])
    top = top.sort_values('signed_importance')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))

    colors1 = ['green' if v > 0 else ('red' if v < 0 else 'gray') for v in top['signed_importance']]
    ax1.barh(range(top_n), top['signed_importance'], color=colors1, alpha=0.8)
    ax1.set_yticks(range(top_n))
    ax1.set_yticklabels(top['gene_family'], fontsize=8)
    ax1.set_xlabel('Signed Importance', fontsize=12)
    ax1.set_title('Feature Importance with Directionality', fontsize=12)
    ax1.axvline(0, color='black', linestyle='--', linewidth=1)
    ax1.grid(axis='x', alpha=0.3)

    colors2 = ['green' if v > 0 else ('red' if v < 0 else 'gray') for v in top['effect_size']]
    ax2.barh(range(top_n), top['effect_size'], color=colors2, alpha=0.8)
    ax2.set_yticks(range(top_n))
    ax2.set_yticklabels(top['gene_family'], fontsize=8)
    ax2.set_xlabel('Effect Size', fontsize=12)
    ax2.set_title('Feature Effect Sizes', fontsize=12)
    ax2.axvline(0, color='black', linestyle='--', linewidth=1)
    ax2.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved feature effects plot to: {output_file}")

def step5_feature_directionality():
    """Step 5: Analyze feature directionality"""
    print("\n" + "="*80)
    print("STEP 5: FEATURE DIRECTIONALITY ANALYSIS")
    print("="*80)

    print("\nLoading data...")
    X_full = pd.read_csv('gene_family_matrix.csv', index_col=0)

    if os.path.exists('complete_genomes_with_proteins.csv'):
        meta = pd.read_csv('complete_genomes_with_proteins.csv')
    else:
        meta = pd.read_csv('complete_genomes_labeled.csv')

    if meta.duplicated('assembly_accession').any():
        meta = meta.drop_duplicates('assembly_accession', keep='first')

    common = sorted(set(meta['assembly_accession']) & set(X_full.index))
    meta = (meta[meta['assembly_accession'].isin(common)]
            .set_index('assembly_accession')
            .loc[common])
    X_full = X_full.loc[common]

    purity = pd.read_csv('gene_family_purity_stats.csv')
    informative = purity.loc[~purity['is_pure'], 'gene_family'].tolist()
    X = X_full[informative]
    y = meta['is_diazotroph'].astype(int)

    print(f"  Analyzing {X.shape[1]} gene families across {X.shape[0]} genomes")

    # Load importances
    imports = []
    for model_tag in ['Random_Forest', 'Gradient_Boosting', 'Logistic_Regression', 'XGBoost']:
        fn = f'feature_importance_{model_tag}.csv'
        if os.path.exists(fn):
            df = pd.read_csv(fn)
            if {'gene_family', 'importance'}.issubset(df.columns):
                df['model'] = model_tag
                imports.append(df)
    feature_importance_df = pd.concat(imports, ignore_index=True) if imports else None

    direction_df = calculate_feature_directionality(X, y, feature_importance_df)
    direction_df.to_csv('feature_directionality_full.csv', index=False)
    print("Saved full directionality results to: feature_directionality_full.csv")

    plot_feature_effects(direction_df, top_n=30)

    # Summary
    top = direction_df.head(50)
    top.to_csv('feature_directionality_summary.csv', index=False)
    print("Saved summary to: feature_directionality_summary.csv")

    print("\n" + "="*80)
    print("STEP 5 COMPLETE")
    print("="*80)

    return True

# ============================================================================
# STEP 6: MERGE ANNOTATIONS (simplified - skip if files not available)
# ============================================================================

def step6_merge_annotations():
    """Step 6: Merge Swiss-Prot annotations (optional)"""
    print("\n" + "="*80)
    print("STEP 6: MERGE ANNOTATIONS (OPTIONAL)")
    print("="*80)

    # Check if Swiss-Prot results exist
    if not os.path.exists('gf_reps_vs_sprot.tsv'):
        print("\nSwiss-Prot results not found - skipping annotation merge")
        print("This step is optional and can be run separately if needed")
        print("\n" + "="*80)
        print("STEP 6 SKIPPED")
        print("="*80)
        return True

    # If 06 script exists, run it
    if os.path.exists('06_merge_annotations_selected.py'):
        print("\nRunning annotation merge script...")
        result = subprocess.run(['python3', '06_merge_annotations_selected.py'])
        if result.returncode == 0:
            print("\n" + "="*80)
            print("STEP 6 COMPLETE")
            print("="*80)
            return True

    print("\nAnnotation merge script not found - skipping")
    print("\n" + "="*80)
    print("STEP 6 SKIPPED")
    print("="*80)
    return True

# ============================================================================
# STEP 7: EXPAND GENE FAMILIES (simplified)
# ============================================================================

def step7_expand_gene_families():
    """Step 7: Expand gene family information"""
    print("\n" + "="*80)
    print("STEP 7: EXPAND GENE FAMILY INFORMATION")
    print("="*80)

    if os.path.exists('07_expand_gene_families.py'):
        print("\nRunning gene family expansion...")
        result = subprocess.run(['python3', '07_expand_gene_families.py'])
        if result.returncode == 0:
            print("\n" + "="*80)
            print("STEP 7 COMPLETE")
            print("="*80)
            return True
        else:
            print("\nWARNING: Step 7 encountered errors but continuing...")
            return True
    else:
        print("\nGene family expansion script not found - skipping")
        print("\n" + "="*80)
        print("STEP 7 SKIPPED")
        print("="*80)
        return True

# ============================================================================
# STEP 8: BUILD NARRATIVE TABLE (simplified)
# ============================================================================

def step8_build_narrative_table():
    """Step 8: Build narrative table"""
    print("\n" + "="*80)
    print("STEP 8: BUILD NARRATIVE TABLE")
    print("="*80)

    if os.path.exists('08_build_narrative_table.py'):
        print(f"\nBuilding narrative table for top {config.top_n_features} features...")
        env = os.environ.copy()
        env['TOP_N'] = str(config.top_n_features)
        result = subprocess.run(['python3', '08_build_narrative_table.py'], env=env)
        if result.returncode == 0:
            print("\n" + "="*80)
            print("STEP 8 COMPLETE")
            print("="*80)
            return True
        else:
            print("\nWARNING: Step 8 encountered errors but continuing...")
            return True
    else:
        print("\nNarrative table script not found - skipping")
        print("\n" + "="*80)
        print("STEP 8 SKIPPED")
        print("="*80)
        return True

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Consolidated Pangenome Analysis Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', type=str,
                        default='nif_hdk_hits_enriched_with_quality_checkm.csv',
                        help='Input CSV file with genome data')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads for MMseqs2')
    parser.add_argument('--min-genomes', type=int, default=40,
                        help='Minimum genomes per gene family')
    parser.add_argument('--start-step', type=int, default=1,
                        help='Start from this step (1-8)')
    parser.add_argument('--end-step', type=int, default=8,
                        help='End at this step (1-8)')
    parser.add_argument('--skip-download', action='store_true',
                        help='Skip protein download (assume already done)')

    args = parser.parse_args()

    # Update config
    config.input_csv = args.input
    config.threads = args.threads
    config.min_genomes_per_family = args.min_genomes

    print("\n" + "="*80)
    print("CONSOLIDATED PANGENOME ANALYSIS PIPELINE")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Input: {config.input_csv}")
    print(f"  Threads: {config.threads}")
    print(f"  Min genomes per family: {config.min_genomes_per_family}")
    print(f"  Steps: {args.start_step} to {args.end_step}")
    print("="*80)

    start_time = time.time()

    # Define pipeline steps
    steps = [
        (1, "Prepare Data", step1_prepare_data),
        (2, "Download Proteins", step2_download_proteins),
        (3, "Create Gene Families", step3_create_gene_families),
        (4, "Classification", step4_classification),
        (5, "Feature Directionality", step5_feature_directionality),
        (6, "Merge Annotations", step6_merge_annotations),
        (7, "Expand Gene Families", step7_expand_gene_families),
        (8, "Build Narrative Table", step8_build_narrative_table),
    ]

    # Run selected steps
    for step_num, step_name, step_func in steps:
        if step_num < args.start_step or step_num > args.end_step:
            continue

        if args.skip_download and step_num == 2:
            print(f"\nSkipping Step {step_num}: {step_name}")
            continue

        try:
            success = step_func()
            if not success:
                print(f"\nERROR: Step {step_num} failed. Stopping pipeline.")
                sys.exit(1)
        except Exception as e:
            print(f"\nERROR in Step {step_num}: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    elapsed = time.time() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)

    print("\n" + "="*80)
    print("PIPELINE COMPLETE!")
    print("="*80)
    print(f"\nTotal time: {hours}h {minutes}m")
    print("\nKey output files:")
    print("  - complete_genomes_labeled.csv")
    print("  - gene_family_matrix.csv")
    print("  - classification_summary.csv")
    print("  - roc_curves.png")
    print("  - feature_importance_*.csv")
    print("  - feature_directionality_full.csv")
    print("  - feature_effects.png")
    if os.path.exists('narrative_top_features.csv'):
        print("  - narrative_top_features.csv")
    print("="*80)

if __name__ == '__main__':
    main()
