#!/usr/bin/env python3
"""
Step 5: Analyze feature directionality
Determine if features are positively or negatively associated with diazotrophy
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use('Agg')  # HPC-safe
import matplotlib.pyplot as plt

# ---------------------- helpers ----------------------

def calculate_feature_directionality(X, y, feature_importance_df):
    """
    Positive = higher diazotroph rate when feature is present.
    Negative = higher non-diazotroph rate when feature is present.
    """
    print("\nCalculating feature directionality...")
    results = []

    # Ensure y is 0/1 Series aligned to X
    y = y.astype(int).loc[X.index]

    for feature in X.columns:
        present = (X[feature] == 1)
        n_present = int(present.sum())
        n_absent  = int((~present).sum())

        diaz_rate_present = y[present].mean() if n_present > 0 else 0.0
        diaz_rate_absent  = y[~present].mean() if n_absent  > 0 else 0.0

        effect_size = float(diaz_rate_present - diaz_rate_absent)
        direction = 'Positive' if effect_size > 0 else ('Negative' if effect_size < 0 else 'Neutral')

        # 2x2 contingency and chi-square
        # rows: present/absent; cols: diaz=1/non=0
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

    # Merge with feature importance (average across models)
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
    """Plot signed importance and raw effect sizes for top-N features."""
    top_n = min(top_n, len(direction_df))
    if top_n == 0:
        print("No features to plot.")
        return

    top = direction_df.head(top_n).copy()
    # signed importance: importance * sign(effect)
    top['signed_importance'] = top['importance'] * np.sign(top['effect_size'])
    top = top.sort_values('signed_importance')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))

    # Signed importance
    colors1 = ['green' if v > 0 else ('red' if v < 0 else 'gray') for v in top['signed_importance']]
    ax1.barh(range(top_n), top['signed_importance'], color=colors1, alpha=0.8)
    ax1.set_yticks(range(top_n))
    ax1.set_yticklabels(top['gene_family'], fontsize=8)
    ax1.set_xlabel('Signed Importance', fontsize=12)
    ax1.set_title('Feature Importance with Directionality\n(Green=Diazotroph, Red=Non-diazotroph)', fontsize=12)
    ax1.axvline(0, color='black', linestyle='--', linewidth=1)
    ax1.grid(axis='x', alpha=0.3)

    # Raw effect sizes
    colors2 = ['green' if v > 0 else ('red' if v < 0 else 'gray') for v in top['effect_size']]
    ax2.barh(range(top_n), top['effect_size'], color=colors2, alpha=0.8)
    ax2.set_yticks(range(top_n))
    ax2.set_yticklabels(top['gene_family'], fontsize=8)
    ax2.set_xlabel('Effect Size (Delta Diazotroph Rate)', fontsize=12)
    ax2.set_title('Feature Effect Sizes\n(present minus absent)', fontsize=12)
    ax2.axvline(0, color='black', linestyle='--', linewidth=1)
    ax2.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved feature effects plot to: {output_file}")


def create_feature_summary(direction_df, top_n=50):
    """Create a CSV summary for top features with quick interpretations."""
    top_n = min(top_n, len(direction_df))
    top = direction_df.head(top_n).copy()

    def interp(row):
        sign = '+' if row['effect_size'] > 0 else ''
        pct = row['effect_size'] * 100.0
        if row['direction'] == 'Neutral':
            return "No difference in diazotroph rate"
        return f"Presence associated with {('diazotrophy' if row['effect_size']>0 else 'non-diazotrophy')} ({sign}{pct:.1f}%-pts)"

    top['interpretation'] = top.apply(interp, axis=1)

    cols = [
        'gene_family', 'importance', 'direction', 'effect_size',
        'diaz_rate_present', 'diaz_rate_absent', 'n_genomes_with_feature',
        'p_value', 'interpretation'
    ]
    top[cols].to_csv('feature_directionality_summary.csv', index=False)
    print("Saved feature directionality summary to: feature_directionality_summary.csv")
    return top


# ---------------------- main ----------------------

def main():
    print("="*80)
    print("STEP 5: FEATURE DIRECTIONALITY ANALYSIS")
    print("="*80)

    # Load presence/absence
    print("\nLoading data...")
    X_full = pd.read_csv('gene_family_matrix.csv', index_col=0)

    # Metadata, prefer filtered set with proteins (like step 4)
    if os.path.exists('complete_genomes_with_proteins.csv'):
        meta = pd.read_csv('complete_genomes_with_proteins.csv')
    else:
        meta = pd.read_csv('complete_genomes_labeled.csv')

    # --- De-duplicate BEFORE intersecting ---
    if meta.duplicated('assembly_accession').any():
        dup = int(meta.duplicated('assembly_accession').sum())
        print(f"  Note: de-duplicating metadata on assembly_accession (dropping {dup} duplicates)")
        meta = meta.drop_duplicates('assembly_accession', keep='first')

    # Intersection and deterministic order
    common = sorted(set(meta['assembly_accession']) & set(X_full.index))
    meta = (meta[meta['assembly_accession'].isin(common)]
            .set_index('assembly_accession')
            .loc[common])
    X_full = X_full.loc[common]

    # Filter to informative families (same as step 4)
    purity = pd.read_csv('gene_family_purity_stats.csv')
    informative = purity.loc[~purity['is_pure'], 'gene_family'].tolist()
    X = X_full[informative]
    y = meta['is_diazotroph'].astype(int)

    print(f"  Analyzing {X.shape[1]} informative gene families across {X.shape[0]} genomes")

    # Load per-model importances if available
    imports = []
    for model_tag in ['Random_Forest', 'Gradient_Boosting', 'Logistic_Regression', 'XGBoost']:
        fn = f'feature_importance_{model_tag}.csv'
        if os.path.exists(fn):
            df = pd.read_csv(fn)
            if {'gene_family', 'importance'}.issubset(df.columns):
                df['model'] = model_tag
                imports.append(df)
    feature_importance_df = pd.concat(imports, ignore_index=True) if imports else None
    if feature_importance_df is None:
        print("Warning: No feature importance files found (directionality will still be computed).")

    # Directionality
    direction_df = calculate_feature_directionality(X, y, feature_importance_df)
    direction_df.to_csv('feature_directionality_full.csv', index=False)
    print("Saved full directionality results to: feature_directionality_full.csv")

    # Plots and summary
    plot_feature_effects(direction_df, top_n=30, output_file='feature_effects.png')
    create_feature_summary(direction_df, top_n=50)

    # Quick console highlights
    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)
    pos = direction_df[direction_df['direction'] == 'Positive'].head(10)
    if len(pos):
        print(f"\nTop {len(pos)} Diazotroph-Associated Features:")
        for _, r in pos.iterrows():
            print(f"  {r['gene_family']}: +{r['effect_size']*100:.1f}%-pts when present")
    neg = direction_df[direction_df['direction'] == 'Negative'].head(10)
    if len(neg):
        print(f"\nTop {len(neg)} Non-Diazotroph-Associated Features:")
        for _, r in neg.iterrows():
            print(f"  {r['gene_family']}: {r['effect_size']*100:.1f}%-pts toward non-diazotrophy when present")

    print("\n" + "="*80)
    print("STEP 5 COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
