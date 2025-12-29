#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def load_rep_to_gf(gene_family_info_csv: str) -> dict:
    rep_to_gf = {}
    with open(gene_family_info_csv, newline="") as f:
        r = csv.DictReader(f)
        required = {"gene_family", "representative"}
        if not required.issubset(set(r.fieldnames or [])):
            raise ValueError(f"{gene_family_info_csv} must contain columns: {sorted(required)}")

        for row in r:
            rep = row["representative"]
            gf  = row["gene_family"]
            rep_to_gf[rep] = gf
    return rep_to_gf

def member_to_protein_id(member: str, keep_genome_prefix: bool) -> str:
    # member looks like: genome|protein_accession (because you wrote it that way in all_proteins.faa)
    if keep_genome_prefix:
        return member
    if "|" in member:
        return member.split("|", 1)[1]
    return member

def main():
    ap = argparse.ArgumentParser(description="Make protein->gene_family mapping from MMseqs clusters.")
    ap.add_argument("--clusters", default="gene_families_clusters.tsv",
                    help="MMseqs createtsv output (rep<TAB>member).")
    ap.add_argument("--family-info", default="gene_family_info.csv",
                    help="CSV with columns gene_family, representative (from Step 3).")
    ap.add_argument("--out", default="protein_family_map.tsv",
                    help="Output TSV with 2 columns: protein_id<TAB>gene_family")
    ap.add_argument("--keep-genome-prefix", action="store_true",
                    help="Use full member id (genome|protein) instead of only protein accession.")
    args = ap.parse_args()

    rep_to_gf = load_rep_to_gf(args.family_info)

    clusters_path = Path(args.clusters)
    if not clusters_path.exists():
        raise FileNotFoundError(f"Clusters file not found: {args.clusters}")

    n_in = 0
    n_out = 0
    n_skipped = 0

    with open(args.out, "w", newline="") as out_f, open(args.clusters, "r") as in_f:
        # Write header (comment out if you prefer headerless)
        out_f.write("protein_id\tgene_family\n")

        for line in in_f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            rep, member = parts[0], parts[1]
            n_in += 1

            gf = rep_to_gf.get(rep)
            if gf is None:
                # This means the rep was filtered out (e.g., < min_genomes_per_family) and is not in gene_family_info.csv
                n_skipped += 1
                continue

            protein_id = member_to_protein_id(member, args.keep_genome_prefix)
            out_f.write(f"{protein_id}\t{gf}\n")
            n_out += 1

    print(f"Done.")
    print(f"  Input lines:   {n_in}")
    print(f"  Output pairs:  {n_out}")
    print(f"  Skipped reps (not in gene_family_info.csv): {n_skipped}")
    print(f"  Wrote: {args.out}")

if __name__ == "__main__":
    main()
