#!/usr/bin/env python3
from util import discover_refs
import sys, os, csv, gzip
from pathlib import Path
import yaml
from tqdm import tqdm
import pandas as pd
from Bio import Entrez
from util import (
    ensure_dir, ftp_to_https, http_download, run_cmd,
    hmmsearch_available, hmmbuild_available, parse_tblout,
    HMMSEARCH_CMD
)

ROOT = Path(__file__).resolve().parent

def load_cfg():
    with open(ROOT / "config.yaml","r") as fh:
        return yaml.safe_load(fh)

def entrez_setup():
    email = os.getenv("ENTREZ_EMAIL", "")
    if not email:
        print("ERROR: Please set ENTREZ_EMAIL (required by NCBI).", file=sys.stderr)
        sys.exit(2)
    Entrez.email = email
    api_key = os.getenv("NCBI_API_KEY", None)
    if api_key:
        Entrez.api_key = api_key

def query_ncbi():
    cfg = load_cfg()
    entrez_setup()
    out_dir = ensure_dir(ROOT / cfg["paths"]["results"] / "assemblies")
    term = cfg["taxonomy"]["taxon_query"]

    def _assembly_levels():
        env_levels = os.getenv("NIF_ASSEMBLY_LEVELS", "").strip()
        if env_levels:
            parts = [p.strip() for p in env_levels.split(",") if p.strip()]
            if parts:
                print(f"Overriding assembly levels via NIF_ASSEMBLY_LEVELS -> {parts}", file=sys.stderr)
                return parts
        return cfg["taxonomy"]["assembly_level"]

    levels = _assembly_levels()
    refseq_only = cfg["taxonomy"]["refseq_only"]

    level_term = " OR ".join([f'"{lvl}"[Assembly Level]' for lvl in levels])
    term_full = f"{term} AND ({level_term})"
    if refseq_only:
        term_full += " AND refseq[filter]"
    print(f"NCBI Entrez query: {term_full}", file=sys.stderr)

    with Entrez.esearch(db="assembly", term=term_full, retmax=100000) as h:
        sr = Entrez.read(h)
    ids = sr.get("IdList", [])
    print(f"Found {len(ids)} assembly IDs", file=sys.stderr)

    rows = []
    B = 250
    for i in tqdm(range(0, len(ids), B), desc="Fetch docsum"):
        chunk = ids[i:i+B]
        with Entrez.esummary(db="assembly", id=",".join(chunk), report="full") as h:
            summ = Entrez.read(h, validate=False)
        for doc in summ["DocumentSummarySet"]["DocumentSummary"]:
            rows.append({
                "assembly_accession": doc["AssemblyAccession"],
                "organism":          doc["Organism"],
                "taxid":             doc["Taxid"],
                "ftp_refseq":        doc.get("FtpPath_RefSeq",""),
                "ftp_genbank":       doc.get("FtpPath_GenBank","")
            })

    df = pd.DataFrame(rows).drop_duplicates("assembly_accession")
    # avoid NaN -> strings
    for col in ("ftp_refseq","ftp_genbank"):
        if col in df.columns:
            df[col] = df[col].fillna("")
    df.to_csv(out_dir / "assemblies.tsv", sep='\t', index=False)
    print(f"Wrote {len(df)} assemblies -> {out_dir/'assemblies.tsv'}")

def _pick_ftp(row: dict) -> str:
    # return a usable string or empty
    for key in ("ftp_refseq","ftp_genbank"):
        val = row.get(key, "")
        if isinstance(val, str) and val.strip() and val.lower() != "na":
            return val
    return ""

def download_proteomes():
    cfg = load_cfg()
    asm_path = ROOT / cfg["paths"]["assemblies_tsv"]
    if not asm_path.exists():
        print(f"ERROR: {asm_path} not found. Run: python pipeline.py query", file=sys.stderr)
        sys.exit(2)

    asm = pd.read_csv(asm_path, sep='\t', dtype=str)
    out_dir = ensure_dir(ROOT / cfg["download"]["out_dir"])
    timeout  = int(cfg["download"]["timeout_sec"])
    retries  = int(cfg["download"]["max_retries"])

    manifest_rows = []
    for _, r in tqdm(asm.iterrows(), total=len(asm), desc="Download protein.faa.gz"):
        acc  = r["assembly_accession"]
        base = _pick_ftp(r)
        if not base:
            manifest_rows.append({"assembly_accession": acc, "status": "no_ftp"})
            continue
        url  = ftp_to_https(base) + "/" + Path(base).name + "_protein.faa.gz"
        dest = out_dir / f"{acc}_protein.faa.gz"
        if not dest.exists() or dest.stat().st_size == 0:
            try:
                http_download(url, dest, timeout=timeout, max_retries=retries)
                status = "ok"
            except Exception as e:
                status = f"error:{e}"
        else:
            status = "cached"
        manifest_rows.append({"assembly_accession": acc, "url": url, "path": str(dest), "status": status})

    mdf = pd.DataFrame(manifest_rows)
    mdf.to_csv(out_dir / "download_manifest.tsv", sep='\t', index=False)
    print(f"Wrote manifest -> {out_dir/'download_manifest.tsv'}")

def combine_proteins():
    import os, gzip, glob, shutil
    from pathlib import Path
    import pandas as pd

    asm_tsv = Path("results/assemblies/assemblies.tsv")
    out_dir = Path("results/proteomes")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_faa = out_dir / "combined_proteins.faa"
    manifest = out_dir / "combined_manifest.tsv"
    missing = out_dir / "missing_proteomes.tsv"

    # Always read as strings; kill NaNs
    df = pd.read_csv(asm_tsv, sep="\t", dtype=str).fillna("")

    # Try to detect an existing path column; otherwise we resolve by accession
    cand_cols = [c for c in df.columns if c.lower() in
                 ("protein_faa_gz", "protein_faa", "faa_path", "proteome_faa")]
    if cand_cols:
        path_col = cand_cols[0]
    else:
        path_col = None

    used = []
    missing_rows = []

    def first_existing(acc, hinted_path):
        cands = []
        if hinted_path:
            cands.append(hinted_path)
        if acc:
            # common patterns
            cands += glob.glob(f"results/downloads/**/*{acc}*_protein.faa.gz", recursive=True)
            cands += glob.glob(f"results/downloads/**/*{acc}*_protein.faa", recursive=True)
        for p in cands:
            try:
                if os.path.getsize(p) > 0:
                    return p
            except OSError:
                pass
        return ""

    with open(out_faa, "wb") as w, open(manifest, "w") as man:
        print("assembly_accession\tpath", file=man)
        for _, row in df.iterrows():
            acc = (row.get("assembly_accession", "") or "").strip()
            hinted = (row.get(path_col, "") if path_col else "").strip()
            f = first_existing(acc, hinted)
            if not f:
                missing_rows.append((acc, hinted))
                continue
            print(f"{acc}\t{f}", file=man)
            if f.endswith(".gz"):
                with gzip.open(f, "rb") as r:
                    shutil.copyfileobj(r, w)
            else:
                with open(f, "rb") as r:
                    shutil.copyfileobj(r, w)
            used.append(f)

    with open(missing, "w") as m:
        print("assembly_accession\trequested_path", file=m)
        for acc, hint in missing_rows:
            if acc or hint:
                print(f"{acc}\t{hint}", file=m)

    print(f"[combine] wrote {out_faa} from {len(used)} proteomes; "
          f"{len(missing_rows)} missing ? {missing}")

def build_missing_hmms(refs):
    """
    If HMMs are absent but FASTA seeds are present, build HMMs with hmmbuild.
    Returns: dict{subunit: [hmm_paths]}
    """
    if not hmmbuild_available():
        return {k: [r["path"] for r in v if r["type"]=="hmm"] for k,v in refs.items()}

    built = {"nifH": [], "nifD": [], "nifK": []}
    out_base = ROOT / "refs/hmms"
    for sub, lst in refs.items():
        # existing HMMs
        for r in lst:
            if r["type"]=="hmm":
                built[sub].append(r["path"])
        # build from FASTA if needed
        fasta_only = [r["path"] for r in lst if r["type"]=="fasta"]
        for fa in fasta_only:
            hmm_out = out_base / sub / (Path(fa).stem + ".hmm")
            hmm_out.parent.mkdir(parents=True, exist_ok=True)
            if not hmm_out.exists() or hmm_out.stat().st_size==0:
                run_cmd(["hmmbuild", str(hmm_out), str(fa)], check=True)
            built[sub].append(hmm_out)
    return built

def run_hmmsearch_all(hmms_by_subunit):
    cfg = load_cfg()
    if not hmmsearch_available():
        print(
            "ERROR: hmmsearch is not available on PATH. Install HMMER or load your cluster module (e.g., 'module load hmmer/3.4').",
            file=sys.stderr,
        )
        sys.exit(2)
    combined = ROOT / cfg["paths"]["combined_proteins"]
    if not combined.exists():
        print(f"ERROR: {combined} not found. Run combine first.", file=sys.stderr)
        sys.exit(2)
    out_dir = ensure_dir(ROOT / cfg["scan"]["hmm_out_dir"])
    cpu = str(cfg["scan"]["cpu"])

    for sub, hmm_list in hmms_by_subunit.items():
        for hmm in hmm_list:
            stem = Path(hmm).stem
            tbl   = out_dir / f"{sub}__{stem}.tblout"
            domtb = out_dir / f"{sub}__{stem}.domtblout"
            cmd = [
                HMMSEARCH_CMD,"--cpu", cpu, "--noali",
                "--tblout", str(tbl),
                "--domtblout", str(domtb),
                str(hmm), str(combined)
            ]
            print("RUN:", " ".join(cmd), file=sys.stderr)
            run_cmd(cmd, check=True)

def summarize_hits():
    cfg = load_cfg()
    ecut   = float(cfg["scan"]["evalue_threshold"])
    hmm_dir= ROOT / cfg["scan"]["hmm_out_dir"]
    asm    = pd.read_csv(ROOT / cfg["paths"]["assemblies_tsv"], sep='\t', dtype=str)
    accs   = list(asm["assembly_accession"])

    subs   = ["nifH","nifD","nifK"]
    hits   = {s:{} for s in subs}

    for tbl in sorted(hmm_dir.glob("*.tblout")):
        base = tbl.name
        if "__" not in base:
            continue
        sub = base.split("__",1)[0]
        if sub not in subs:
            continue
        refname = base.split("__",1)[1].replace(".tblout","")
        for r in parse_tblout(tbl):
            tgt = r["target"]       # "{assembly}|{protein_id}"
            if "|" not in tgt:
                continue
            asm_acc, prot_id = tgt.split("|",1)
            i_eval = r["i_evalue"]
            cur = hits[sub].get(asm_acc)
            if (cur is None) or (i_eval < cur[0]):
                hits[sub][asm_acc] = (i_eval, prot_id, refname)

    out_rows = []
    for acc in accs:
        row = {"assembly_accession": acc}
        for s in subs:
            best = hits[s].get(acc)
            if best:
                best_eval, best_prot, best_ref = best
                row[f"{s}_present"]     = 1 if best_eval <= ecut else 0
                row[f"{s}_best_evalue"] = best_eval
                row[f"{s}_best_ref"]    = best_ref
            else:
                row[f"{s}_present"]     = 0
                row[f"{s}_best_evalue"] = ""
                row[f"{s}_best_ref"]    = ""
        ids = []
        for s in subs:
            best = hits[s].get(acc)
            ids.append(f"{s[-1].upper()}={best[1]}" if best else f"{s[-1].upper()}=")
        row["best_hit_ids"] = ";".join(ids)
        out_rows.append(row)

    out = pd.DataFrame(out_rows)
    ensure_dir(ROOT / "results")
    out.to_csv(ROOT / "results/nif_hdk_hits.csv", index=False)
    print(f"Wrote -> {ROOT/'results/nif_hdk_hits.csv'}")

def run_all():
    query_ncbi()
    download_proteomes()
    combine_proteins()
    refs = discover_refs()
    hmms = build_missing_hmms(refs)
    run_hmmsearch_all(hmms)
    summarize_hits()

def scan():
    refs = discover_refs()
    hmms = build_missing_hmms(refs)
    run_hmmsearch_all(hmms)

def main():
    if len(sys.argv) < 2:
        print("Usage: python pipeline.py [query|download|combine|scan|summarize|run-all]")
        sys.exit(1)
    cmd = sys.argv[1]
    if   cmd == "query":      query_ncbi()
    elif cmd == "download":   download_proteomes()
    elif cmd == "combine":    combine_proteins()
    elif cmd == "scan":       scan()
    elif cmd == "summarize":  summarize_hits()
    elif cmd == "run-all":    run_all()
    else:
        print(f"Unknown command: {cmd}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
