#!/usr/bin/env python3
"""
02_download_proteins.py  (HOTFIX)
Robust, HPC-friendly downloader for NCBI protein FASTA per assembly accession.

Key changes in this hotfix:
- **Correct file paths**: Uses NCBI *assembly_summary* tables to build the exact
  per-assembly FTP/HTTPS paths, including the assembly name suffix in filenames
  (e.g., GCF_000005845.2_ASM584v2_protein.faa.gz).
- **Prefers NCBI Datasets** when available (CLI or REST) and gracefully falls back
  to HTTPS to ftp.ncbi.nlm.nih.gov using assembly_summary.
- **Handles GCF (RefSeq) vs GCA (GenBank)** and falls back to *_translated_cds.faa.gz
  when *_protein.faa.gz is absent.
- **HPC-friendly**: single outbound HTTPS only, retry logic, rate limiting, logs.
- **Input flexible**: takes either a one-per-line accession list, a CSV with an
  `assembly_accession` column, **or the canonical nif_hdk_hits CSV** with
  reproducible filtering to `Complete Genome` rows (and `download_ok=1`).
- **Output**: writes FASTA files to `proteins/` as `<accession>.faa` (gunzipped).

Usage examples
--------------
# A) Use nif_hdk_hits_enriched_with_quality_checkm.csv directly
python 02_download_proteins.py --nif-csv nif_hdk_hits_enriched_with_quality_checkm.csv

# B) Using accessions.txt (one per line)
python 02_download_proteins.py --accessions accessions.txt --outdir downloads --prefer datasets-cli

# C) Using a generic CSV and filtering to Complete Genome
python 02_download_proteins.py --csv labeled_complete_genomes.csv --complete-only

# D) Limit to 8 concurrent downloads and slow the rate
python 02_download_proteins.py --nif-csv nif_hdk_hits_enriched_with_quality_checkm.csv -j 8 --sleep 0.5

Notes
-----
- If your cluster blocks outbound traffic from compute nodes, run this on the
  **login** or **data transfer node (DTN)**, then point step 3 to the resulting
  `proteins/` directory.
- As of mid-2024, NCBI encourages using **NCBI Datasets** rather than hand-rolled
  FTP paths. This script supports both approaches.
"""

import argparse
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

USER_AGENT = "DiazotrophPipeline/1.0 ({email})"
DEFAULT_EMAIL = "james.young@jacks.sdstate.edu"

# Assembly summary sources
REFSEQ_SUMMARY_URL  = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
GENBANK_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"

def sh(cmd: List[str], check=False, capture=True, text=True, env=None, timeout=None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=check, capture_output=capture, text=text, env=env, timeout=timeout)

def have(cmd: str) -> bool:
    try:
        r = sh([cmd, "--help"], capture=True)
        return r.returncode in (0, 1) or (r.stdout or r.stderr)
    except Exception:
        return False

def read_accessions_from_file(path: Path) -> List[str]:
    accs = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"): continue
            accs.append(s)
    return accs

def read_accessions_from_csv(path: Path, accession_col: str = "assembly_accession") -> List[str]:
    accs = []
    with open(path, newline="") as fh:
        sniffer = csv.Sniffer()
        sample = fh.read(2048)
        fh.seek(0)
        dialect = sniffer.sniff(sample) if sample else csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        if accession_col not in reader.fieldnames:
            raise SystemExit(f"ERROR: column '{accession_col}' not found in CSV header: {reader.fieldnames}")
        for row in reader:
            acc = (row.get(accession_col) or "").strip()
            if not acc:
                continue
            # Optional filter: Complete Genome
            level = (row.get("assembly_level") or "").strip()
            yield (acc, level)
    # not reached


def read_complete_genomes_from_nif(path: Path) -> List[str]:
    """Parse nif_hdk_hits CSV to reproducibly capture target Complete Genomes.

    Selection rules (production-safe):
    - assembly_level == "Complete Genome"
    - download_ok is truthy ("1", "true", "yes", or empty/None)
    """

    required_cols = {"assembly_accession", "assembly_level"}
    with open(path, newline="") as hdr:
        header_reader = csv.DictReader(hdr)
        missing = required_cols - set(header_reader.fieldnames or [])
    if missing:
        raise SystemExit(f"ERROR: missing required columns in {path}: {missing}")

    accessions: List[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            level = (row.get("assembly_level") or "").strip()
            if level != "Complete Genome":
                continue

            download_ok_raw = (row.get("download_ok") or "").strip().lower()
            if download_ok_raw not in ("", "1", "true", "yes", "y"):
                continue

            acc = (row.get("assembly_accession") or "").strip()
            if acc:
                accessions.append(acc)

    deduped = sorted(dict.fromkeys(accessions))
    if not deduped:
        raise SystemExit(
            f"ERROR: no Complete Genome rows with download_ok in {path}. Check the input file."
        )

    return deduped

def filter_complete_only(items: Iterable[Tuple[str, str]]) -> List[str]:
    out = []
    for acc, level in items:
        if level == "Complete Genome":
            out.append(acc)
    return out

def curl(url: str, dest: Path, email: str, tries: int = 3, sleep: float = 0.5) -> bool:
    # Use curl in a conservative, HPC-friendly way
    ua = USER_AGENT.format(email=email or DEFAULT_EMAIL)
    cmd = [
        "curl", "-L", "--fail", "--retry", str(tries), "--retry-delay", "2",
        "-A", ua, "-o", str(dest), url
    ]
    r = sh(cmd)
    if r.returncode != 0:
        return False
    # sanity
    try:
        if dest.stat().st_size == 0:
            return False
    except FileNotFoundError:
        return False
    time.sleep(sleep)
    return True

def gunzip_to(src_gz: Path, dest_faa: Path) -> bool:
    try:
        with gzip.open(src_gz, "rb") as fin, open(dest_faa, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        return True
    except Exception:
        return False

def load_assembly_summaries(email: str) -> Tuple[Dict[str,str], Dict[str,str]]:
    """
    Returns two dicts:
      refseq_map[ACCESSION] -> ftp_path
      genbank_map[ACCESSION] -> ftp_path
    ftp_path examples end with the assembly name directory, e.g.:
      https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2
    """
    tmp = Path(".tmp_download")
    tmp.mkdir(exist_ok=True)
    refseq_file  = tmp / "assembly_summary_refseq.txt"
    genbank_file = tmp / "assembly_summary_genbank.txt"
    ok1 = curl(REFSEQ_SUMMARY_URL, refseq_file, email=email, tries=3, sleep=0.1)
    ok2 = curl(GENBANK_SUMMARY_URL, genbank_file, email=email, tries=3, sleep=0.1)
    if not (ok1 and ok2):
        print("WARNING: failed to fetch one or both assembly_summary files. Will try Datasets first.")
    def parse(path: Path) -> Dict[str,str]:
        m = {}
        if not path.exists():
            return m
        with open(path) as fh:
            for line in fh:
                if not line or line.startswith("#"): 
                    continue
                # assembly_summary columns: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 20:
                    continue
                acc = parts[0].strip()   # assembly_accession
                ftp = parts[19].strip()  # ftp_path
                if ftp and ftp != "na":
                    m[acc] = ftp.replace("ftp://", "https://")  # prefer HTTPS
        return m
    return parse(refseq_file), parse(genbank_file)

def path_candidates(ftp_dir: str) -> List[str]:
    """
    Given ftp_dir: .../GCF_000005845.2_ASM584v2
    Return likely protein files in priority order.
    """
    tail = ftp_dir.rstrip("/").split("/")[-1]
    return [
        f"{ftp_dir}/{tail}_protein.faa.gz",
        f"{ftp_dir}/{tail}_translated_cds.faa.gz",
        f"{ftp_dir}/{tail}_protein.gpff.gz"
    ]

def download_one_via_summary(acc: str, outdir_prot: Path, maps: Tuple[Dict[str,str],Dict[str,str]], email: str, sleep: float, tries: int) -> Tuple[str, bool, str]:
    refseq_map, genbank_map = maps
    ftp = refseq_map.get(acc) if acc.startswith("GCF_") else genbank_map.get(acc)
    if not ftp:
        return acc, False, "no_ftp_path_in_assembly_summary"
    for url in path_candidates(ftp):
        tmp_gz = outdir_prot / f"{acc}.tmp.gz"
        ok = curl(url, tmp_gz, email=email, tries=tries, sleep=sleep)
        if not ok:
            if tmp_gz.exists():
                tmp_gz.unlink(missing_ok=True)
            continue
        # gunzip to .faa (or .gpff -> we still name .faa for pipeline compatibility)
        out_faa = outdir_prot / f"{acc}.faa"
        if gunzip_to(tmp_gz, out_faa):
            tmp_gz.unlink(missing_ok=True)
            return acc, True, "ok"
        else:
            tmp_gz.unlink(missing_ok=True)
    return acc, False, "not_found_or_failed"

def download_one_via_datasets_cli(acc: str, outdir_prot: Path, email: str, sleep: float) -> Tuple[str, bool, str]:
    """
    datasets download genome accession <acc> --include protein --filename <acc>.zip
    then unzip and copy ncbi_dataset/data/<acc>/**/protein.faa -> proteins/<acc>.faa
    """
    ua = USER_AGENT.format(email=email or DEFAULT_EMAIL)
    ztmp = outdir_prot / f"{acc}.zip"
    cmd = ["datasets", "download", "genome", "accession", acc, "--include", "protein", "--filename", str(ztmp), "--no-progressbar"]
    r = sh(cmd)
    if r.returncode != 0 or not ztmp.exists() or ztmp.stat().st_size == 0:
        return acc, False, "datasets_cli_failed"
    # Unzip shallowly and pick protein.faa
    import zipfile
    try:
        with zipfile.ZipFile(ztmp) as zf:
            candidates = [n for n in zf.namelist() if n.endswith("/protein.faa")]
            if not candidates:
                return acc, False, "datasets_zip_no_protein"
            # if more than one pick first
            data = zf.read(candidates[0])
            out_faa = outdir_prot / f"{acc}.faa"
            with open(out_faa, "wb") as out:
                out.write(data)
    except Exception as e:
        return acc, False, f"unzip_error:{e}"
    finally:
        ztmp.unlink(missing_ok=True)
    time.sleep(sleep)
    return acc, True, "ok"

def download_one_via_datasets_api(acc: str, outdir_prot: Path, email: str, sleep: float) -> Tuple[str, bool, str]:
    """
    Try the REST API using curl without requiring the CLI.
    Endpoint (v2): GET /datasets/v2/genome/accession/{acc}/download?include=protein
    """
    ztmp = outdir_prot / f"{acc}.zip"
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include=protein"
    ok = curl(url, ztmp, email=email, tries=3, sleep=sleep)
    if not ok:
        return acc, False, "datasets_api_failed"
    # Unzip and copy
    import zipfile
    try:
        with zipfile.ZipFile(ztmp) as zf:
            candidates = [n for n in zf.namelist() if n.endswith("/protein.faa")]
            if not candidates:
                return acc, False, "datasets_zip_no_protein"
            data = zf.read(candidates[0])
            out_faa = outdir_prot / f"{acc}.faa"
            with open(out_faa, "wb") as out:
                out.write(data)
    except Exception as e:
        return acc, False, f"unzip_error:{e}"
    finally:
        ztmp.unlink(missing_ok=True)
    return acc, True, "ok"

def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--accessions", type=Path, help="Text file with one assembly accession per line")
    src.add_argument("--csv", type=Path, help="CSV containing an 'assembly_accession' column (optionally 'assembly_level')")
    src.add_argument("--nif-csv", dest="nif_csv", type=Path,
                    help="Use nif_hdk_hits CSV; automatically filters to Complete Genome rows with download_ok=1")
    p.add_argument("--complete-only", action="store_true", help="If using --csv, keep only rows with assembly_level == 'Complete Genome'")
    p.add_argument("--outdir", type=Path, default=Path("downloads"), help="Directory to place downloads and logs")
    p.add_argument("--email", default=DEFAULT_EMAIL, help="Email to include in User-Agent (recommended by NCBI)")
    p.add_argument("-j", "--jobs", type=int, default=8, help="Parallel downloads")
    p.add_argument("--sleep", type=float, default=0.25, help="Sleep between requests (seconds)")
    p.add_argument("--prefer", choices=["datasets-cli","datasets-api","assembly-summary"], default="datasets-cli",
                   help="Preferred primary method; script falls back automatically")
    p.add_argument("--retries", type=int, default=3, help="Retries per file (curl)")
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    proteins_dir = args.outdir / "proteins"
    proteins_dir.mkdir(exist_ok=True)

    # Build accession list
    source_desc = ""
    if args.nif_csv:
        accessions = read_complete_genomes_from_nif(args.nif_csv)
        source_desc = f"{args.nif_csv} (Complete Genome + download_ok=1)"
    elif args.accessions:
        accessions = read_accessions_from_file(args.accessions)
        source_desc = str(args.accessions)
    else:
        pairs = list(read_accessions_from_csv(args.csv))
        accessions = filter_complete_only(pairs) if args.complete_only else [a for a,_ in pairs]
        source_desc = f"{args.csv} ({'Complete Genome only' if args.complete_only else 'all rows'})"

    accessions = sorted(dict.fromkeys(accessions))

    if not accessions:
        print("No accessions found. Exiting.")
        sys.exit(0)

    manifest = args.outdir / "target_complete_genomes.txt"
    with open(manifest, "w") as mf:
        for acc in accessions:
            mf.write(f"{acc}\n")
    print(f"Target genomes: {len(accessions)} (source: {source_desc})")
    print(f"Manifest: {manifest}")

    # Tool detection
    have_datasets = have("datasets")
    have_curl = have("curl")

    # Load assembly summaries for the fallback
    maps = load_assembly_summaries(args.email)

    log_path = args.outdir / "download_log.tsv"
    with open(log_path, "w") as log:
        print("accession\tstatus\tmethod\tmessage", file=log)

        def task(acc: str) -> Tuple[str, str, str]:
            # Try the chosen method, then fall back
            methods = []
            if args.prefer == "datasets-cli" and have_datasets:
                methods = [("datasets-cli", download_one_via_datasets_cli)]
            elif args.prefer == "datasets-api" and have_curl:
                methods = [("datasets-api", lambda a, d, e, s: download_one_via_datasets_api(a, d, e, s))]
            else:
                methods = []
            # Always include the other options as fallbacks if available
            if have_curl:
                methods.append(("datasets-api", lambda a, d, e, s: download_one_via_datasets_api(a, d, e, s)))
            methods.append(("assembly-summary", lambda a, d, e, s: download_one_via_summary(a, d, maps, e, s, args.retries)))

            for mname, mfun in methods:
                acc_, ok, msg = mfun(acc, proteins_dir, args.email, args.sleep) if "datasets" in mname else mfun(acc, proteins_dir, args.email, args.sleep)
                if ok:
                    return acc, "ok", mname
            return acc, "failed", msg

        futures = {}
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
            for acc in accessions:
                futures[ex.submit(task, acc)] = acc
            ok_count = 0
            fail_count = 0
            for fut in as_completed(futures):
                acc, status, info = fut.result()
                if status == "ok":
                    ok_count += 1
                    print(f"{acc}\tok\t{info}\t", file=log)
                else:
                    fail_count += 1
                    print(f"{acc}\tfail\t{info}\t", file=log)
                print(f"[{ok_count}/{ok_count+fail_count}] {acc}: {status} ({info})")

    print("\nSummary: wrote proteins to", proteins_dir)
    print("Log:", log_path)
    print("Done.")

if __name__ == "__main__":
    main()
