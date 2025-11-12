#!/usr/bin/env python3
"""
debug_download.py
Quick diagnostics for HPC download failures.

What it checks:
1) DNS + HTTPS to NCBI FTP host (small file)
2) HTTPS to NCBI Datasets API (status endpoint)
3) Presence of tools (datasets, curl, wget)
4) Ability to fetch assembly_summary files
5) Try a single known path with and without assembly suffix to show the difference
"""

import os, shutil, subprocess, sys, tempfile, pathlib

TEST_SMALL = "https://ftp.ncbi.nlm.nih.gov/robots.txt"
TEST_API   = "https://api.ncbi.nlm.nih.gov/datasets/v2/status"
# A classic assembly with clear names
TEST_ASM = "GCF_000005845.2"  # E. coli K-12
TEST_RIGHT = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_protein.faa.gz"
TEST_WRONG = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2/GCF_000005845.2_protein.faa.gz"

def have(cmd):
    try:
        r = subprocess.run([cmd,"--help"], capture_output=True, text=True)
        return r.returncode in (0,1) or r.stdout or r.stderr
    except Exception:
        return False

def curl(url, dest):
    return subprocess.run(["curl","-L","--fail","-o",dest,url]).returncode == 0

def main():
    print("=== Tool availability ===")
    for t in ("datasets","curl","wget","python","conda","mamba"):
        print(f"{t:12s}:", "YES" if have(t) else "no")
    print("\n=== Network smoke tests ===")
    ok1 = curl(TEST_SMALL, "robots.txt")
    print("HTTPS to ftp.ncbi.nlm.nih.gov:", "OK" if ok1 else "FAIL")
    ok2 = curl(TEST_API, "datasets_status.json")
    print("HTTPS to api.ncbi.nlm.nih.gov :", "OK" if ok2 else "FAIL")
    if ok1: os.remove("robots.txt")
    if ok2: os.remove("datasets_status.json")

    print("\n=== Fetch assembly summaries ===")
    ok3 = curl("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt","assembly_summary_refseq.txt")
    ok4 = curl("https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt","assembly_summary_genbank.txt")
    print("refseq summary:", "OK" if ok3 else "FAIL")
    print("genbank summary:", "OK" if ok4 else "FAIL")
    for f in ("assembly_summary_refseq.txt","assembly_summary_genbank.txt"):
        if os.path.exists(f): os.remove(f)

    print("\n=== Illustrate correct vs incorrect file path ===")
    okR = curl(TEST_RIGHT, "ok.faa.gz")
    okW = curl(TEST_WRONG, "wrong.faa.gz")
    print("Correct (with assembly name suffix):", "OK" if okR else "FAIL")
    print("Incorrect (missing suffix)        :", "OK" if okW else "FAIL")
    for f in ("ok.faa.gz","wrong.faa.gz"):
        if os.path.exists(f): os.remove(f)

    print("\nDone. If API is OK but FTP is blocked, install NCBI Datasets CLI and use --prefer datasets-cli.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
