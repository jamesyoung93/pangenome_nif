import os
import time, shutil, subprocess, urllib.request
from pathlib import Path

HMMSEARCH_CMD = os.environ.get("HMMSEARCH_CMD", "hmmsearch")

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def ftp_to_https(ftp_url: str) -> str:
    if not ftp_url:
        return ""
    return ftp_url.replace("ftp://", "https://")

def http_download(url: str, out_path: Path, timeout: int = 60, max_retries: int = 3):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp, open(out_path, "wb") as f:
                shutil.copyfileobj(resp, f)
            return out_path
        except Exception as e:
            last_err = e
            time.sleep(2 * attempt)
    raise RuntimeError(f"Failed to download after {max_retries} tries: {url} ; last error: {last_err}")

def run_cmd(cmd, check=True, capture=False):
    if capture:
        cp = subprocess.run(cmd, check=check, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return cp.stdout
    else:
        subprocess.run(cmd, check=check)

def hmmsearch_available():
    try:
        subprocess.run([HMMSEARCH_CMD, "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except Exception:
        return False

def hmmbuild_available():
    try:
        subprocess.run(["hmmbuild", "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except Exception:
        return False

def open_text(path: Path):
    return open(path, "rt", encoding="utf-8", errors="ignore")

def parse_tblout(tblout_path: Path):
    """
    Parse HMMER --tblout (3.x). Return list of dicts with keys:
    target, query, i_evalue, hmm_name, description.
    """
    rows = []
    if not tblout_path.exists():
        return rows
    hmm_name = tblout_path.stem
    with open_text(tblout_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 18:
                continue
            target = parts[0]          # target name
            query_name = parts[3]      # query name
            i_eval = parts[12]         # i-Evalue
            desc = " ".join(parts[19:]) if len(parts) > 19 else ""
            try:
                i_eval_val = float(i_eval) if i_eval not in ("*", "NA") else float("inf")
            except Exception:
                i_eval_val = float("inf")
            rows.append({
                "target": target,
                "query": query_name,
                "i_evalue": i_eval_val,
                "hmm_name": hmm_name,
                "description": desc
            })
    return rows


def discover_refs(ref_root="refs"):
    from pathlib import Path
    import os
    out = {"nifH": [], "nifD": [], "nifK": []}

    # Preferred layout: refs/hmms/<sub>/*.hmm
    for sub in ("nifH","nifD","nifK"):
        primary = sorted((Path(ref_root)/"hmms"/sub).glob("*.hmm"))
        if primary:
            files = primary
        else:
            # Fallback: any *.hmm under refs that contains the subunit tag in filename
            files = [p for p in Path(ref_root).rglob("*.hmm") if sub.lower() in p.name.lower()]
        out[sub] = [{"subunit": sub, "name": p.stem, "hmm": str(p), "seed": None} for p in files]

    # Absolute last resort: dump any *.hmm into nifH to avoid crash (still warns)
    if not any(out.values()):
        any_h = list(Path(ref_root).rglob("*.hmm"))
        if any_h:
            out["nifH"] = [{"subunit":"nifH","name":p.stem,"hmm":str(p),"seed":None} for p in any_h]
    if not any(out.values()):
        raise RuntimeError(f"No .hmm files under {ref_root}. Run setup_refs.sh")

    return out
# --- BEGIN: discover_refs for nifH/D/K (appended) ---
def discover_refs(ref_root="refs"):
    from pathlib import Path
    ref_root = Path(ref_root)
    out = {"nifH": [], "nifD": [], "nifK": []}

    # Preferred layout
    for sub in ("nifH","nifD","nifK"):
        hmm_dir  = ref_root / "hmms"  / sub
        seed_dir = ref_root / "seeds" / sub

        # Existing HMMs
        if hmm_dir.is_dir():
            for h in sorted(hmm_dir.glob("*.hmm")):
                out[sub].append({
                    "type": "hmm",
                    "subunit": sub,
                    "name": h.stem,
                    "hmm": str(h),
                    "seed": None
                })

        # Optional seeds (allow pipeline to build missing HMMs)
        if seed_dir.is_dir():
            seeds = list(seed_dir.glob("*.fa")) + list(seed_dir.glob("*.fasta"))
            for s in sorted(seeds):
                target_hmm = str((ref_root / "hmms" / sub / (s.stem + ".hmm")))
                out[sub].append({
                    "type": "seed",
                    "subunit": sub,
                    "name": s.stem,
                    "hmm": target_hmm,
                    "seed": str(s)
                })

    # Fallback: any *.hmm under refs whose filename mentions the subunit tag
    if not any(out.values()):
        for sub in ("nifH","nifD","nifK"):
            for h in ref_root.rglob("*.hmm"):
                if sub.lower() in h.name.lower():
                    out[sub].append({
                        "type": "hmm",
                        "subunit": sub,
                        "name": h.stem,
                        "hmm": str(h),
                        "seed": None
                    })

    # If absolutely nothing found, abort with a clear error
    if not any(out.values()):
        raise RuntimeError(f"No HMMs or seeds found under {ref_root}. Run setup_refs.sh with --source and --commit.")

    return out
# --- END ---
# --- BEGIN: discover_refs for nifH/D/K (appended) ---
def discover_refs(ref_root="refs"):
    from pathlib import Path
    ref_root = Path(ref_root)
    out = {"nifH": [], "nifD": [], "nifK": []}

    # Preferred layout
    for sub in ("nifH","nifD","nifK"):
        hmm_dir  = ref_root / "hmms"  / sub
        seed_dir = ref_root / "seeds" / sub

        # Existing HMMs
        if hmm_dir.is_dir():
            for h in sorted(hmm_dir.glob("*.hmm")):
                out[sub].append({
                    "type": "hmm",
                    "subunit": sub,
                    "name": h.stem,
                    "hmm": str(h),
                    "seed": None
                })

        # Optional seeds (allow pipeline to build missing HMMs)
        if seed_dir.is_dir():
            seeds = list(seed_dir.glob("*.fa")) + list(seed_dir.glob("*.fasta"))
            for s in sorted(seeds):
                target_hmm = str((ref_root / "hmms" / sub / (s.stem + ".hmm")))
                out[sub].append({
                    "type": "seed",
                    "subunit": sub,
                    "name": s.stem,
                    "hmm": target_hmm,
                    "seed": str(s)
                })

    # Fallback: any *.hmm under refs whose filename mentions the subunit tag
    if not any(out.values()):
        for sub in ("nifH","nifD","nifK"):
            for h in ref_root.rglob("*.hmm"):
                if sub.lower() in h.name.lower():
                    out[sub].append({
                        "type": "hmm",
                        "subunit": sub,
                        "name": h.stem,
                        "hmm": str(h),
                        "seed": None
                    })

    # If absolutely nothing found, abort with a clear error
    if not any(out.values()):
        raise RuntimeError(f"No HMMs or seeds found under {ref_root}. Run setup_refs.sh with --source and --commit.")

    return out
# --- END ---

# --- BEGIN: compat discover_refs with 'path' + 'hmm' ---
def discover_refs(ref_root="refs"):
    # Returns dict { 'nifH': [...], 'nifD': [...], 'nifK': [...] }.
    # Each item has at least:
    #   type: 'hmm' or 'seed'
    #   subunit: 'nifH'|'nifD'|'nifK'
    #   name: file stem
    #   path: path to HMM (for type=='hmm', existing HMM; for seed, target HMM to be built)
    #   hmm:  same as 'path' (alias for code that expects 'hmm')
    #   seed: (only for type=='seed') seed FASTA path if present
    from pathlib import Path as _P
    ref_root = _P(ref_root)
    out = {"nifH": [], "nifD": [], "nifK": []}

    def add_hmm(sub, p):
        p = _P(p)
        out[sub].append({
            "type": "hmm",
            "subunit": sub,
            "name": p.stem,
            "path": str(p),   # what build_missing_hmms() uses
            "hmm":  str(p),   # alias used elsewhere
        })

    # Preferred layout: refs/hmms/<sub>/*.hmm and optional refs/seeds/<sub>/*
    for sub in ("nifH","nifD","nifK"):
        hmm_dir  = ref_root / "hmms"  / sub
        seed_dir = ref_root / "seeds" / sub

        if hmm_dir.is_dir():
            for h in sorted(hmm_dir.glob("*.hmm")):
                add_hmm(sub, h)

        if seed_dir.is_dir():
            for s in sorted(list(seed_dir.glob("*.fa")) + list(seed_dir.glob("*.fasta"))):
                target = hmm_dir / (s.stem + ".hmm")
                out[sub].append({
                    "type": "seed",
                    "subunit": sub,
                    "name": s.stem,
                    "seed": str(s),
                    "path": str(target),  # hmmbuild output
                    "hmm":  str(target),
                })

    # Fallback: any *.hmm under refs/** that matches subunit in filename
    if not any(out.values()):
        for h in ref_root.rglob("*.hmm"):
            hn = h.name.lower()
            for sub in ("nifH","nifD","nifK"):
                if sub.lower() in hn:
                    add_hmm(sub, h)

    if not any(out.values()):
        raise RuntimeError(f"No HMMs or seeds under {ref_root}. Put HMMs in refs/hmms/nifH|nifD|nifK/*.hmm (or seeds under refs/seeds/<sub>/).")

    return out
# --- END: compat discover_refs ---
