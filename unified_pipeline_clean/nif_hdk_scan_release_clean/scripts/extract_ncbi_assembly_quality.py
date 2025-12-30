#!/usr/bin/env python3
import argparse, json, re, csv
def find_first(obj, names):
    want={n.lower() for n in names}; stack=[obj]
    while stack:
        x=stack.pop()
        if isinstance(x,dict):
            for k,v in x.items():
                if k.lower() in want: return v
                if isinstance(v,(dict,list)): stack.append(v)
        elif isinstance(x,list):
            stack.extend(x)
    return None
ap=argparse.ArgumentParser()
ap.add_argument("--in",  dest="inp", required=True)
ap.add_argument("--out", dest="out", required=True)
a=ap.parse_args()
rows=[]
with open(a.inp, encoding="utf-8") as f:
    for line in f:
        if not line.strip(): continue
        j=json.loads(line)
        acc = j.get("accession") or find_first(j,{"assembly_accession","accession"}) or ""
        org = find_first(j,{"organism_name","organism-name","organism"}) or ""
        lvl = find_first(j,{"assembly_level","assembly-level"}) or ""
        contigs   = find_first(j,{"number_of_contigs","number-of-contigs","contig_count"})
        contig_n50= find_first(j,{"contig_n50","contig-n50"})
        scaffs    = find_first(j,{"number_of_scaffolds","number-of-scaffolds","scaffold_count"})
        scaffold_n50=find_first(j,{"scaffold_n50","scaffold-n50"})
        total_len = find_first(j,{"total_ungapped_length","ungapped_sequence_length","total_sequence_length"})
        gc        = find_first(j,{"gc_percent","gc-content","gcpercent"})
        cov       = find_first(j,{"coverage","genome_coverage","read_coverage"})
        ann       = find_first(j,{"annotation","annotation_info","annotation-metadata"}) or {}
        ann_name  = (ann.get("name") if isinstance(ann,dict) else None) or find_first(ann,{"name","annotation_name"})
        ann_date  = (ann.get("date") if isinstance(ann,dict) else None) or find_first(ann,{"date","annotation_date"})
        chk       = find_first(j,{"checkm","checkm_analysis","checkm-analysis"}) or {}
        chk_comp  = (chk.get("completeness") if isinstance(chk,dict) else None) or find_first(chk,{"completeness"})
        chk_cont  = (chk.get("contamination") if isinstance(chk,dict) else None) or find_first(chk,{"contamination"})
        notes     = j.get("genome_notes") or j.get("genomeNotes") or []
        if isinstance(notes,dict): notes=[notes]
        txt="; ".join({(n.get("message") or n.get("text") or "") if isinstance(n,dict) else str(n) for n in notes if n})
        rows.append({"accession":acc,"organism_name":org,"assembly_level":lvl,
            "number_of_contigs":contigs,"contig_n50":contig_n50,
            "number_of_scaffolds":scaffs,"scaffold_n50":scaffold_n50,
            "total_ungapped_length":total_len,"gc_percent":gc,"coverage":cov,
            "annotation_name":ann_name,"annotation_date":ann_date,
            "checkm_completeness":chk_comp,"checkm_contamination":chk_cont,
            "genome_notes":txt,
            "derived_from_metagenome": int(bool(re.search(r"derived from metagenome", txt or "", re.I))),
            "annotation_fails_MAG_check": int(bool(re.search(r"annotation fails MAG completeness check", txt or "", re.I))),
        })
cols=["accession","organism_name","assembly_level","number_of_contigs","contig_n50","number_of_scaffolds","scaffold_n50","total_ungapped_length","gc_percent","coverage","annotation_name","annotation_date","checkm_completeness","checkm_contamination","genome_notes","derived_from_metagenome","annotation_fails_MAG_check"]
with open(a.out,"w",newline="",encoding="utf-8") as o:
    w=csv.DictWriter(o,fieldnames=cols,extrasaction="ignore",delimiter="\t"); w.writeheader(); w.writerows(rows)
print("Wrote ->", a.out)
