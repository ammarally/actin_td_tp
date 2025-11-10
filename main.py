#!/usr/bin/env python3
"""
Viractin TD/TP end-to-end runner
Author: Ali Ammar (@ammarally)
Runs: DB build (if needed), PersoHMM scan, HMMER extraction, MSA, trimming, tree, and QC summary.

Usage:
  python scripts/main.py --threads 4
"""

import os, sys, shutil, subprocess, pathlib, argparse, re

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
MYDATA = ROOT / "myData"
HMMS = ROOT / "hmms" / "PersoHMM"
RESULTS = ROOT / "results"
LOGS = ROOT / "logs"
AA_DIR = RESULTS / "aa"
HMMER_DIR = RESULTS / "hmmer"
HITS_DIR = RESULTS / "hmm_hits"
QC_DIR = RESULTS / "qc"
REF_FASTA = ROOT / "ref" / "Actin_ARPs_reference.fasta"  # optional

MAG_LIST = DATA / "listeMAGs.txt"
PERSO_HMM = HMMS / "genes.hmm"

# default input contigs suffixes to try for each MAG
CONTIG_SUFFIXES = ["-contigs.fa", ".fa", ".fasta", ".fna"]

REQ_CMDS = [
    "anvi-gen-contigs-database",
    "anvi-run-hmms",
    "anvi-get-sequences-for-gene-calls",
    "hmmscan",
    "seqkit",
    # at least one aligner:
    # "muscle" or "mafft"
    "trimal",
    "iqtree",
]

def sh(cmd, log=None, check=True):
    """Run a shell command with streaming and optional log file."""
    print(f"\n$ {cmd}")
    proc = subprocess.Popen(cmd, shell=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out_lines = []
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        f = open(log, "a")
    else:
        f = None
    try:
        for line in proc.stdout:
            sys.stdout.write(line)
            out_lines.append(line)
            if f:
                f.write(line)
    finally:
        if f:
            f.close()
    code = proc.wait()
    if check and code != 0:
        raise RuntimeError(f"Command failed ({code}): {cmd}")
    return "".join(out_lines)

def which(cmd):
    return shutil.which(cmd)

def check_deps():
    missing = []
    for c in REQ_CMDS:
        if c == "anvi-run-hmms" and not PERSO_HMM.exists():
            print(f"[WARN] PersoHMM not found at {PERSO_HMM}. Download/index required (genes.hmm + .h3*).")
        if not which(c):
            missing.append(c)
    # Need at least one aligner
    aligner = "muscle" if which("muscle") else ("mafft" if which("mafft") else None)
    if aligner is None:
        missing.append("muscle_or_mafft")
    if missing:
        print("[ERROR] Missing dependencies:", ", ".join(missing))
        sys.exit(1)
    return "muscle" if which("muscle") else "mafft"

def mag_names():
    with open(MAG_LIST) as f:
        for line in f:
            name = line.strip()
            if name:
                yield name

def find_contigs_path(mag):
    for suff in CONTIG_SUFFIXES:
        p = DATA / f"{mag}{suff}"
        if p.exists():
            return p
    return None

def ensure_dirs():
    for d in [RESULTS, LOGS, AA_DIR, HMMER_DIR, HITS_DIR, QC_DIR]:
        d.mkdir(parents=True, exist_ok=True)

def build_contigs_dbs(overwrite=False):
    for mag in mag_names():
        cdb = MYDATA / f"{mag}-DB.db"
        if cdb.exists() and not overwrite:
            print(f"[SKIP] Contigs DB exists: {cdb.name}")
            continue
        contigs = find_contigs_path(mag)
        if not contigs:
            print(f"[WARN] No contigs file for {mag} (tried {CONTIG_SUFFIXES}) — skipping DB build.")
            continue
        log = LOGS / f"01_gen_contigs_{mag}.log"
        sh(f'anvi-gen-contigs-database -f "{contigs}" -o "{cdb}"', log=log)

def run_perso_hmm_on_each_db(overwrite=False, threads=1):
    if not PERSO_HMM.exists():
        print(f"[WARN] PersoHMM DB not found at {PERSO_HMM}. Ensure files genes.hmm + .h3* are present.")
    for mag in mag_names():
        cdb = MYDATA / f"{mag}-DB.db"
        if not cdb.exists():
            print(f"[SKIP] No DB for {mag}.")
            continue
        # anvi-run-hmms updates the contigs DB; no explicit output to detect, so allow re-run
        log = LOGS / "02_run_hmms.log"
        sh(f'anvi-run-hmms -c "{cdb}" -H "{HMMS}" -T {threads}', log=log)

def qc_completeness(overwrite=False):
    # generate per-MAG stats
    for mag in mag_names():
        cdb = MYDATA / f"{mag}-DB.db"
        if not cdb.exists():
            continue
        out = QC_DIR / f"{mag}_stats.txt"
        if out.exists() and not overwrite:
            print(f"[SKIP] QC exists: {out.name}")
            continue
        log = LOGS / "qc" / f"{mag}.log"
        sh(f'anvi-display-contigs-stats "{cdb}" --report-as-text -o "{out}"', log=log)
    # summary table
    summary = QC_DIR / "SCG_summary.tsv"
    with open(summary, "w") as S:
        S.write("MAG\tSCG_source\tSCG_found\tCompleteness\tRedundancy\n")
        for mag in mag_names():
            f = QC_DIR / f"{mag}_stats.txt"
            if not f.exists():
                S.write(f"{mag}\tNA\tNA\tNA\tNA\n")
                continue
            txt = f.read_text()
            # Parse lines if present in your anvio version
            src = re.search(r"Single-copy core gene source\s*:\s*(.+)", txt)
            found = re.search(r"Number of single-copy genes found\s*:\s*([0-9/]+)", txt)
            comp = re.search(r"Estimated completeness\s*:\s*([\d\.]+)\s*%", txt)
            red  = re.search(r"Estimated redundancy\s*:\s*([\d\.]+)\s*%", txt)
            S.write(f"{mag}\t{(src.group(1) if src else 'NA')}\t{(found.group(1) if found else 'NA')}\t{(comp.group(1) if comp else 'NA')}\t{(red.group(1) if red else 'NA')}\n")
    print(f"[OK] Wrote QC summary: {summary}")

def export_aa(overwrite=False, threads=1):
    for mag in mag_names():
        cdb = MYDATA / f"{mag}-DB.db"
        aa = AA_DIR / f"{mag}.aa.fa"
        if aa.exists() and not overwrite:
            print(f"[SKIP] AA exists: {aa.name}")
        else:
            sh(f'anvi-get-sequences-for-gene-calls -c "{cdb}" --get-aa-sequences -o "{aa}" --overwrite-output-destinations -v')

def hmmer_extract_hits(evalue="1e-12"):
    # ensure HMM is pressed (hmmscan uses .h3*; if not, try hmmpress)
    if not (HMMS / "genes.hmm.h3i").exists():
        if which("hmmpress"):
            sh(f'hmmpress "{PERSO_HMM}"')
        else:
            print("[WARN] hmmpress not available; ensure genes.hmm.h3* exist.")

    all_hits = HITS_DIR / "all_viractin_hits.fa"
    # truncate combined file
    open(all_hits, "w").close()

    for mag in mag_names():
        aa = AA_DIR / f"{mag}.aa.fa"
        tbl = HMMER_DIR / f"{mag}.tbl"
        ids = HMMER_DIR / f"{mag}.ids"
        out = HITS_DIR / f"{mag}_viractin.fa"
        if not aa.exists():
            print(f"[SKIP] No AA for {mag}")
            continue
        # run hmmscan
        sh(f'hmmscan --cpu 1 --tblout "{tbl}" -E {evalue} "{PERSO_HMM}" "{aa}" > /dev/null')
        # collect sequence IDs from column 4 (query)
        # exclude comments
        sh(f"awk '!/^#/{{print $4}}' '{tbl}' | sort -u > '{ids}'", check=True)
        # extract
        seq_count = sh(f'seqkit grep -f "{ids}" "{aa}" > "{out}" && grep -c "^>" "{out}" || echo 0')
        # append to combined if non-empty
        if os.path.getsize(out) > 0:
            sh(f'cat "{out}" >> "{all_hits}"')

    print(f"[OK] Combined hits: {all_hits}  (seqs: {sh(f"grep -c ^> {all_hits} || echo 0", check=False).strip()})")

def build_alignment_and_tree(threads=1):
    raw = HITS_DIR / "all_viractin_hits.fa"
    if not raw.exists() or os.path.getsize(raw) == 0:
        print("[ERROR] No viractin hits found. Aborting MSA/Tree.")
        sys.exit(2)

    # Optionally add references
    actin_raw = RESULTS / "actin_plus_refs.raw.fa"
    if REF_FASTA.exists():
        sh(f'cat "{REF_FASTA}" "{raw}" > "{actin_raw}"')
    else:
        shutil.copy(raw, actin_raw)

    # Clean FASTA (remove '*', trailing '-' and normalize)
    cleaned = RESULTS / "actin_plus_refs.fa"
    awk = (
        r"""awk 'BEGIN{RS=">"; ORS=""} NR>1{hdr=$1; sub(/[ \t\r\n]+$/, "", hdr); """
        r"""seq=$0; sub(/^.*\n/,"",seq); gsub(/\n/,"",seq); gsub(/\*/,"",seq); gsub(/-+$/,"",seq); """
        r"""print ">" hdr "\n" seq "\n"}' """
        f'"{actin_raw}" > "{cleaned}"'
    )
    sh(awk)

    # choose aligner
    aligner = "muscle" if which("muscle") else "mafft"
    aln = RESULTS / ("actin_plus_refs.muscle.fa" if aligner == "muscle" else "actin_plus_refs.mafft.fa")
    if aligner == "muscle":
        sh(f'muscle -in "{cleaned}" -out "{aln}"')
    else:
        sh(f'mafft --thread {threads} --auto "{cleaned}" > "{aln}"')

    # trim with trimAl -> PHYLIP
    trimmed = RESULTS / "actin_plus_refs.trimal.phy"
    sh(f'trimal -automated1 -in "{aln}" -out "{trimmed}" -phylip')

    # IQ-TREE
    iqlog = LOGS / "05_iqtree.log"
    sh(f'iqtree -s "{trimmed}" -m LG+R7 -alrt 1000 -bb 1000 -nt {threads}', log=iqlog)

    treefile = RESULTS / "actin_plus_refs.trimal.phy.treefile"
    iqreport = RESULTS / "actin_plus_refs.trimal.phy.iqtree"
    print(f"[OK] Tree: {treefile}")
    print(f"[OK] IQ-TREE report: {iqreport}")

def main():
    ap = argparse.ArgumentParser(description="Run Viractin TD/TP pipeline end-to-end.")
    ap.add_argument("--threads", type=int, default=1, help="Threads for IQ-TREE/MAFFT.")
    ap.add_argument("--overwrite", action="store_true", help="Rebuild even if files exist.")
    args = ap.parse_args()

    ensure_dirs()
    aligner = check_deps()

    print("\n=== STEP 1: Build contigs DBs (if needed) ===")
    build_contigs_dbs(overwrite=args.overwrite)

    print("\n=== STEP 2: Run PersoHMM on each contigs DB ===")
    run_perso_hmm_on_each_db(overwrite=args.overwrite, threads=args.threads)

    print("\n=== STEP 3: QC completeness & redundancy summary ===")
    qc_completeness(overwrite=args.overwrite)

    print("\n=== STEP 4: Export AA sequences per MAG ===")
    export_aa(overwrite=args.overwrite, threads=args.threads)

    print("\n=== STEP 5: HMMER scan & extract Viractins ===")
    hmmer_extract_hits(evalue="1e-12")

    print("\n=== STEP 6: Alignment, trimming, and ML tree ===")
    build_alignment_and_tree(threads=args.threads)

    print("\nAll done ✅ — see results/ and logs/ for outputs.")
    print("Key files:")
    print("  - results/hmm_hits/all_viractin_hits.fa")
    print("  - results/actin_plus_refs.trimal.phy (trimmed alignment)")
    print("  - results/actin_plus_refs.trimal.phy.treefile (ML tree)")
    print("  - results/actin_plus_refs.trimal.phy.iqtree (model, supports)")

if __name__ == "__main__":
    main()

