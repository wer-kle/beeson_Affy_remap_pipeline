#!/usr/bin/env python3
"""
01_make_flanks_Affy_MNEc670.py

Generate reference-context flanking sequences for each Affymetrix MNEc670 SNP,
analogous in purpose to Beeson et al.'s 01_make_probeseq.py.

For each SNP in the NetAffx annotation:
  - Take the 'Flank' field: e.g.  TGTT...ACATA[A/C]GGTA...TATA
  - Take 'Ref Allele' (EquCab3 base at that coordinate).
  - Construct a single reference-strand context sequence:
        flank_ref = left + RefAllele + right
  - Write to FASTA with ID based on Probe Set ID (AX-...).

We restrict to:
  - biallelic SNPs with Ref Allele and Alt Allele in {A, C, G, T}
  - flanks that contain exactly one [X/Y] bracket.

Output FASTA will be used as the "probe" input for Beeson-style BLAST
remapping to EquCab3 (and beyond).
"""

import os
import sys
import re
import pandas as pd


def make_affy_flanks(annot_path, out_fasta_path):
    """Read NetAffx CSV and write reference-context flanks as FASTA."""
    sys.stderr.write(f"[INFO] Reading Affy annotation: {annot_path}\n")
    df = pd.read_csv(annot_path, comment="#")

    required_cols = ["Probe Set ID", "Flank", "Ref Allele", "Alt Allele"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        sys.stderr.write(f"[ERR] Missing required columns: {', '.join(missing)}\n")
        sys.exit(1)

    n_rows = len(df)
    sys.stderr.write(f"[INFO] Annotation rows: {n_rows}\n")

    out_dir = os.path.dirname(out_fasta_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Regex to parse Flank: left [A/B] right
    # Require exactly one [X/Y] with A/C/G/T alleles.
    flank_pattern = re.compile(r"^(.*)\[([ACGT])\/([ACGT])\](.*)$", re.IGNORECASE)

    n_written = 0
    n_skipped_nonnuc = 0
    n_skipped_flank = 0
    n_skipped_ref_mismatch = 0

    with open(out_fasta_path, "w") as out_f:
        for _, row in df.iterrows():
            probe_id = str(row["Probe Set ID"]).strip()
            flank_raw = str(row["Flank"]).strip()
            ref = str(row["Ref Allele"]).strip().upper()
            alt = str(row["Alt Allele"]).strip().upper()

            # Require standard nucleotide ref/alt
            if ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"}:
                n_skipped_nonnuc += 1
                continue

            m = flank_pattern.match(flank_raw)
            if not m:
                n_skipped_flank += 1
                continue

            left = m.group(1).upper()
            a1 = m.group(2).upper()
            a2 = m.group(3).upper()
            right = m.group(4).upper()

            # Optional sanity: Ref Allele should match one of the bracket alleles.
            if ref not in {a1, a2}:
                n_skipped_ref_mismatch += 1
                continue

            flank_ref = left + ref + right

            # FASTA ID: use Probe Set ID (AX-...), keep it simple
            fasta_id = probe_id

            out_f.write(f">{fasta_id}\n{flank_ref}\n")
            n_written += 1

    sys.stderr.write(f"[INFO] Wrote flanks for {n_written} SNPs to {out_fasta_path}\n")
    sys.stderr.write(f"[INFO] Skipped {n_skipped_nonnuc} rows (non-ACGT ref/alt)\n")
    sys.stderr.write(f"[INFO] Skipped {n_skipped_flank} rows (Flank not parseable)\n")
    sys.stderr.write(f"[INFO] Skipped {n_skipped_ref_mismatch} rows (Ref not in [A/B])\n")


def main():
    project_root = "/Volumes/Public/SNP_array_data/horse_data"

    annot_rel = "ref/manifest_files/Affy_annotation_file.csv"
    out_rel = "metadata/remap_beeson/Affy_MNEc670_flanks_ref.fa"

    annot_path = os.path.join(project_root, annot_rel)
    out_path = os.path.join(project_root, out_rel)

    if not os.path.exists(annot_path):
        sys.stderr.write(f"[ERR] Annotation file not found: {annot_path}\n")
        sys.exit(1)

    make_affy_flanks(annot_path, out_path)


if __name__ == "__main__":
    main()
