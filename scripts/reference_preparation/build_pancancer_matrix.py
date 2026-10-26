#!/usr/bin/env python3
"""
Combine per-cancer-type/tissue TPM matrices from pull_gdc_rnaseq.py and
pull_gtex_rnaseq.py's output/ into two separate gene x sample matrices --
one for GDC (TCGA+CPTAC tumor/adjacent-normal) and one for GTEx (normal
donor tissue). Kept separate rather than joined into one pan-cancer matrix:
GDC's `tpm_unstranded` (pull_gdc_rnaseq.py's parse_star_counts()) and GTEx's
own RNA-SeQC-pipeline TPM are two different quantification pipelines with no
cross-source renormalization, and GDC's 60,660-gene set vs GTEx's 74,628-gene
set only partially overlap -- outer-joining them into one frame meant ~56%
of the combined gene axis was NaN for one side or the other, which is both
scientifically iffy (mixing pipelines under one column) and, in practice,
made the final gzip write take hours for a >4GB file. Splitting means each
matrix's gene index is uniform across all its own labels, so no join/NaN
inflation and a much smaller, much faster write.

GTEx's two cultured-cell-line labels (cells_cultured_fibroblasts,
cells_ebv-transformed_lymphocytes) are excluded -- they're genotyping/eQTL
reference cell lines, not primary tissue, and don't belong in a "normal
tissue" baseline.

Column label format: <cancer_type_or_tissue>|<patient_or_donor_id>|<sample_type_code>
  e.g. TCGA-BRCA|TCGA-A1-A0SD|TP   (TP = Primary Tumor, NT = Solid Tissue
  Normal, see SAMPLE_TYPE_CODES below for the rest; GTEx samples are all GTEX)

Only combines whatever labels have already finished (checks for
<label>.tpm_matrix.tsv.gz under output/<source>/) -- safe to run before
pull_gdc_rnaseq.py/pull_gtex_rnaseq.py are fully done, and to re-run again
once more labels land.

Output: output/gdc_tpm_matrix.tsv.gz + output/gdc_manifest.tsv
        output/gtex_tpm_matrix.tsv.gz + output/gtex_manifest.tsv
Pushed to blob at cureseek/public-data/<name>.

Usage:
    source .venv/bin/activate
    python3 build_pancancer_matrix.py
"""
import argparse
import subprocess
import time
from pathlib import Path

import pandas as pd

OUT_ROOT = Path(__file__).resolve().parent / "output"
BLOB_ROOT = "cureseekblob:cureseek-blob/cureseek/public-data"

# Faster than pandas' to_csv(compression="gzip") default (gzip level 9) --
# level 9 was the main reason the single combined matrix took hours to
# write; level 3 is a large speedup for a small size cost.
GZIP_LEVEL = 3

SAMPLE_TYPE_CODES = {
    "Primary Tumor": "TP",
    "Solid Tissue Normal": "NT",
    "Recurrent Tumor": "TR",
    "Metastatic": "TM",
    "Primary Blood Derived Cancer - Bone Marrow": "TBM",
    "Primary Blood Derived Cancer - Peripheral Blood": "TPB",
    "Normal Tissue": "GTEX",  # GTEx donor samples, not tumor-adjacent
}

GTEX_EXCLUDED_LABELS = {"cells_cultured_fibroblasts", "cells_ebv-transformed_lymphocytes"}


def log(msg):
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}", flush=True)


def find_finished_labels(source_dirs):
    paths = []
    for source_dir in source_dirs:
        found = sorted(OUT_ROOT.glob(f"{source_dir}/*/*.tpm_matrix.tsv.gz"))
        if source_dir == "GTEx":
            found = [p for p in found if p.parent.name not in GTEX_EXCLUDED_LABELS]
        paths.extend(found)
    return paths


def column_label(cancer_type, case_submitter_id, sample_type):
    code = SAMPLE_TYPE_CODES.get(sample_type, "NA")
    return f"{cancer_type}|{case_submitter_id}|{code}"


def build_source_matrix(source_dirs, out_name):
    tag = "+".join(source_dirs)
    matrix_paths = find_finished_labels(source_dirs)
    if not matrix_paths:
        log(f"[{tag}] no finished labels found under output/ -- nothing to combine")
        return
    log(f"[{tag}] combining {len(matrix_paths)} finished label matrices")

    per_label_frames = []
    manifest_rows = []
    for tpm_path in matrix_paths:
        label_dir = tpm_path.parent
        label = label_dir.name
        source = label_dir.parent.name
        manifest = pd.read_csv(label_dir / f"{label}.sample_manifest.tsv", sep="\t")

        tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
        rename = {}
        for _, row in manifest.iterrows():
            barcode = row["sample_barcode"]
            if barcode not in tpm.columns:
                continue
            new_col = column_label(label, row["case_submitter_id"], row["sample_type"])
            rename[barcode] = new_col
            manifest_rows.append({
                "column_label": new_col,
                "source": source,
                "cancer_type": label,
                "patient_id": row["case_submitter_id"],
                "sample_barcode": barcode,
                "sample_type": row["sample_type"],
                "file_id": row["file_id"],
            })
        tpm = tpm.rename(columns=rename)
        per_label_frames.append(tpm)
        log(f"  {source}/{label}: {tpm.shape[1]} samples")

    log(f"[{tag}] joining on gene_id (labels share one gene set -- no NaN inflation expected)")
    combined = pd.concat(per_label_frames, axis=1, join="outer")
    n_incomplete_genes = combined.isna().any(axis=1).sum()
    if n_incomplete_genes:
        log(f"  WARNING: {n_incomplete_genes} genes are missing from at least one label's file "
            f"within {tag} (expected to be 0) -- left as NaN, not dropped")

    manifest_df = pd.DataFrame(manifest_rows)
    dupes = manifest_df["column_label"].duplicated()
    if dupes.any():
        log(f"  WARNING: {dupes.sum()} duplicate column labels (same patient/cancer_type/sample_type "
            f"appearing twice, e.g. two aliquots of the same sample) -- both columns are kept side by "
            f"side under the identical label, not merged or dropped; disambiguate via the manifest's "
            f"file_id/sample_barcode if you need to tell them apart")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    matrix_out = OUT_ROOT / f"{out_name}_tpm_matrix.tsv.gz"
    manifest_out = OUT_ROOT / f"{out_name}_manifest.tsv"
    combined.to_csv(matrix_out, sep="\t", compression={"method": "gzip", "compresslevel": GZIP_LEVEL})
    manifest_df.to_csv(manifest_out, sep="\t", index=False)
    log(f"[{tag}] wrote {combined.shape[1]} samples x {combined.shape[0]} genes -> {matrix_out}")

    subprocess.run(["rclone", "copyto", str(matrix_out), f"{BLOB_ROOT}/{matrix_out.name}"], check=True)
    subprocess.run(["rclone", "copyto", str(manifest_out), f"{BLOB_ROOT}/{manifest_out.name}"], check=True)
    log(f"[{tag}] uploaded -> {BLOB_ROOT}/")


GROUPS = {
    "gdc": (["TCGA", "CPTAC"], "gdc"),
    "gtex": (["GTEx"], "gtex"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", choices=sorted(GROUPS), help="build only this group (default: both, sequentially)")
    args = ap.parse_args()

    groups = [GROUPS[args.group]] if args.group else list(GROUPS.values())
    for source_dirs, out_name in groups:
        build_source_matrix(source_dirs, out_name)


if __name__ == "__main__":
    main()
