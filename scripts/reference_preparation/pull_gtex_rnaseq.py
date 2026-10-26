#!/usr/bin/env python3
"""
Fresh pull of GTEx normal-tissue RNA-seq, labeled by tissue -- the normal-
sample counterpart to pull_gdc_rnaseq.py's TCGA/CPTAC tumor pull.

Unlike GDC, GTEx already publishes one pre-split gene expression file PER
TISSUE directly from the source (no per-sample file assembly needed):
  https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/
    tpms-by-tissue/gene_tpm_adult_gtex_v11_<tissue>.gct.gz
    counts-by-tissue/gene_reads_adult_gtex_v11_<tissue>.gct.gz
54 tissues, release GTEx v11 (2026-05-19 build). Sample -> subject/tissue
labels come from the SampleAttributesDS.txt annotation file (SAMPID, SMTS).

Output layout matches pull_gdc_rnaseq.py exactly, so both feed the same
build_pancancer_matrix.py unchanged:
  output/GTEx/<tissue>/<tissue>.tpm_matrix.tsv.gz
  output/GTEx/<tissue>/<tissue>.counts_matrix.tsv.gz
  output/GTEx/<tissue>/<tissue>.sample_manifest.tsv
pushed to blob at cureseek/public-data/gtex-rnaseq/<tissue>/.

Usage:
    source .venv/bin/activate
    python3 pull_gtex_rnaseq.py [--only bladder,liver] [--dry-run]
"""
import argparse
import io
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

GCS_ROOT = "https://storage.googleapis.com/adult-gtex"
RELEASE = "v11"
GCT_TAG = "adult_gtex_v11"
OUT_ROOT = Path(__file__).resolve().parent / "output"
BLOB_ROOT = "cureseekblob:cureseek-blob/cureseek/public-data"
SAMPLE_ATTRS_URL = f"{GCS_ROOT}/annotations/{RELEASE}/metadata-files/GTEx_Analysis_{RELEASE}_Annotations_SampleAttributesDS.txt"


def log(msg):
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}", flush=True)


def list_tissues():
    r = requests.get("https://storage.googleapis.com/storage/v1/b/adult-gtex/o", params={
        "prefix": f"bulk-gex/{RELEASE}/rna-seq/tpms-by-tissue/gene_tpm_{GCT_TAG}_",
        "fields": "items(name)",
    }, timeout=60)
    r.raise_for_status()
    prefix = f"bulk-gex/{RELEASE}/rna-seq/tpms-by-tissue/gene_tpm_{GCT_TAG}_"
    return sorted(item["name"][len(prefix):-len(".gct.gz")] for item in r.json()["items"])


def load_sample_attrs():
    r = requests.get(SAMPLE_ATTRS_URL, timeout=120)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t", usecols=["SAMPID", "SMTS", "SMTSD"])
    return df.set_index("SAMPID")


def download_gct(url):
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=300)
            r.raise_for_status()
            return pd.read_csv(io.BytesIO(r.content), sep="\t", skiprows=2,
                                compression="gzip", index_col=0).drop(columns=["Description"])
        except Exception as e:
            log(f"  download attempt {attempt + 1} failed ({url}): {e}")
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"failed to download {url} after 3 attempts")


def build_matrix_for_tissue(tissue, attrs):
    out_dir = OUT_ROOT / "GTEx" / tissue
    tpm_path = out_dir / f"{tissue}.tpm_matrix.tsv.gz"
    if tpm_path.exists():
        log(f"[GTEx/{tissue}] already built, skipping")
        return

    tpm_url = f"{GCS_ROOT}/bulk-gex/{RELEASE}/rna-seq/tpms-by-tissue/gene_tpm_{GCT_TAG}_{tissue}.gct.gz"
    reads_url = f"{GCS_ROOT}/bulk-gex/{RELEASE}/rna-seq/counts-by-tissue/gene_reads_{GCT_TAG}_{tissue}.gct.gz"
    log(f"[GTEx/{tissue}] downloading tpm + reads")
    tpm = download_gct(tpm_url)
    counts = download_gct(reads_url)

    manifest_rows = []
    for sampid in tpm.columns:
        subject_id = "-".join(sampid.split("-")[:2])
        smts = attrs["SMTS"].get(sampid)
        manifest_rows.append({
            "file_id": f"gene_tpm_{GCT_TAG}_{tissue}.gct.gz",
            "case_id": subject_id,
            "case_submitter_id": subject_id,
            "sample_barcode": sampid,
            "sample_type": "Normal Tissue",
            "project_id": "GTEx",
            "primary_site": smts if pd.notna(smts) else tissue,
        })
    manifest = pd.DataFrame(manifest_rows)

    out_dir.mkdir(parents=True, exist_ok=True)
    tpm.to_csv(tpm_path, sep="\t", compression="gzip")
    counts.to_csv(out_dir / f"{tissue}.counts_matrix.tsv.gz", sep="\t", compression="gzip")
    manifest.to_csv(out_dir / f"{tissue}.sample_manifest.tsv", sep="\t", index=False)
    log(f"[GTEx/{tissue}] wrote {tpm.shape[1]} samples x {tpm.shape[0]} genes")

    blob_dest = f"{BLOB_ROOT}/gtex-rnaseq/{tissue}/"
    subprocess.run(["rclone", "copy", str(out_dir), blob_dest, "--stats-one-line"], check=True)
    log(f"[GTEx/{tissue}] uploaded -> {blob_dest}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="comma-separated tissue list, for testing a subset")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    tissues = list_tissues()
    if args.only:
        only = set(args.only.split(","))
        tissues = [t for t in tissues if t in only]
    log(f"{len(tissues)} tissues queued: {tissues}")
    if args.dry_run:
        return

    log("loading sample attributes (subject/tissue metadata)")
    attrs = load_sample_attrs()

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(build_matrix_for_tissue, t, attrs): t for t in tissues}
        for future in as_completed(futures):
            tissue = futures[future]
            try:
                future.result()
            except Exception as e:
                log(f"[GTEx/{tissue}] FAILED: {e}")

    log("done")


if __name__ == "__main__":
    main()
