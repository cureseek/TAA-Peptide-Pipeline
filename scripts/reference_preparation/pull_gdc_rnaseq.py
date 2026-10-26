#!/usr/bin/env python3
"""
Fresh pull of TCGA + CPTAC RNA-seq expression, from GDC, as gene x sample
matrices labeled by cancer type -- not raw per-sample files.

Both cohorts are harmonized by GDC through the same pipeline (STAR aligned,
GENCODE-annotated gene counts), so one code path covers both:
  TCGA: 33 projects, each project_id IS the cancer type (TCGA-BRCA, ...).
  CPTAC: all cases sit under two umbrella projects (CPTAC-2, CPTAC-3) that
         each span several organs, so cases.primary_site is used as the
         cancer-type label instead. Head & neck sub-sites (larynx, tongue,
         mouth, pharynx, ...) are merged into one Head_and_Neck bucket --
         GDC records those as ~10 separate ICD-O sites, but they're all the
         same CPTAC-HNSCC cohort. Everything else with no clear bucket (small
         n, ambiguous site) falls into "Other" rather than getting its own
         near-empty file.

Per cancer type, three files land in `output/<source>/<label>/` and get
pushed to blob at cureseek/public-data/<source-lower>-rnaseq/<label>/:
  <label>.tpm_matrix.tsv.gz     genes (rows) x samples (cols), tpm_unstranded
  <label>.counts_matrix.tsv.gz  genes (rows) x samples (cols), raw STAR counts
  <label>.sample_manifest.tsv   sample_barcode, case_id, sample_type, file_id

Idempotent: a label already present in output/ (tpm_matrix.tsv.gz exists) is
skipped, so a killed/interrupted run can just be restarted.

Usage:
    source .venv/bin/activate
    python3 pull_gdc_rnaseq.py [--only TCGA-BRCA,Kidney] [--dry-run]
"""
import argparse
import io
import json
import shutil
import subprocess
import tarfile
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

GDC_API = "https://api.gdc.cancer.gov"
OUT_ROOT = Path(__file__).resolve().parent / "output"
BLOB_ROOT = "cureseekblob:cureseek-blob/cureseek/public-data"
BATCH_SIZE = 200  # file ids per bulk /data POST
GENE_EXPR_FILTERS = [
    {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
    {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["RNA-Seq"]}},
    {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": ["STAR - Counts"]}},
    {"op": "in", "content": {"field": "files.access", "value": ["open"]}},
]

HEAD_AND_NECK_SITES = {
    "Larynx", "Other and unspecified parts of tongue", "Floor of mouth",
    "Other and ill-defined sites in lip, oral cavity and pharynx",
    "Oropharynx", "Gum", "Lip", "Other and unspecified parts of mouth",
    "Base of tongue", "Tonsil",
}
CPTAC_SITE_LABELS = {
    "Bronchus and lung": "Lung",
    "Kidney": "Kidney",
    "Uterus, NOS": "Uterus",
    "Pancreas": "Pancreas",
    "Skin": "Skin",
    "Brain": "Brain",
    "Stomach": "Stomach",
    "Hematopoietic and reticuloendothelial systems": "Hematopoietic",
    "Breast": "Breast",
    "Colon": "Colon",
    "Ovary": "Ovary",
}


def log(msg):
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}", flush=True)


def gdc_post(endpoint, payload):
    r = requests.post(f"{GDC_API}/{endpoint}", json=payload, timeout=120)
    r.raise_for_status()
    return r


def query_files(project_ids, primary_sites=None):
    """Page through GDC /files, return one row per file with the metadata
    needed to name and group matrix columns."""
    filters = [{"op": "in", "content": {"field": "cases.project.project_id", "value": project_ids}}]
    filters += GENE_EXPR_FILTERS
    if primary_sites:
        filters.append({"op": "in", "content": {"field": "cases.primary_site", "value": primary_sites}})
    payload_filters = {"op": "and", "content": filters}

    fields = ",".join([
        "file_id", "cases.case_id", "cases.submitter_id",
        "cases.samples.submitter_id", "cases.samples.sample_type",
        "cases.project.project_id", "cases.primary_site",
    ])

    rows = []
    frm = 1
    size = 500
    while True:
        r = gdc_post("files", {
            "filters": payload_filters, "fields": fields,
            "size": size, "from": frm, "format": "json",
        })
        data = r.json()["data"]
        hits = data["hits"]
        for h in hits:
            case = (h.get("cases") or [{}])[0]
            samples = case.get("samples") or [{}]
            rows.append({
                "file_id": h["file_id"],
                "case_id": case.get("case_id"),
                "case_submitter_id": case.get("submitter_id"),
                "sample_barcode": samples[0].get("submitter_id") or case.get("submitter_id"),
                "sample_type": samples[0].get("sample_type"),
                "project_id": (case.get("project") or {}).get("project_id"),
                "primary_site": case.get("primary_site"),
            })
        total = data["pagination"]["total"]
        if frm + size > total:
            break
        frm += size
    return pd.DataFrame(rows).drop_duplicates(subset="file_id")


def download_batch(file_ids, dest_dir):
    """POST a batch of file ids to GDC's bulk /data endpoint, extract the
    returned tarball into dest_dir/<file_id>/<name>.tsv. Retries on failure
    since a multi-hour run over ~15k files will hit transient network blips."""
    for attempt in range(3):
        try:
            r = requests.post(f"{GDC_API}/data", json={"ids": file_ids}, timeout=600)
            r.raise_for_status()
            with tarfile.open(fileobj=io.BytesIO(r.content)) as tf:
                tf.extractall(dest_dir)
            return
        except Exception as e:
            log(f"  batch download attempt {attempt + 1} failed: {e}")
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"batch of {len(file_ids)} files failed after 3 attempts")


def parse_star_counts(path):
    """Return (gene_id-indexed tpm Series, gene_id-indexed counts Series)."""
    df = pd.read_csv(path, sep="\t", skiprows=1)
    df = df[df["gene_id"].str.startswith("ENSG", na=False)]
    df = df.set_index("gene_id")
    return df["tpm_unstranded"], df["unstranded"]


def build_matrix_for_label(source, label, manifest):
    out_dir = OUT_ROOT / source / label
    tpm_path = out_dir / f"{label}.tpm_matrix.tsv.gz"
    manifest_path = out_dir / f"{label}.sample_manifest.tsv"
    if tpm_path.exists() and manifest_path.exists():
        log(f"[{source}/{label}] already built, skipping ({len(manifest)} files in manifest)")
        return

    log(f"[{source}/{label}] {len(manifest)} files to fetch")
    out_dir.mkdir(parents=True, exist_ok=True)
    tpm_cols, count_cols = {}, {}

    file_ids = manifest["file_id"].tolist()
    with tempfile.TemporaryDirectory(prefix=f"gdc_{label}_") as td:
        for i in range(0, len(file_ids), BATCH_SIZE):
            batch = file_ids[i:i + BATCH_SIZE]
            log(f"[{source}/{label}] batch {i // BATCH_SIZE + 1}/{-(-len(file_ids) // BATCH_SIZE)} ({len(batch)} files)")
            download_batch(batch, td)
            for fid in batch:
                fdir = Path(td) / fid
                tsvs = list(fdir.glob("*.tsv"))
                if not tsvs:
                    log(f"  WARNING: no tsv extracted for {fid}, skipping")
                    continue
                barcode = manifest.loc[manifest["file_id"] == fid, "sample_barcode"].iloc[0]
                tpm, counts = parse_star_counts(tsvs[0])
                tpm_cols[barcode] = tpm
                count_cols[barcode] = counts
            shutil.rmtree(td, ignore_errors=True)
            Path(td).mkdir(exist_ok=True)

    if not tpm_cols:
        log(f"[{source}/{label}] WARNING: nothing parsed, no output written")
        return

    tpm_matrix = pd.concat(tpm_cols, axis=1)
    count_matrix = pd.concat(count_cols, axis=1)
    tpm_matrix.to_csv(out_dir / f"{label}.tpm_matrix.tsv.gz", sep="\t", compression="gzip")
    count_matrix.to_csv(out_dir / f"{label}.counts_matrix.tsv.gz", sep="\t", compression="gzip")
    manifest.to_csv(out_dir / f"{label}.sample_manifest.tsv", sep="\t", index=False)
    log(f"[{source}/{label}] wrote {tpm_matrix.shape[1]} samples x {tpm_matrix.shape[0]} genes")

    blob_dest = f"{BLOB_ROOT}/{source.lower()}-rnaseq/{label}/"
    subprocess.run(["rclone", "copy", str(out_dir), blob_dest, "--stats-one-line"], check=True)
    log(f"[{source}/{label}] uploaded -> {blob_dest}")


def tcga_labels():
    r = requests.get(f"{GDC_API}/projects", params={
        "filters": json.dumps({"op": "in", "content": {"field": "program.name", "value": ["TCGA"]}}),
        "fields": "project_id", "size": "100", "format": "json",
    }, timeout=60)
    r.raise_for_status()
    project_ids = sorted(h["project_id"] for h in r.json()["data"]["hits"])
    return [(pid, [pid], None) for pid in project_ids]


def cptac_labels():
    manifest = query_files(["CPTAC-2", "CPTAC-3"])
    labels = {}
    for site in manifest["primary_site"].dropna().unique():
        label = "Head_and_Neck" if site in HEAD_AND_NECK_SITES else CPTAC_SITE_LABELS.get(site, "Other")
        labels.setdefault(label, []).append(site)
    return [(label, ["CPTAC-2", "CPTAC-3"], sites) for label, sites in sorted(labels.items())]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="comma-separated label list, for testing a subset")
    ap.add_argument("--dry-run", action="store_true", help="print what would be pulled, don't download")
    ap.add_argument("--workers", type=int, default=4, help="cancer-type labels to fetch concurrently")
    args = ap.parse_args()
    only = set(args.only.split(",")) if args.only else None

    jobs = [("TCGA", *t) for t in tcga_labels()] + [("CPTAC", *c) for c in cptac_labels()]
    if only:
        jobs = [j for j in jobs if j[1] in only]

    log(f"{len(jobs)} cancer-type labels queued: {[j[1] for j in jobs]}")
    labeled_manifests = []
    for source, label, project_ids, primary_sites in jobs:
        manifest = query_files(project_ids, primary_sites)
        if args.dry_run:
            log(f"[{source}/{label}] would fetch {len(manifest)} files")
            continue
        labeled_manifests.append((source, label, manifest))

    if args.dry_run:
        return

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(build_matrix_for_label, source, label, manifest): label
            for source, label, manifest in labeled_manifests
        }
        for future in as_completed(futures):
            label = futures[future]
            try:
                future.result()
            except Exception as e:
                log(f"[{label}] FAILED: {e}")

    log("done")


if __name__ == "__main__":
    main()
