#!/usr/bin/env python3
"""
Parquet cache for the wide gdc_tpm_matrix.tsv.gz / gtex_tpm_matrix.tsv.gz
built by build_pancancer_matrix.py.

Two problems with reading the .tsv.gz directly, on every analysis run:
  1. Slow: a full-width gzip CSV parse just to read a few hundred columns
     out of 14,592 (GDC) / 18,637 (GTEx).
  2. Fragile: pandas' `usecols` validates against the file's literal header
     names, and those aren't even consistent between the two files --
     GDC's index column is "gene_id" with GENCODE version suffixes
     (ENSG00000000003.15); GTEx's is "Name" (inherited straight from GTEx's
     own GCT format, never normalized) with its own version suffixes.

This module reads each .tsv.gz ONCE, canonicalizes the index to a
version-stripped "gene_id" column, and writes a .parquet cache next to it.
Every subsequent read (from here or any other script) can then select just
the needed sample columns via pyarrow's schema-based column pruning --
fast, and immune to the header-naming inconsistency above since callers
only ever see "gene_id".

Usage as a library:
    from matrix_cache import ensure_parquet, read_matrix_columns
    parquet_path = ensure_parquet(GDC_MATRIX_TSV_GZ)
    df = read_matrix_columns(parquet_path, ["TCGA-GBM|TCGA-06-0125|TP", ...])
"""
from pathlib import Path

import pandas as pd


def _strip_version(value):
    return str(value).split(".")[0]


def ensure_parquet(tsv_gz_path):
    """Build <name>.parquet next to <name>.tsv.gz if missing or stale.

    Reads the index positionally (index_col=0) rather than by name -- GDC's
    own index column is called "gene_id", GTEx's is "Name" (inherited from
    GTEx's own GCT format, never normalized). Both get canonicalized to a
    version-stripped "gene_id" here, so callers never need to know which.
    """
    tsv_gz_path = Path(tsv_gz_path)
    parquet_path = tsv_gz_path.with_name(tsv_gz_path.name.replace(".tsv.gz", ".parquet"))

    if parquet_path.exists() and not tsv_gz_path.exists():
        # e.g. an Azure Batch node that only ever gets the prebuilt .parquet
        # staged, never the source .tsv.gz -- nothing to compare staleness
        # against, so just trust the cache.
        return parquet_path
    if parquet_path.exists() and parquet_path.stat().st_mtime >= tsv_gz_path.stat().st_mtime:
        return parquet_path

    print(f"[matrix_cache] building {parquet_path.name} from {tsv_gz_path.name} (one-time)...", flush=True)
    df = pd.read_csv(tsv_gz_path, sep="\t", index_col=0)
    df.index = df.index.map(_strip_version)
    df.index.name = "gene_id"
    if df.index.duplicated().any():
        # GENCODE PAR_Y-style duplicates after version stripping are rare but
        # possible -- keep the first occurrence rather than silently summing
        # or erroring, and say so.
        n_dupes = df.index.duplicated().sum()
        print(f"[matrix_cache] WARNING: {n_dupes} duplicate gene_id after version-stripping, keeping first", flush=True)
        df = df[~df.index.duplicated(keep="first")]
    df = df.reset_index()

    tmp_path = parquet_path.with_suffix(".parquet.tmp")
    df.to_parquet(tmp_path, engine="pyarrow", compression="snappy", index=False)
    tmp_path.rename(parquet_path)
    print(f"[matrix_cache] wrote {parquet_path} ({parquet_path.stat().st_size / 1e9:.2f} GB)", flush=True)
    return parquet_path


def read_matrix_columns(parquet_path, sample_columns):
    """Read only gene_id + the requested sample columns, indexed by gene_id."""
    df = pd.read_parquet(parquet_path, columns=["gene_id", *sample_columns], engine="pyarrow")
    return df.set_index("gene_id")
