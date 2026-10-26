#!/usr/bin/env python3
"""
TAA reference-protein peptide binding screen.

Inputs:
  /data4/NewOutput/<rna_patient>/rnaseq-results/TAA.xml
  /data3/PipelineOutput/<hla_patient>/HLAclassification/final_ClassI.HLAgenotype4digits
  /data3/PipelineOutput/<hla_patient>/HLAclassification/final_ClassII.HLAgenotype4digits

Outputs:
  /data4/NewOutput/<rna_patient>/rnaseq-results/TAA/Targets_Gene_Sequence.txt
  /data4/NewOutput/<rna_patient>/rnaseq-results/TAA/TAA_Peptide.xml
  /data4/NewOutput/<rna_patient>/rnaseq-results/TAA/protein_manifest.json
  /data4/NewOutput/<rna_patient>/rnaseq-results/TAA/peptide_prediction/

Class I:
  Natural 9-mers
  MHCflurry affinity prediction
  IC50 < 500 nM
  Rank by:
    1. Number of binding HLA alleles (descending)
    2. Mean IC50 (ascending)
    3. Minimum IC50 (ascending)
    4. Peptide sequence
  Top 10 per gene

Class II:
  Natural 15-mers
  MixMHC2pred --no_context
  Rank by %Rank_best ascending
  Top 10 per gene

No anchor engineering.
No theoretical gene explanations.

Note:
Predicted binding does not demonstrate endogenous antigen processing,
MHC presentation, immunogenicity, or clinical suitability.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlencode, quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET

import pandas as pd


# ============================================================
# Fixed configuration
# ============================================================

DATA4_ROOT = Path('/data4/NewOutput')
DATA3_ROOT = Path('/data3/PipelineOutput')

MIX_EXE = Path(
    '/data3/tools/MixMHC2pred-2.0/MixMHC2pred_unix'
)

MIX_PWM = Path(
    '/data3/tools/MixMHC2pred-2.0/PWMdef'
)

MIX_CONDA_ENV = 'mhcii'

UNIPROT = 'https://rest.uniprot.org/uniprotkb'

AA = set('ACDEFGHIKLMNPQRSTVWY')

I_LENGTH = 9
II_LENGTH = 15

I_CUTOFF = 500.0

TOP_N = 10


# ============================================================
# General helpers
# ============================================================

def patient_folder(value):
    name = value.strip()

    if not re.fullmatch(
        r'[A-Za-z0-9][A-Za-z0-9_-]*',
        name
    ):
        raise ValueError(
            'Invalid patient folder name: ' + repr(name)
        )

    return name


def required(path):
    if not path.is_file():
        raise FileNotFoundError(str(path))

    return path


def atomic_text(path, value):
    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    tmp = path.with_name(
        path.name + '.tmp'
    )

    tmp.write_text(
        value,
        encoding='utf-8'
    )

    tmp.replace(path)


# ============================================================
# Read TAA_analysis_01.py XML
# ============================================================

def read_taa(xml_path):

    root = ET.parse(
        required(xml_path)
    ).getroot()

    if root.tag != 'TAAOverexpression':
        raise ValueError(
            f'Unexpected TAA.xml root: {root.tag}'
        )

    result = []
    seen = set()

    for node in root.findall('./Gene'):

        if (
            node.attrib.get(
                'overexpressed', ''
            ).lower() != 'true'
        ):
            raise ValueError(
                'TAA.xml contains a Gene not marked overexpressed'
            )

        symbol = node.attrib.get(
            'name', ''
        ).strip().upper()

        gid = node.attrib.get(
            'ensembl_id', ''
        ).split('.')[0].strip()

        if (
            not symbol
            or not re.fullmatch(
                r'ENSG\d{11}',
                gid
            )
        ):
            raise ValueError(
                f'Invalid Gene symbol/Ensembl ID in XML: '
                f'{node.attrib}'
            )

        if symbol in seen:
            raise ValueError(
                f'Duplicate gene symbol in TAA.xml: {symbol}'
            )

        seen.add(symbol)

        result.append({
            'gene': symbol,
            'ensembl': gid,
            'tpm': node.attrib.get('TPM', '')
        })

    return result


# ============================================================
# HLA Class I
# ============================================================

def parse_hla_i(path):

    text = required(path).read_text(
        encoding='utf-8'
    )

    pairs = re.findall(
        r'(?:HLA-)?([ABC])\*(\d{2}:\d{2})',
        text,
        re.I
    )

    result = list(
        dict.fromkeys(
            f'HLA-{a.upper()}*{code}'
            for a, code in pairs
        )
    )

    if not result:
        raise ValueError(
            'No Class I HLA alleles parsed from '
            + str(path)
        )

    return result


# ============================================================
# HLA Class II
# ============================================================

def parse_hla_ii(path):

    text = required(path).read_text(
        encoding='utf-8'
    )

    pattern = (
        r'(?:HLA-)?'
        r'(DPA1|DPB1|DQA1|DQB1|DRB1|DRB3|DRB4|DRB5)'
        r'\*(\d{2,3}:\d{2,3})'
    )

    loci = defaultdict(list)

    for locus, code in re.findall(
        pattern,
        text,
        re.I
    ):

        allele = (
            f'HLA-{locus.upper()}*{code}'
        )

        if allele not in loci[locus.upper()]:
            loci[locus.upper()].append(allele)

    if not loci:
        raise ValueError(
            'No Class II HLA alleles parsed from '
            + str(path)
        )

    return dict(loci)


def token(allele):

    return (
        allele
        .replace('HLA-', '')
        .replace('*', '_')
        .replace(':', '_')
    )


def norm_token(value):

    s = (
        value.upper()
        .replace('HLA-', '')
        .replace('*', '_')
        .replace(':', '_')
        .replace('/', '__')
        .replace('-', '_')
    )

    s = re.sub(
        r'_+',
        '_',
        s
    )

    # Restore alpha/beta delimiter.
    s = re.sub(
        r'(DPA1_\d+_\d+)_(DPB1_\d+_\d+)',
        r'\1__\2',
        s
    )

    s = re.sub(
        r'(DQA1_\d+_\d+)_(DQB1_\d+_\d+)',
        r'\1__\2',
        s
    )

    return s


def display_token(value):

    v = norm_token(value)

    pair = re.fullmatch(
        r'(DPA1|DQA1)_(\d+)_(\d+)'
        r'__(DPB1|DQB1)_(\d+)_(\d+)',
        v
    )

    if pair:

        a, x, y, b, z, w = pair.groups()

        return (
            f'HLA-{a}*{x}:{y}-'
            f'HLA-{b}*{z}:{w}'
        )

    single = re.fullmatch(
        r'(DRB[1-5])_(\d+)_(\d+)',
        v
    )

    if single:

        a, x, y = single.groups()

        return f'HLA-{a}*{x}:{y}'

    return v


def supported_ii(loci):

    if not MIX_PWM.is_dir():
        raise FileNotFoundError(
            str(MIX_PWM)
        )

    available = {
        norm_token(p.stem): p.stem
        for p in MIX_PWM.glob('*.txt')
    }

    requested = []

    # DR alleles
    for locus in (
        'DRB1',
        'DRB3',
        'DRB4',
        'DRB5'
    ):

        requested.extend(
            token(a)
            for a in loci.get(locus, [])
        )

    # DP and DQ heterodimers
    for alpha, beta in [
        ('DPA1', 'DPB1'),
        ('DQA1', 'DQB1')
    ]:

        requested.extend(
            token(a) + '__' + token(b)
            for a in loci.get(alpha, [])
            for b in loci.get(beta, [])
        )

    matched = []
    unsupported = []
    displays = {}

    for r in dict.fromkeys(requested):

        key = norm_token(r)

        if key in available:

            actual = available[key]

            if actual not in matched:
                matched.append(actual)

            displays[key] = display_token(actual)

        else:
            unsupported.append(r)

    return matched, unsupported, displays


# ============================================================
# UniProt API
# ============================================================

def get_json(url, timeout=35):

    req = Request(
        url,
        headers={
            'User-Agent':
                'CureScience-TAA-Research/1.0',
            'Accept':
                'application/json'
        }
    )

    with urlopen(req, timeout=timeout) as r:
        return json.load(r)


def query_protein(gene, ensembl):

    # Search reviewed human proteins using exact gene symbol.

    q = (
        f'(gene_exact:{gene}) '
        f'AND (organism_id:9606) '
        f'AND (reviewed:true)'
    )

    url = (
        UNIPROT
        + '/search?'
        + urlencode({
            'query': q,
            'size': 50,
            'format': 'json'
        })
    )

    search = get_json(url)

    matches = []

    for item in search.get('results', []):

        names = []

        for g in item.get('genes', []):

            primary = g.get('geneName', {})
            primary_name = primary.get('value', '').strip().upper()

            if primary_name:
                names.append(primary_name)

            for synonym in g.get('synonyms', []):
                synonym_name = synonym.get('value', '').strip().upper()

                if synonym_name:
                    names.append(synonym_name)

        if gene not in names:
            continue

        refs = [
            r.get('id', '').split('.')[0]

            for r in item.get(
                'uniProtKBCrossReferences', []
            )

            if r.get('database') == 'Ensembl'
        ]

        accession = item.get(
            'primaryAccession'
        )

        if accession:

            matches.append((
                accession,
                ensembl in refs,
                gene in [
                    g.get('geneName', {}).get(
                        'value', ''
                    ).upper()

                    for g in item.get(
                        'genes', []
                    )
                ]
            ))

    # Prefer exact Ensembl match.
    exact_id = [
        m for m in matches
        if m[1]
    ]

    # Fallback to primary gene name.
    pool = (
        exact_id
        if exact_id
        else [
            m for m in matches
            if m[2]
        ]
    )

    pool = list(
        dict.fromkeys(
            m[0] for m in pool
        )
    )

    if len(pool) != 1:

        raise ValueError(
            f'UniProt mapping unresolved/ambiguous '
            f'for {gene} ({ensembl}): {pool}; '
            f'manual review needed'
        )

    accession = pool[0]

    doc = get_json(
        UNIPROT
        + '/'
        + quote(accession)
        + '.json'
    )

    primary_names = [
        g.get('geneName', {}).get(
            'value', ''
        ).upper()

        for g in doc.get(
            'genes', []
        )
    ]

    if gene not in primary_names:

        raise ValueError(
            f'UniProt primary gene mismatch: '
            f'{gene} vs {primary_names} '
            f'({accession})'
        )

    sequence = doc.get(
        'sequence', {}
    ).get(
        'value', ''
    ).upper()

    if (
        not sequence
        or not re.fullmatch(
            r'[A-Z]+',
            sequence
        )
    ):

        raise ValueError(
            f'Invalid UniProt protein sequence: '
            f'{gene}/{accession}'
        )

    if len(sequence) != doc.get(
        'sequence', {}
    ).get('length'):

        raise ValueError(
            f'UniProt sequence length mismatch: '
            f'{gene}/{accession}'
        )

    return {
        'gene': gene,
        'ensembl': ensembl,
        'accession': accession,
        'isoform': 'canonical',
        'length': len(sequence),
        'sequence': sequence,
        'sha256': hashlib.sha256(
            sequence.encode('ascii')
        ).hexdigest(),
        'downloaded_utc': datetime.now(
            timezone.utc
        ).isoformat(),
        'sequence_version': doc.get(
            'entryAudit', {}
        ).get('sequenceVersion')
    }


# ============================================================
# Protein sequence storage
# ============================================================

def wrap_fasta(seq, width=70):

    return '\n'.join(
        seq[i:i + width]
        for i in range(
            0,
            len(seq),
            width
        )
    )


def get_sequences(
    targets,
    work,
    refresh=False
):

    records_path = (
        work / 'protein_manifest.json'
    )

    cached = {}

    if (
        records_path.exists()
        and not refresh
    ):

        loaded = json.loads(
            records_path.read_text(
                encoding='utf-8'
            )
        )

        cached = {
            x['gene']: x

            for x in loaded.get(
                'records', []
            )

            if 'gene' in x
        }

    records = []
    failures = {}

    for item in targets:

        gene = item['gene']
        gid = item['ensembl']

        try:

            previous = cached.get(gene)

            if (
                previous
                and previous.get('ensembl') == gid
                and previous.get('sequence')
            ):

                seq = previous['sequence']

                actual_hash = hashlib.sha256(
                    seq.encode('ascii')
                ).hexdigest()

                if actual_hash != previous.get(
                    'sha256'
                ):

                    raise ValueError(
                        'Cached sequence hash mismatch; '
                        'use --refresh-proteins '
                        'after review'
                    )

                record = previous

                print(
                    f'[UniProt] cached {gene}: '
                    f'{record["accession"]}',
                    flush=True
                )

            else:

                record = query_protein(
                    gene,
                    gid
                )

                print(
                    f'[UniProt] downloaded {gene}: '
                    f'{record["accession"]} '
                    f'({record["length"]} aa)',
                    flush=True
                )

            records.append(record)

        except (
            HTTPError,
            URLError,
            TimeoutError,
            OSError,
            ValueError,
            KeyError
        ) as exc:

            failures[gene] = str(exc)

            print(
                f'[WARN] No validated protein '
                f'for {gene}: {exc}',
                file=sys.stderr,
                flush=True
            )

    payload = {
        'source':
            'UniProtKB reviewed human canonical',
        'records': records,
        'errors': failures
    }

    atomic_text(
        records_path,
        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2
        ) + '\n'
    )

    fasta = '\n'.join(

        (
            f'>{r["gene"]}'
            f'|UniProt={r["accession"]}'
            f'|Ensembl={r["ensembl"]}'
            f'|isoform=canonical'
            f'|length={r["length"]}'
            f'|sha256={r["sha256"]}\n'
            f'{wrap_fasta(r["sequence"])}'
        )

        for r in records
    )

    atomic_text(
        work / 'Targets_Gene_Sequence.txt',
        (fasta + '\n') if fasta else ''
    )

    return records, failures


# ============================================================
# Natural peptide generation
# ============================================================

def windows(sequence, size):

    # Skip peptides containing non-standard amino acids.
    # Never remove internal residues, because doing so
    # would create artificial peptide junctions.

    found = {}

    for i in range(
        len(sequence) - size + 1
    ):

        peptide = sequence[
            i:i + size
        ]

        if set(peptide) <= AA:

            # Keep first occurrence in protein.
            found.setdefault(
                peptide,
                i + 1
            )

    return found


# ============================================================
# MHC Class I prediction
# ============================================================

def run_mhc_i(seq_by_gene, alleles):

    from mhcflurry import (
        Class1AffinityPredictor
    )

    predictor = (
        Class1AffinityPredictor.load()
    )

    unique = sorted({
        p

        for seq in seq_by_gene.values()

        for p in seq
    })

    values = {}

    if not unique:

        return {
            g: []
            for g in seq_by_gene
        }

    for allele in alleles:

        # Process chunks to limit memory usage.
        for start in range(
            0,
            len(unique),
            3000
        ):

            chunk = unique[
                start:start + 3000
            ]

            predictions = (
                predictor.predict_to_dataframe(
                    allele=allele,
                    peptides=chunk
                )
            )

            if not {
                'peptide',
                'prediction'
            } <= set(predictions.columns):

                raise ValueError(
                    f'Unexpected MHCflurry '
                    f'columns for {allele}'
                )

            for row in predictions.itertuples(
                index=False
            ):

                ic50 = float(
                    row.prediction
                )

                if not math.isfinite(ic50):

                    raise ValueError(
                        f'Non-finite MHCflurry IC50 '
                        f'for {row.peptide} / {allele}'
                    )

                values[
                    (str(row.peptide), allele)
                ] = ic50

    output = {}

    for gene, peptides in (
        seq_by_gene.items()
    ):

        candidates = []

        for peptide, position in (
            peptides.items()
        ):

            hits = {
                a: values[(peptide, a)]

                for a in alleles

                if (
                    (peptide, a) in values
                    and values[
                        (peptide, a)
                    ] < I_CUTOFF
                )
            }

            if not hits:
                continue

            vals = list(
                hits.values()
            )

            candidates.append({
                'sequence': peptide,
                'position': position,
                'binding': hits,
                'mean_ic50':
                    sum(vals) / len(vals),
                'min_ic50':
                    min(vals)
            })

        # Same ranking logic as TME:
        # 1. More binding HLA alleles
        # 2. Lower mean IC50
        # 3. Lower minimum IC50

        candidates.sort(
            key=lambda x: (
                x['min_ic50'],
                x['sequence']
            )
        )

        output[gene] = candidates[
            :TOP_N
        ]

    return output


# ============================================================
# MHC Class II prediction
# ============================================================

def run_mhc_ii(
    seq_by_gene,
    allele_tokens,
    displays,
    out_dir
):

    all_peptides = sorted({
        p

        for seq in seq_by_gene.values()

        for p in seq
    })

    result = {
        gene: []

        for gene in seq_by_gene
    }

    if (
        not all_peptides
        or not allele_tokens
    ):
        return result

    required(MIX_EXE)

    inp = (
        out_dir
        / 'mixmhc2pred_input_15mers.txt'
    )

    out = (
        out_dir
        / 'mixmhc2pred_output.tsv'
    )

    atomic_text(
        inp,
        '\n'.join(all_peptides) + '\n'
    )

    command = [
        'conda',
        'run',
        '-n',
        MIX_CONDA_ENV,
        str(MIX_EXE),
        '--no_context',
        '--input',
        str(inp),
        '--output',
        str(out),
        '--alleles',
        *allele_tokens
    ]

    subprocess.run(
        command,
        check=True
    )

    required(out)

    # Read MixMHC2pred output.
    df = pd.read_csv(
        out,
        sep='\t',
        comment='#',
        dtype=str
    )

    if (
        df.empty
        or len(df.columns) == 1
    ):

        df = pd.read_csv(
            out,
            sep=r'\s+',
            comment='#',
            engine='python',
            dtype=str
        )

    def col(exact, contains):

        normalized = {
            str(x).lower(): x
            for x in df.columns
        }

        for x in exact:

            if x.lower() in normalized:
                return normalized[
                    x.lower()
                ]

        for column in df.columns:

            if any(
                x.lower()
                in column.lower()

                for x in contains
            ):
                return column

        raise ValueError(
            'Missing MixMHC2pred '
            'output column. Columns: '
            + str(list(df.columns))
        )

    pcol = col(
        ['Peptide', 'Sequence'],
        ['peptide', 'sequence']
    )

    rcol = col(
        ['%Rank_best', 'Rank_best'],
        ['rank_best']
    )

    acol = col(
        [
            'BestAllele',
            'Allele_best',
            'Best_allele'
        ],
        [
            'bestallele',
            'allele_best',
            'best_allele'
        ]
    )

    ccol = col(
        [
            'Core_best',
            'BestCore',
            'Core'
        ],
        [
            'core_best',
            'bestcore'
        ]
    )

    lookup = defaultdict(list)

    all_peptides_set = set(
        all_peptides
    )

    for row in df.itertuples(
        index=False,
        name=None
    ):

        cells = dict(
            zip(
                df.columns,
                row
            )
        )

        pep = str(
            cells[pcol]
        ).upper().strip()

        if pep not in all_peptides_set:
            continue

        try:

            rank = float(
                str(cells[rcol])
                .replace('%', '')
                .strip()
            )

        except ValueError:
            continue

        if not math.isfinite(rank):
            continue

        allele = str(
            cells[acol]
        ).strip()

        core = str(
            cells[ccol]
        ).strip().upper()

        if (
            not allele
            or allele.lower() == 'nan'
            or not core
            or core.lower() == 'nan'
        ):
            continue

        key = norm_token(
            allele
        )

        if key not in displays:

            raise ValueError(
                f'Unexpected best HLA allele '
                f'reported by MixMHC2pred: '
                f'{allele}'
            )

        lookup[pep].append((
            rank,
            displays[key],
            core
        ))

    # Map every peptide back to its source gene.
    for gene, peptide_positions in (
        seq_by_gene.items()
    ):

        ranked = []

        for peptide, position in (
            peptide_positions.items()
        ):

            if peptide not in lookup:
                continue

            rank, allele, core = min(
                lookup[peptide],
                key=lambda x: (
                    x[0],
                    x[1],
                    x[2]
                )
            )

            ranked.append({
                'sequence': peptide,
                'position': position,
                'percent_rank': rank,
                'hla': allele,
                'binding_core': core
            })

        ranked.sort(
            key=lambda x: (
                x['percent_rank'],
                x['sequence'],
                x['hla']
            )
        )

        result[gene] = ranked[
            :TOP_N
        ]

    return result


# ============================================================
# XML output
# ============================================================

def write_output(
    dest,
    targets,
    records,
    errors,
    hla_i,
    hla_ii,
    unsupported,
    i_results,
    ii_results,
    data4_name,
    data3_name
):

    by_gene = {
        r['gene']: r
        for r in records
    }

    root = ET.Element(
        'TAAPeptideTargets',
        {
            'rna_patient':
                data4_name,

            'hla_patient':
                data3_name,

            'protein_reference':
                'UniProtKB_reviewed_canonical',

            'classI_model':
                'MHCflurry_Class1AffinityPredictor',

            'classII_model':
                'MixMHC2pred-2.0',

            'classI_ic50_cutoff_nM':
                '500',

            'max_peptides_per_gene_per_class':
                str(TOP_N),

            'anchor_engineering':
                'false',

            'candidate_status':
                'binding_prediction_only'
        }
    )

    # HLA information
    hla = ET.SubElement(
        root,
        'HLATyping'
    )

    ET.SubElement(
        hla,
        'ClassI',
        alleles=';'.join(hla_i)
    )

    ET.SubElement(
        hla,
        'ClassII',
        alleles=';'.join(hla_ii),
        unsupported=';'.join(unsupported)
    )

    # Gene targets
    for target in targets:

        gene = target['gene']

        r = by_gene.get(gene)

        node = ET.SubElement(
            root,
            'Gene',
            {
                'name': gene,
                'ensembl_id':
                    target['ensembl'],
                'TPM':
                    target['tpm'],
                'protein_status':
                    (
                        'resolved'
                        if r
                        else 'unresolved'
                    )
            }
        )

        if r:

            ET.SubElement(
                node,
                'Protein',
                {
                    'uniprot_accession':
                        r['accession'],

                    'isoform':
                        r['isoform'],

                    'sequence_length':
                        str(r['length']),

                    'sequence_sha256':
                        r['sha256'],

                    'sequence_version':
                        str(
                            r.get(
                                'sequence_version'
                            ) or ''
                        )
                }
            )

        else:

            ET.SubElement(
                node,
                'UnresolvedReason'
            ).text = errors.get(
                gene,
                'Unknown protein-mapping error'
            )

        # ----------------------------------------
        # Class I
        # ----------------------------------------

        cl1 = ET.SubElement(
            node,
            'MHCClassI',
            count=str(
                len(
                    i_results.get(
                        gene, []
                    )
                )
            )
        )

        for idx, x in enumerate(
            i_results.get(
                gene, []
            ),
            1
        ):

            ep = ET.SubElement(
                cl1,
                'Peptide',
                rank=str(idx),
                sequence=x['sequence'],
                protein_start_1based=str(
                    x['position']
                ),
                mean_ic50_nM=(
                    f'{x["mean_ic50"]:.6g}'
                ),
                minimum_ic50_nM=(
                    f'{x["min_ic50"]:.6g}'
                ),
                binding_hla_count=str(
                    len(x['binding'])
                )
            )

            for allele, affinity in sorted(
                x['binding'].items(),
                key=lambda y: (
                    y[1],
                    y[0]
                )
            ):

                ET.SubElement(
                    ep,
                    'Binding',
                    hla=allele,
                    ic50_nM=f'{affinity:.6g}'
                )

        # ----------------------------------------
        # Class II
        # ----------------------------------------

        cl2 = ET.SubElement(
            node,
            'MHCClassII',
            count=str(
                len(
                    ii_results.get(
                        gene, []
                    )
                )
            )
        )

        for idx, x in enumerate(
            ii_results.get(
                gene, []
            ),
            1
        ):

            ET.SubElement(
                cl2,
                'Peptide',
                rank=str(idx),
                sequence=x['sequence'],
                protein_start_1based=str(
                    x['position']
                ),
                hla=x['hla'],
                percent_rank=(
                    f'{x["percent_rank"]:.6g}'
                ),
                binding_core=x[
                    'binding_core'
                ]
            )

    ET.indent(
        root,
        space='  '
    )

    dest.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    temp = dest.with_suffix(
        '.xml.tmp'
    )

    ET.ElementTree(root).write(
        temp,
        encoding='utf-8',
        xml_declaration=True
    )

    temp.replace(dest)


# ============================================================
# Main
# ============================================================

def main(argv=None):

    ap = argparse.ArgumentParser(
        description=__doc__
    )

    ap.add_argument(
        '--rna-patient',
        help=(
            'Folder under /data4/NewOutput '
            '(prompt if omitted)'
        )
    )

    ap.add_argument(
        '--hla-patient',
        help=(
            'Folder under /data3/PipelineOutput '
            '(prompt if omitted)'
        )
    )

    ap.add_argument(
        '--refresh-proteins',
        action='store_true',
        help=(
            'Re-query UniProt rather than reuse '
            'validated local protein cache'
        )
    )

    args = ap.parse_args(argv)

    # Ask for data4 patient folder.
    rna_patient = patient_folder(
        args.rna_patient
        or input(
            'Enter data4 RNA patient folder name: '
        )
    )

    # Ask for data3 HLA patient folder.
    hla_patient = patient_folder(
        args.hla_patient
        or input(
            'Enter data3 HLA patient folder name: '
        )
    )

    base = (
        DATA4_ROOT
        / rna_patient
        / 'rnaseq-results'
    )

    inp = (
        base
        / 'TAA.xml'
    )

    work = (
        base
        / 'TAA'
    )

    prediction_dir = (
        work
        / 'peptide_prediction'
    )

    prediction_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    hla_root = (
        DATA3_ROOT
        / hla_patient
        / 'HLAclassification'
    )

    # ----------------------------------------
    # Step 1: Load overexpressed TAA targets
    # ----------------------------------------

    targets = read_taa(inp)

    # ----------------------------------------
    # Step 2: Load patient HLA
    # ----------------------------------------

    hla_i = parse_hla_i(
        hla_root
        / 'final_ClassI.HLAgenotype4digits'
    )

    hla_ii = parse_hla_ii(
        hla_root
        / 'final_ClassII.HLAgenotype4digits'
    )

    ii_tokens, unsupported, displays = (
        supported_ii(hla_ii)
    )

    print(
        '[INPUT] TAA genes:',
        len(targets),
        'HLA-I:',
        hla_i,
        'HLA-II supported:',
        ii_tokens,
        flush=True
    )

    if unsupported:

        print(
            '[WARN] Unsupported Class II '
            'HLA alleles:',
            unsupported,
            file=sys.stderr
        )

    if not ii_tokens:

        print(
            '[WARN] None of the Class II alleles '
            'supported; output empty Class II results',
            file=sys.stderr
        )

    # ----------------------------------------
    # Step 3: Protein retrieval
    # ----------------------------------------

    if not targets:

        print(
            '[INFO] No overexpressed targets in '
            'TAA.xml; producing empty XML'
        )

        records = []
        errors = {}
        i_results = {}
        ii_results = {}

        atomic_text(
            work / 'Targets_Gene_Sequence.txt',
            ''
        )

        atomic_text(
            work / 'protein_manifest.json',
            json.dumps(
                {
                    'records': [],
                    'errors': {}
                },
                indent=2
            ) + '\n'
        )

    else:

        records, errors = get_sequences(
            targets,
            work,
            args.refresh_proteins
        )

        # ----------------------------------------
        # Step 4: Generate natural peptides
        # ----------------------------------------

        seq_i = {
            r['gene']: windows(
                r['sequence'],
                I_LENGTH
            )
            for r in records
        }

        seq_ii = {
            r['gene']: windows(
                r['sequence'],
                II_LENGTH
            )
            for r in records
        }

        # ----------------------------------------
        # Step 5: Predict MHC binding
        # ----------------------------------------

        if records:

            i_results = run_mhc_i(
                seq_i,
                hla_i
            )

            ii_results = run_mhc_ii(
                seq_ii,
                ii_tokens,
                displays,
                prediction_dir
            )

        else:

            i_results = {}
            ii_results = {}

    # ----------------------------------------
    # Step 6: Write final XML
    # ----------------------------------------

    dest = (
        work
        / 'TAA_Peptide.xml'
    )

    write_output(
        dest,
        targets,
        records,
        errors,
        hla_i,
        ii_tokens,
        unsupported,
        i_results,
        ii_results,
        rna_patient,
        hla_patient
    )

    print(
        '[DONE] sequences:',
        work / 'Targets_Gene_Sequence.txt'
    )

    print(
        '[DONE] peptide XML:',
        dest
    )

    print(
        '[DONE] validated protein records:',
        len(records),
        'unresolved:',
        len(errors)
    )

    print(
        '[DONE] Class I peptides:',
        sum(
            map(
                len,
                i_results.values()
            )
        )
    )

    print(
        '[DONE] Class II peptides:',
        sum(
            map(
                len,
                ii_results.values()
            )
        )
    )


if __name__ == '__main__':

    try:
        main()

    except KeyboardInterrupt:
        sys.exit(
            '\n[STOPPED]'
        )

    except Exception as exc:

        sys.exit(
            f'[ERROR] {type(exc).__name__}: {exc}'
        )
