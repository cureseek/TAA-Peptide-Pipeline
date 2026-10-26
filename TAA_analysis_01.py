import os
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMEXPR_MAX_THREADS'):
    os.environ[_v] = '40'
import argparse
import gzip
import re
import sys
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
import pandas as pd

REF_ROOT = Path('/data4/local-bioinformatics/RNA-CNV')
GDC_MATRIX = REF_ROOT / 'public-rnaseq/gdc_tpm_matrix.parquet'
GDC_MANIFEST = REF_ROOT / 'public-rnaseq/gdc_manifest.tsv'
GTF = REF_ROOT / 'GenCode/gencode.v50.primary_assembly.annotation.gtf.gz'
GTEX_TPM = Path('/data3/Reference/GTEx/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz')
GTEX_META = Path('/data3/Reference/GTEx/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt')
OUTPUT_ROOT = Path('/data4/NewOutput')
MIN_COHORT_TPM = 1.0  
CUTOFF = 80.0
MIN_PATIENT_TPM = 1.0
GENES_CSV = '''
WT1,
CTAG1B,
MUC1,
TERT,
ERBB2,
MAGEA3,
MSLN,
BIRC5,
CEACAM5,
KLK3,
ACPP,
FOLH1,
PMEL,
PRAME
'''

GENES = list(dict.fromkeys(
    s.strip().upper()
    for s in GENES_CSV.split(',')
    if s.strip()
))



CANCER_TYPE_MAP = dict(zip('LAML ACC BLCA LGG GBM BRCA CESC CHOL COAD READ ESCA HNSC KICH KIRC KIRP LIHC LUAD LUSC DLBC MESO OV PAAD PCPG PRAD SARC SKCM STAD TGCT THYM THCA UCEC UCS UVM'.split(),
                            ('TCGA-'+s for s in 'LAML ACC BLCA LGG GBM BRCA CESC CHOL COAD READ ESCA HNSC KICH KIRC KIRP LIHC LUAD LUSC DLBC MESO OV PAAD PCPG PRAD SARC SKCM STAD TGCT THYM THCA UCEC UCS UVM'.split())))


def require(path):
    if not path.is_file():
        raise FileNotFoundError(f'Required file missing: {path}')
    return path


def select_number(label, choices):
    print('\nSelect '+label+':')
    for i, name in enumerate(choices, 1):
        print(f'  {i:2}. {name}')
    while True:
        answer=input(f'{label} number (1-{len(choices)}): ').strip()
        if answer.isdecimal() and 1 <= int(answer) <= len(choices):
            return choices[int(answer)-1]
        print('Please enter a listed number.')


def symbol_mapping():
    result = {}
    all_pc = set()
    with gzip.open(require(GTF), 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            parts=line.rstrip('\n').split('\t')
            if len(parts)<9 or parts[2]!='gene': continue
            attrs=parts[8]
            ids=re.search(r'gene_id "([^"]+)',attrs)
            name=re.search(r'gene_name "([^"]+)',attrs)
            typ=re.search(r'gene_type "([^"]+)',attrs)
            if ids and name and typ and typ.group(1)=='protein_coding':
                gid=ids.group(1).split('.')[0]
                all_pc.add(gid)
                result.setdefault(name.group(1).upper(),[]).append(gid)
    return result, all_pc


def load_salmon(patient):
    path=OUTPUT_ROOT/patient/'rnaseq-results'/'star_salmon'/f'{patient}_tumor'/'quant.genes.sf'
    require(path)
    tab=pd.read_csv(path,sep='\t',usecols=['Name','TPM'])
    tab['gene_id']=tab['Name'].astype(str).str.split('.').str[0]
    tab['TPM']=pd.to_numeric(tab['TPM'],errors='coerce')
    if tab['TPM'].dropna().lt(0).any(): raise ValueError('Negative patient TPM found')
    tab=tab.dropna(subset=['TPM']).sort_values('TPM',ascending=False).drop_duplicates('gene_id')
    tpm=tab.set_index('gene_id')['TPM']
    print(f'[INFO] Existing Salmon TPM: {path}; {len(tpm)} gene IDs')
    return tpm


def gdc_read(columns):
    sys.path.insert(0,str(REF_ROOT/'public-rnaseq'))
    from matrix_cache import read_matrix_columns
    return read_matrix_columns(require(GDC_MATRIX),columns)


def factor(patient_tpm, matrix, pc_ids, label):
    # Exact median-of-ratios rule from Michael's first script.
    ids=sorted(pc_ids & set(patient_tpm.index) & set(matrix.index))
    if not ids: raise ValueError('No shared protein-coding genes for '+label)
    cohort_median=matrix.loc[ids].apply(pd.to_numeric,errors='coerce').median(axis=1)
    tpm=patient_tpm.loc[ids]
    use=cohort_median.ge(MIN_COHORT_TPM)&tpm.gt(0)
    result=float((tpm[use]/cohort_median[use]).median())
    if not np.isfinite(result) or result<=0: raise ValueError(f'Invalid {label} scaling factor: {result}')
    print(f'[INFO] {label} median-of-ratios factor: {result:.5f} (n={int(use.sum())})')
    return result


def gdc_percentile(scaled_tpm, matrix):
    # Michael's PANEL comparison: <= ties count fully; log2(TPM+1) monotonic transform.
    ref=np.log2(matrix.apply(pd.to_numeric,errors='coerce')+1.0)
    value=np.log2(scaled_tpm+1.0)
    n=ref.notna().sum(axis=1)
    return (ref.le(value,axis=0).sum(axis=1)/n*100).where(n>0,np.nan),n


def gtex_percentile(patient_tpm, ids, tissue, meta):
    # Exact second-script rule: raw patient TPM, <= ties fully; NO FC filter, NO GDC scaling.
    samples=(meta.loc[meta['SMTSD'].astype(str).str.strip().eq(tissue),'SAMPID']
             .dropna().astype(str).drop_duplicates().tolist())
    if not samples: raise ValueError('GTEx tissue has no samples: '+tissue)
    with gzip.open(require(GTEX_TPM),'rt') as fh:
        header=fh.readline(); fh.readline()
        columns=fh.readline().rstrip('\n').split('\t')
    selected=[s for s in samples if s in set(columns)]
    if not selected: raise ValueError('GTEx samples absent from GCT: '+tissue)
    # Read only chosen tissue columns to avoid loading the enormous full GTEx expression matrix.
    with gzip.open(GTEX_TPM,'rt') as fh:
        tab=pd.read_csv(fh,sep='\t',skiprows=2,usecols=['Name',*selected],low_memory=False)
    tab['gene_id']=tab['Name'].astype(str).str.split('.').str[0]
    tab=tab[tab['gene_id'].isin(ids)].copy()
    if tab.empty: raise ValueError('No selected gene IDs present in GTEx reference')
    # Second script converts to numeric, fillna(0), and sums duplicate unversioned gene IDs.
    numeric=tab[selected].apply(pd.to_numeric,errors='coerce').fillna(0.0)
    numeric.index=tab['gene_id'].to_numpy()
    ref=numeric.groupby(level=0).sum()
    ref=ref.reindex(ids)
    rank=ref.le(patient_tpm.reindex(ids),axis=0).sum(axis=1)/ref.notna().sum(axis=1)*100
    rank=rank.where(ref.notna().sum(axis=1)>0,np.nan)
    print(f'[INFO] GTEx tissue={tissue}; sample count={len(selected)}')
    return rank,ref.notna().sum(axis=1),len(selected)


def xml_write(rows,patient,cancer,tissue,matched_n,pan_n,gtex_n,factors,destination):
    root=ET.Element('TAAOverexpression',{
        'patient':patient,'cancer_type':cancer,'normal_tissue':tissue,
        'matched_samples':str(matched_n),'pancancer_samples':str(pan_n),
        'normal_samples':str(gtex_n),'percentile_cutoff':'80',
        'matched_scaling_factor':f'{factors[0]:.8g}',
        'pancancer_scaling_factor':f'{factors[1]:.8g}',
        'passed_genes':str(len(rows))})
    for row in rows.itertuples(index=False):
        ET.SubElement(root,'Gene',{
            'name':row.gene_name,'ensembl_id':row.gene_id,
            'TPM':f'{row.TPM:.6g}',
            'matched_cohort_percentile':f'{row.matched_percentile:.4f}',
            'PanCancer_percentile':f'{row.pancancer_percentile:.4f}',
            'GTEx_percentile':f'{row.gtex_percentile:.4f}',
            'overexpressed':'true'})
    tree=ET.ElementTree(root)
    ET.indent(tree,space='  ')
    destination.parent.mkdir(parents=True,exist_ok=True)
    tree.write(destination,encoding='utf-8',xml_declaration=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--patient',help='Patient folder; otherwise prompted')
    parser.add_argument('--cancer-type',help='TCGA project short code; otherwise select number')
    parser.add_argument('--tissue',help='Exact GTEx SMTSD tissue; otherwise select number')
    args=parser.parse_args()
    patient=(args.patient or input('Enter patient folder name: ')).strip()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',patient):
        raise ValueError('Invalid patient folder name')
    cancer_options=sorted(CANCER_TYPE_MAP)
    cancer=(args.cancer_type or select_number('TCGA cancer type',cancer_options)).upper()
    if cancer not in CANCER_TYPE_MAP: raise ValueError('Invalid cancer code: '+cancer)
    meta=pd.read_csv(require(GTEX_META),sep='\t',low_memory=False)
    if not {'SAMPID','SMTSD'}.issubset(meta.columns): raise ValueError('GTEx metadata missing SAMPID/SMTSD')
    tissues=sorted(meta['SMTSD'].dropna().astype(str).str.strip().unique())
    tissue=args.tissue or select_number('GTEx normal tissue',tissues)
    if tissue not in tissues: raise ValueError('Unsupported GTEx SMTSD tissue: '+tissue)
    patient_tpm=load_salmon(patient)
    mapping,pc_ids=symbol_mapping()
    entries=[]; unresolved=[]
    for gene in GENES:
        ids=mapping.get(gene,[])
        if not ids: unresolved.append(gene)
        else:
            for gid in ids: entries.append((gene,gid))
    if unresolved: print('[WARN] Unresolved protein-coding symbols:',', '.join(unresolved))
    if not entries: raise ValueError('No TAA genes mapped to GENCODE v50')
    # One gene per Ensembl ID, consistent with the first script's drop_duplicates('gene_id').
    panel=pd.DataFrame(entries,columns=['gene_name','gene_id']).drop_duplicates('gene_id').set_index('gene_id')
    manifest=pd.read_csv(require(GDC_MANIFEST),sep='\t')
    for col in ('column_label','cancer_type','sample_type'):
        if col not in manifest: raise ValueError('GDC manifest missing '+col)
    matched_cols=manifest.loc[manifest['cancer_type'].eq(CANCER_TYPE_MAP[cancer]) & manifest['sample_type'].eq('Primary Tumor'),'column_label'].tolist()
    pan_cols=manifest.loc[manifest['sample_type'].eq('Primary Tumor'),'column_label'].tolist()
    if not matched_cols or not pan_cols: raise ValueError('No GDC reference samples for requested cohorts')
    print(f'[INFO] GDC matched={len(matched_cols)}, PanCancer={len(pan_cols)}; reading cached reference...')
    matched=gdc_read(matched_cols)
    pan=gdc_read(pan_cols)
    fm=factor(patient_tpm,matched,pc_ids,'matched')
    fp=factor(patient_tpm,pan,pc_ids,'PanCancer')
    ids=sorted(set(panel.index)&set(patient_tpm.index)&set(matched.index)&set(pan.index))
    if not ids: raise ValueError('No shared TAA gene IDs in patient/GDC')
    print(f'[INFO] Evaluating {len(ids)} TAA gene IDs')
    tpm=patient_tpm.reindex(ids)
    m_rank,m_n=gdc_percentile(tpm/fm,matched.reindex(ids))
    p_rank,p_n=gdc_percentile(tpm/fp,pan.reindex(ids))
    g_rank,g_n,gtex_samples=gtex_percentile(patient_tpm,ids,tissue,meta)
    result=panel.reindex(ids).copy()
    result['TPM']=tpm
    result['matched_percentile']=m_rank
    result['pancancer_percentile']=p_rank
    result['gtex_percentile']=g_rank
    result['is_hit'] = (
    (result['TPM'] >= MIN_PATIENT_TPM) &
    (result['matched_percentile'] >= CUTOFF) &
    (result['pancancer_percentile'] >= CUTOFF) &
    (result['gtex_percentile'] >= CUTOFF)
)
    result=result.reset_index().sort_values(['is_hit','matched_percentile','pancancer_percentile'],ascending=[False,False,False])
    hits=result[result.is_hit]
    dest=OUTPUT_ROOT/patient/'rnaseq-results'/'TAA.xml'
    xml_write(hits,patient,CANCER_TYPE_MAP[cancer],tissue,len(matched_cols),len(pan_cols),gtex_samples,(fm,fp),dest)
    # Optional diagnostic: does not change selection or final XML.
    diagnostic=dest.with_name('TAA_all_genes.csv')
    result.to_csv(diagnostic,index=False)
    print(f'[DONE] {len(hits)}/{len(result)} evaluated TAA gene IDs pass all three percentiles >=80')
    print('[DONE] XML:',dest)
    print('[INFO] QC CSV:',diagnostic)

if __name__=='__main__':
    try: main()
    except (OSError,ValueError,ImportError,KeyError) as exc:
        sys.exit(f'[ERROR] {exc}')
