# TAA-Peptide-Pipeline

A patient-specific bioinformatics framework for **tumor-associated antigen (TAA) selection and peptide-binding prediction**, integrating tumor RNA-seq expression, **TCGA/GDC** cancer reference cohorts, **GTEx** normal-tissue expression, **GENCODE v50**, **UniProtKB/Swiss-Prot**, and patient-specific **MHC class I and class II** genotypes.

The project prioritizes a curated **14-gene TAA panel** through relative expression benchmarking and a minimum tumor-expression threshold, then generates naturally occurring protein-derived peptides and ranks their predicted HLA binding using **MHCflurry** and **MixMHC2pred 2.0**.

A central objective is to connect **patient-level transcriptional evidence** with **protein-sequence provenance** and **HLA-specific peptide predictions**, while explicitly distinguishing computational binding predictions from evidence of endogenous presentation, immunogenicity, tumor specificity, or therapeutic efficacy.

**Project status (2026-10-26): Active development and validation.** The two-stage TAA screening and peptide prediction implementation is available. Cross-module patient/HLA provenance checks and independent experimental validation remain important quality-control priorities.

---

## Project Objectives

The main objectives are:

1. Evaluate a literature-/expert-curated panel of 14 tumor-associated antigen genes in each patient's RNA-seq data.
2. Compare patient gene expression with matched-cancer and pan-cancer TCGA/GDC cohorts and a selected normal-tissue GTEx cohort.
3. Exclude genes with undetectable or very low patient expression before downstream peptide prediction.
4. Retrieve reviewed human reference protein sequences from UniProtKB and retain accession, Ensembl ID, sequence version, and checksum information.
5. Generate native protein-derived 9-mers and 15-mers without artificial anchor engineering.
6. Predict class I binding affinity across the patient's available HLA-A/B/C alleles and prioritize peptides by **minimum predicted IC50**.
7. Predict class II binding using supported HLA-DR, -DP, and -DQ models and prioritize peptides by **minimum predicted percent rank**.
8. Produce structured per-patient XML and diagnostic outputs suitable for downstream research reporting and quality review.
9. Identify limitations caused by transcript isoforms, gene-family homology, RNA quantification, HLA inference, and computational presentation models.

This pipeline is designed for **research candidate prioritization**, not as a standalone clinical decision system.

---

## Pipeline Architecture

The workflow consists of an upstream RNA quantification prerequisite followed by two TAA-specific analysis stages.

```text
                  Tumor RNA-seq FASTQ (paired-end)
                                |
                                v
                       nf-core/rnaseq
                       STAR + Salmon
                           GRCh38
                                |
                                v
                         quant.genes.sf
                      Patient gene-level TPM
                                |
                 +--------------+--------------+
                 |              |              |
                 v              v              v
           TCGA matched     TCGA pan-cancer    GTEx selected
           cancer cohort      tumor cohort      normal tissue
                 |              |              |
                 +--------------+--------------+
                                |
                                v
                     TAA_analysis_01.py
                  Curated 14-gene TAA panel
                    GENCODE gene mapping
                 Cohort-specific percentiles
                                |
                        Four eligibility rules
                         Patient TPM >= 1
                   TCGA matched percentile >= 80
                   Pan-cancer percentile >= 80
                     GTEx percentile >= 80
                                |
                                v
                         TAA.xml
                       TAA_all_genes.csv
                                |
                 +--------------+--------------+
                 |                             |
                 v                             v
         UniProtKB reviewed             Patient HLA-I / HLA-II
          canonical proteins             genotype input files
                 |                             |
                 +--------------+--------------+
                                |
                                v
                     TAA_analysis_02.py
                      Native peptide windows
                                |
                   +------------+------------+
                   |                         |
                   v                         v
            MHC Class I                 MHC Class II
             MHCflurry                 MixMHC2pred 2.0
               9-mers                     15-mers
             IC50 < 500 nM            Minimum %Rank
          Minimum IC50 ranking        %Rank ranking
                   |                         |
                   +------------+------------+
                                |
                                v
                 Top 10 peptides per gene / class
                                |
                                v
                       TAA_Peptide.xml
                     protein_manifest.json
                   Targets_Gene_Sequence.txt
                                |
                                v
                    Downstream research report
```

**Important:** The TAA scripts **read existing** STAR–Salmon expression results. `TAA_analysis_01.py` does not itself run RNA alignment or quantification. A separate local RNA-seq workflow must produce `quant.genes.sf` first.

---

## Curated TAA Target Panel

The current implementation evaluates **14 gene symbols representing 12 antigen categories**. The prostate-antigen category consists of three distinct genes (PSA, PAP, and PSMA).

| Antigen | Gene symbol | Notes |
|---|---|---|
| Wilms tumor 1 | `WT1` | Transcription factor; explored across multiple tumor types |
| NY-ESO-1 | `CTAG1B` | Cancer-testis antigen; highly homologous locus warrants RNA-mapping QC |
| Mucin 1 | `MUC1` | Variable tandem-repeat region and isoform considerations |
| Telomerase reverse transcriptase | `TERT` | Telomerase catalytic subunit |
| HER2/neu | `ERBB2` | Receptor tyrosine kinase |
| MAGE-A3 | `MAGEA3` | Cancer-testis antigen |
| Mesothelin | `MSLN` | Cell-surface-associated antigen |
| Survivin | `BIRC5` | Apoptosis/cell-division-associated protein |
| Carcinoembryonic antigen | `CEACAM5` | Cell adhesion-related antigen |
| Prostate-specific antigen (PSA) | `KLK3` | Prostate-associated category |
| Prostatic acid phosphatase (PAP) | `ACPP` | Prostate-associated category |
| Prostate-specific membrane antigen (PSMA) | `FOLH1` | Prostate-associated category |
| gp100 | `PMEL` | Melanocytic-lineage antigen |
| PRAME | `PRAME` | Cancer-testis-associated antigen |

The **same 14 genes are initially evaluated for every patient**; the final eligible set is determined by the patient's expression and reference-cohort comparisons. Inclusion in the panel does **not** establish that an antigen is safe, tumor-exclusive, or clinically actionable.

---

## Module 1: RNA Expression-Based TAA Selection

**Script:** `TAA_analysis_01.py`

### Input datasets

| Input | Function |
|---|---|
| Salmon `quant.genes.sf` | Patient tumor gene-level TPM |
| GENCODE v50 GRCh38 GTF | Maps gene symbols to protein-coding Ensembl gene IDs |
| GDC TPM matrix and manifest | TCGA matched-cancer and pan-cancer tumor cohorts |
| GTEx v10 TPM matrix and sample metadata | Selected normal-tissue reference cohort |

The patient is assigned a TCGA cancer type and an appropriate GTEx normal-tissue category. The normal-tissue choice is a *reference comparison*, not necessarily an exact match to the biopsy site.

### Reference comparison

For the TCGA matched and pan-cancer comparisons, the implementation estimates **separate cohort-specific median-of-ratios scaling factors** from shared protein-coding genes with eligible cohort median TPM and positive patient TPM. Patient TPM is divided by the appropriate factor before determining the corresponding TCGA percentile.

The GTEx comparison uses **unscaled patient TPM** against the selected tissue's GTEx reference TPM values.

Percentiles count ties using `reference <= patient`, so all-zero reference distributions can produce high percentile ranks even when patient TPM is zero. This is why the minimum patient-expression gate is required.

### Current selection thresholds

A gene is eligible only when **all four** conditions hold:

| Criterion | Threshold |
|---|---|
| Patient tumor gene expression | **TPM >= 1.0** |
| TCGA matched-cancer percentile | **>= 80** |
| TCGA pan-cancer percentile | **>= 80** |
| GTEx selected-tissue percentile | **>= 80** |

The **1.0 TPM cutoff is an initial project-specific computational screening threshold**, not a validated threshold for immunogenicity or tumor antigen presentation. The algorithm does not currently establish a tumor-to-normal fold-change minimum or a peptide-level normal-tissue safety margin.

### Outputs

- `TAA.xml`: genes satisfying all current selection rules, including expression and percentile attributes.
- `TAA_all_genes.csv`: evaluated gene IDs, patient TPM, percentile values, and selection status for review.

If no genes pass, an empty target set is a valid analytical outcome.

---

## Module 2: Reference Proteins and HLA Binding Prediction

**Script:** `TAA_analysis_02.py`

### Inputs and sample mapping

This module reads `TAA.xml` from the selected **Data4 RNA patient** and independently reads HLA files from the selected **Data3 HLA patient**:

```text
/data3/PipelineOutput/<HLA_PATIENT>/HLAclassification/
    final_ClassI.HLAgenotype4digits
    final_ClassII.HLAgenotype4digits
```

**Patient identity must be independently confirmed when the Data3 and Data4 folder names differ.** The available command-line interface accepts both names and does not by itself prove that they refer to the same patient.

### UniProt reference selection

The protein retrieval routine:

1. Searches reviewed human UniProtKB entries using the exact gene symbol.
2. Prefers an entry with the corresponding Ensembl cross-reference.
3. Falls back to a unique matching primary gene name when an exact Ensembl cross-reference is unavailable.
4. Rejects unresolved or ambiguous results instead of choosing arbitrarily.
5. Stores the canonical reference sequence, UniProt accession, Ensembl ID, protein length, sequence version, and SHA256 checksum.
6. Reuses a previously validated local sequence cache unless protein refresh is requested.

### MHC Class I

| Parameter | Current configuration |
|---|---|
| Input sequence | Native protein 9-mer windows |
| Predictor | MHCflurry `Class1AffinityPredictor` |
| Alleles | Patient HLA-A, HLA-B, HLA-C alleles parsed from Class I file |
| Candidate binding cutoff | Predicted IC50 **< 500 nM** |
| Ranking | **Minimum predicted IC50 ascending**; sequence as tie-breaker |
| Output | Up to **10 peptides per gene** |

Mean IC50 and the number of HLA alleles with predictions below the cutoff are recorded, but **not used as the primary ranking criterion**. Older comments in the script may still describe a prior ranking strategy; the implemented sort key controls actual behavior.

### MHC Class II

| Parameter | Current configuration |
|---|---|
| Input sequence | Native protein 15-mer windows |
| Predictor | MixMHC2pred **2.0** (`--no_context`) |
| Alleles | Supported DR alleles and modeled DP/DQ alpha–beta combinations |
| Ranking | **Lowest predicted percent rank** |
| Output | Up to **10 peptides per gene** |

The implementation reports a best supported HLA assignment and binding core. **Percent rank is not IC50** and should not be compared numerically with class I binding affinity.

DP/DQ combinations are constructed from observed alpha and beta allele lists; genotype phase and biological pairing compatibility require additional scrutiny. Unsupported model alleles are recorded separately.

### Peptide generation rules

- Peptides are taken from **unmodified human reference protein sequences**.
- No anchor residue substitution or synthetic peptide engineering is performed.
- Peptides containing nonstandard amino acids are excluded.
- Duplicate sequences within a protein are represented by their first occurrence.
- All reported peptide positions are **1-based** relative to the selected reference sequence.

---

## Running the Pipeline

### Prerequisites

This code was developed within a local Linux research environment and uses local reference datasets and tools. It is **not yet packaged as a portable, one-command installer**.

Major dependencies include Python (`pandas`, `numpy`), STAR–Salmon quantification from `nf-core/rnaseq`, MHCflurry with an installed prediction model, MixMHC2pred 2.0, an available `mhcii` Conda environment, UniProt API access, and the GENCODE/GDC/GTEx reference resources configured in the scripts.

Before reuse on another system, replace or parameterize local absolute paths and verify software/model versions and reference-data licenses.

### Step 0 — Generate patient RNA expression

An upstream RNA-seq process must already have produced the gene-level Salmon file:

```text
/data4/NewOutput/<RNA_PATIENT>/rnaseq-results/
    star_salmon/<RNA_PATIENT>_tumor/quant.genes.sf
```

### Step 1 — Run TAA selection

```bash
python TAA_analysis_01.py \
    --patient <RNA_PATIENT> \
    --cancer-type GBM \
    --tissue "Brain - Cortex"
```

`GBM` and `Brain - Cortex` are **examples only**. Choose the correct disease-specific TCGA code and an appropriately justified GTEx reference tissue for the patient being analyzed. Running the script without optional arguments opens interactive prompts.

### Step 2 — Run peptide prediction

```bash
python TAA_analysis_02.py \
    --rna-patient <RNA_PATIENT> \
    --hla-patient <HLA_PATIENT>
```

To force a fresh UniProt retrieval instead of using existing protein cache records:

```bash
python TAA_analysis_02.py \
    --rna-patient <RNA_PATIENT> \
    --hla-patient <HLA_PATIENT> \
    --refresh-proteins
```

Only use the refresh option when updated sequence retrieval is desired and its impact on reproducibility has been considered.

### Recommended execution order

If the expression thresholds, reference matrices, target list, or patient input change, regenerate `TAA.xml` **before** regenerating `TAA_Peptide.xml`. Downstream reporting should consume the newest validated output for the *same* verified patient.

---

## Output Structure

```text
/data4/NewOutput/<RNA_PATIENT>/
  rnaseq-results/
    star_salmon/
      <RNA_PATIENT>_tumor/
        quant.genes.sf
    TAA.xml
    TAA_all_genes.csv
    TAA/
      Targets_Gene_Sequence.txt
      protein_manifest.json
      TAA_Peptide.xml
      peptide_prediction/
        mixmhc2pred_input_15mers.txt
        mixmhc2pred_output.tsv
```

### Output file interpretation

| File | Content |
|---|---|
| `TAA.xml` | Expression-qualified gene targets from module 1 |
| `TAA_all_genes.csv` | Diagnostic table for evaluated genes and percentile filters |
| `Targets_Gene_Sequence.txt` | Retrieved reference protein sequences with identifiers |
| `protein_manifest.json` | Accession, Ensembl ID, reference sequence, length, sequence version, SHA256, and retrieval errors |
| `TAA_Peptide.xml` | Per-gene peptide predictions, HLA assignments, binding scores, and input-patient labels |
| `peptide_prediction/` | MixMHC2pred input/output artifacts |

The peptide XML contains `rna_patient` and `hla_patient` identifiers plus the HLA assignments used in prediction. These fields aid provenance tracking but do **not**, on their own, authenticate patient identity.

---

## Quality Control

| QC area | Check | Reason |
|---|---|---|
| Patient expression | Verify all selected TAA genes have **TPM >= 1** | Prevents all-zero reference ties from admitting unexpressed genes |
| Reference alignment | Validate GENCODE/GRCh38 and Ensembl gene IDs | Prevents gene-ID mismatches |
| Cohort selection | Record TCGA cancer type and GTEx tissue | Supports interpretation of expression percentiles |
| Sample provenance | Confirm RNA and HLA input files belong to the same patient | Prevents cross-patient peptide predictions |
| HLA compatibility | Verify every reported binding HLA is present in the verified genotype or an explicitly supported class II combination | Prevents inconsistent report-level HLA assignment |
| Protein mapping | Check UniProt gene, accession, Ensembl cross-reference, and protein length | Reduces gene-to-protein mismatches |
| Sequence provenance | Check cached SHA256 and UniProt sequence version | Detects corrupted or changed protein sequences |
| Peptide placement | Confirm each reported sequence occurs at its recorded 1-based protein start | Detects indexing or peptide-mapping errors |
| Run freshness | Confirm stage 2 uses stage 1 output from the intended run | Prevents stale TAA targets appearing in reports |
| Model coverage | Review unsupported HLA-II alleles | Makes prediction coverage limitations explicit |

For downstream report generation, a strict preflight should compare **stage 1 and stage 2 gene sets**, TPM values, patient identities, and HLA calls; fail rather than silently incorporate inconsistent results.

---

## Analytical Limitations

### Expression is not antigen presentation

Tumor RNA TPM and binding predictions cannot establish protein abundance, peptide processing, HLA surface presentation, T-cell recognition, vaccine effectiveness, or absence of normal-tissue toxicity.

### Canonical isoforms may differ from patient-expressed proteins

Reference protein 9-mers and 15-mers are derived from curated canonical sequences. Alternative splice isoforms, somatic variants, and variable-length regions—particularly in proteins such as MUC1—may alter the peptide repertoire present in an individual tumor.

### Repetitive and homologous genes

Some cancer-testis antigen loci (including CTAG-family genes) have closely related sequences that may complicate short-read RNA quantification and attribution to a specific genomic gene.

### Reference-cohort interpretation

TCGA and GTEx expression estimates arise from distinct projects and processing contexts. Cohort scaling and percentile ranking reduce some differences but do not completely remove cross-study batch effects. A high percentile is not equivalent to validated tumor specificity.

### HLA model and phase limitations

Not every typed HLA allele is necessarily supported by the installed predictors. DP/DQ chain pairing inferred combinatorially is not a substitute for experimentally phased or otherwise validated HLA haplotypes.

### No systematic normal-organ peptide safety validation

Current selection compares tumor RNA to one chosen GTEx tissue for the patient. It does not comprehensively exclude expression in all essential normal organs or validate off-tumor peptide presentation.

---

## Data Provenance and Reproducibility

A reproducible run should preserve:

- Study/sample identifiers in a protected environment (not in a public repository).
- Patient RNA-seq and HLA input-file paths and hashes.
- GENCODE, GTEx, GDC, and UniProt data versions.
- TCGA/GTEx cohort selection and reference sample counts.
- Patient TPM threshold and cohort-percentile cutoff values.
- MHCflurry and MixMHC2pred versions and model availability.
- Protein accessions, sequence versions, and SHA256 hashes.
- Input/output XML timestamps and hashes.
- Unsupported HLA models and unresolved protein mappings.
- Script commit hash, runtime configuration, and validation results.

No patient-identifiable report, raw FASTQ, genotype file, or clinical record should be published to GitHub without appropriate authorization and de-identification.

---

## Current Development Status

**As of October 26, 2026**

| Component | Status |
|---|---|
| 14-gene curated target panel | Implemented |
| STAR–Salmon TPM ingestion | Implemented (upstream quantification required) |
| GENCODE v50 gene mapping | Implemented |
| TCGA matched/pan-cancer percentile calculation | Implemented |
| GTEx selected-tissue percentile calculation | Implemented |
| Patient TPM >= 1 expression filter | Implemented in stage 1 |
| UniProt reviewed canonical protein retrieval and caching | Implemented |
| MHCflurry class I 9-mer prediction | Implemented |
| MixMHC2pred class II 15-mer prediction | Implemented |
| Top-10-per-gene peptide ranking | Implemented |
| TAA XML outputs and protein manifest | Implemented |
| Automated cross-patient identity and HLA consistency enforcement | Additional validation recommended |
| Prospective immunopeptidomics / T-cell validation | Not established by this repository |
| Portable public installation and regression-test suite | Further packaging and validation needed |

"Implemented" indicates available script logic, not independent clinical validation.

---

## Future Development

Potential enhancements include:

1. Require verified Data3–Data4 patient identity mapping and fail on mismatch.
2. Enforce cross-stage TPM, gene-set, and HLA consistency before report assembly.
3. Make reference paths, thresholds, and program locations configurable by a single YAML or CLI configuration.
4. Add regression tests with fully synthetic, non-identifiable datasets.
5. Review support for alternative transcript isoforms and patient-specific protein sequences.
6. Evaluate tumor/normal expression profiles across a broader safety-oriented tissue atlas.
7. Compare selected peptides with published experimentally observed HLA ligands, preserving evidence levels.
8. Test sensitivity to expression threshold, percentile method, reference scaling, and HLA model availability.
9. Document supported versions and reproducible environments for external researchers.
10. Evaluate candidates through orthogonal proteomic and immunological assays where appropriate.

The long-term objective is a transparent and reproducible research workflow for linking tumor-associated antigen expression with individualized peptide-binding predictions.

---

## Disclaimer

This project is intended for **research use only**. Its output comprises computationally prioritized peptide-binding candidates and must not be interpreted as confirmed tumor antigens, validated immunogenic epitopes, proven safe vaccine components, or clinical treatment recommendations.

All biological interpretation and translational use require appropriate experimental validation, clinical oversight, and applicable institutional review.

## License

No open-source license is asserted by this README. Repository owners should select a license only after confirming the applicable organizational, dependency, and reference-dataset permissions.
