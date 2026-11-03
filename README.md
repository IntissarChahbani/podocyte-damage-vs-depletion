# Podocyte Damage vs Depletion

**Does a low podocyte damage score in a kidney biopsy mean the
podocytes are damaged, or that there are fewer of them?**

Status: planning and pre-registration. No results yet.

## Background
The Podocyte Damage Score (PDS; Padvitski et al., Cell Genomics 2026)
is a 42-gene score derived from mouse models. In bulk tissue a low
score could reflect fewer podocytes (depletion) or sicker podocytes
(damage). The authors did not assess prognosis.

## Objectives
1. Score public human glomerular transcriptomes with the PDS
   (human orthologs).
2. Estimate podocyte abundance per biopsy with at least two methods.
3. Quantify how much of the PDS-kidney function association is
   explained by abundance versus a residual damage component.
4. Test whether the result replicates across cohorts and is robust
   to the choice of estimator.

## Data (public, to be verified)
ERCB glomerular (GSE104948); NEPTUNE glomerular (GSE219185,
GSE197307); human kidney single-nucleus reference (KPMP/HuBMAP).

## Approach
Python (GEOparse, scanpy, decoupler, statsmodels), Snakemake workflow,
pre-registered hypotheses (H1-H4), patient-level permutation null.
Optional R-based estimator (MuSiC) as a robustness check.

## Status
- [x] Repository set up
- [ ] Scoop check logged
- [ ] PDS gene list and ortholog mapping
- [ ] Metadata audit
- [ ] Pre-registration committed