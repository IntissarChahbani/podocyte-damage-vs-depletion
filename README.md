# Podocyte Damage vs Depletion

**Does a low podocyte damage score in a kidney biopsy mean the
podocytes are damaged, or that there are fewer of them? And which
one tracks kidney function?**

Status: planning and pre-registration. No results yet.

## Background
The Podocyte Damage Score (PDS; Padvitski et al., Cell Genomics 2026,
doi:10.1016/j.xgen.2026.101337) is a 42-gene score derived from mouse
models. In bulk tissue a low score could reflect fewer podocytes
(depletion) or sicker podocytes (damage). This distinction is my
inference from how the score is built, and it is what this project
tests. The authors did not assess prognosis. They report that the
score separates controls from disease in bulk data, so that finding
is treated here as a quality check, not a new result.

## Objectives
1. Score public human glomerular transcriptomes with the PDS
   (human orthologs).
2. Estimate podocyte abundance per biopsy with at least two methods.
3. Quantify how much of the PDS-kidney function association is
   explained by abundance versus a residual damage component.
4. Test whether the result replicates across cohorts and is robust
   to the choice of estimator.

Phase 2 (only if outcome data become available): does the damage
component predict kidney failure or eGFR decline beyond abundance
and clinical covariates?

## Data (public, to be verified)
ERCB glomerular (GSE104948); NEPTUNE glomerular (GSE219185,
GSE197307); human kidney single-nucleus reference (KPMP/HuBMAP).
Available clinical variables (eGFR, proteinuria, disease group)
will be confirmed in a metadata audit before any analysis is planned
around them.

## Approach
Python (GEOparse, scanpy, decoupler, statsmodels), Snakemake workflow,
pre-registered hypotheses (H1-H4), patient-level permutation null.
Optional R-based estimator (MuSiC) as a robustness check.

| ID | Hypothesis |
|----|-----------|
| H1 | Bulk PDS is associated with eGFR and proteinuria (quality check) |
| H2 | A damage component remains after adjusting for estimated podocyte abundance |
| H3 | The decomposition replicates across cohorts against a patient-level permutation null |
| H4 | The conclusion holds with at least two abundance estimators and with PDS genes removed from the abundance markers |

Thresholds are fixed in `PREREGISTRATION.md` before any analysis.

## Safeguards
- **Circularity:** abundance markers must not overlap the 42 PDS genes.
- **Rare cell type:** podocytes are poorly captured in single-cell
  references; use more than one estimator, benchmark on pseudobulk
  with known proportions, report uncertainty.
- **Reference mismatch:** healthy references may not fit diseased tissue.
- **Mouse-derived score:** report how many of the 42 genes have usable
  human orthologs.
- **Platform differences:** score within each cohort, compare effect
  sizes, not raw values.

## Limits
- Public bulk cohorts have no direct podocyte counts, so abundance is
  estimated, not measured.
- Damage and depletion may be partly entangled (the score correlates
  with podocyte density in spatial data); an inconclusive split is
  reported as such.
- Phase 2 depends on outcome data that may need an application.
- The project measures how much of a score reflects abundance versus
  damage; it does not establish mechanisms.

## Possible later extensions
Pathway analysis of the residual damage component; threshold or
non-linear relationship between score and kidney function.

## Status
- [x] Repository set up
- [ ] Scoop check: web searches logged; citation checks pending
      (see `notes/scoop_check.md`)
- [ ] PDS gene list and ortholog mapping
- [ ] Metadata audit
- [ ] Pre-registration committed

## Repository structure
```
README.md
PREREGISTRATION.md
notes/            scoop_check.md, log.md
data/             manifest.csv, download script (no raw data committed)
workflow/         Snakefile
envs/             environment.yml
results/  report/
```

## Key references
- Padvitski et al., Cell Genomics 2026: https://doi.org/10.1016/j.xgen.2026.101337
- Preprint: https://doi.org/10.1101/2025.02.20.639269