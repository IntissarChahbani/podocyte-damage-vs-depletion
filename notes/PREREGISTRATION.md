# Pre-registration v1: podocyte damage vs depletion in bulk glomerular transcriptomes

Frozen before any Podocyte Damage Score (PDS) was computed on real data.
Tag: `prereg-v1`. Companion file: `notes/statistical_analysis_plan.md` (full methods).
Anything changed after the tag goes into the Deviations log with date and reason.

## 1. Question
In bulk glomerular transcriptomes, a HIGHER PDS means more podocyte damage
(Padvitski et al., Cell Genomics 2026). But fewer podocytes (depletion) also push
the bulk score up. How much of the score difference between disease and healthy
tissue is explained by estimated podocyte abundance, and how much remains as a
damage component?

## 2. Data (all public, no access request)
- Main: NEPTUNE glomerular RNA-seq, GSE254957 (282 samples, raw counts).
- Replication: ERCB site H7, GSE104948 / GPL22945 (71 samples, log2 microarray).
- Secondary, descriptive: ERCB GPL24120.
- Negative control: NEPTUNE tubulointerstitium, GSE219185.
- Abundance reference: public kidney single-nucleus atlas (access to be verified).
- Phase 1 is diagnosis-based. No kidney function or outcome analysis.

## 3. Score
- 42 mouse-derived genes (39 down, 3 up with damage), mapped to human with Ensembl
  BioMart (all 42 one-to-one). NEPTUNE uses 42/42; ERCB and cross-cohort
  replication use the 34 common genes (32 down, 2 up, N = 34).
- DS = (n_up/N) * DS_up - (n_down/N) * DS_down, AUCell-based, scaled within cohort.
- Implementation: the authors' code (Zenodo 10.5281/zenodo.19916603 and
  10.5281/zenodo.20087796) if usable; otherwise a reimplementation that must
  reproduce a validation input before use on real data.
- FGFBP1 (down gene, detected in 109/282 samples) stays in the primary score;
  sensitivity analysis removes it.

## 4. Samples
- NEPTUNE primary analysis: all 282 samples. Sensitivity analysis excludes the 8
  library-size outliers (Sample_298, 244, 264, 336, 110, 270, 326, 179).
- ERCB H7 primary groups: LD (18, controls), FSGS (10), MCD (5). TN (3) sensitivity
  control only; FSGS&MCD (6) excluded from primary; RPGN (22) and DN (7) exploratory.
- Assumption to verify in the literature before interpretation: podocytes are
  largely preserved in minimal change disease and lost in FSGS.

## 5. Hypotheses (directions fixed in advance)
- H1: PDS is higher in disease groups than in healthy donors.
- H2: estimated podocyte abundance (logit fraction) is lower in FSGS than in MCD.
- H3: a damage component remains after adjusting PDS for abundance.
- H4: H1-H3 replicate in ERCB H7 in direction and size (permutation null).
- H5: the H3 conclusion is stable across estimators, gene sets (FULL / REDUCED /
  without FGFBP1) and marker handling.
- H6: in tubulointerstitium, group effects on PDS and abundance lie within
  +/- 0.3 SD (equivalence test, TOST).

## 6. Fixed decisions
- Random seed: 20261106.
- Omnibus test: ANOVA if residual normality (QQ plot, Shapiro-Wilk) and equal
  variance (Levene) hold, otherwise Kruskal-Wallis; the other is a check.
- Post-hoc and primary family corrected with Holm. Primary family: H1 post-hoc,
  H2, H3 indirect effect, H4. Everything else is exploratory.
- Effect sizes with 10,000-resample bootstrap 95% CIs.
- H3 classification of the proportion mediated (a*b/c): "mostly depletion" if the
  whole CI >= 50%; "mostly damage" if the whole CI <= 25%; "mixed" if the whole
  CI lies in 25-50%; otherwise "inconclusive". The ratio is interpreted only if the
  total-effect CI excludes 0.
- Primary conclusion requires the uncorrected and SIMEX-corrected analyses to give
  the same classification; if they differ, the result is "inconclusive".
- H5: "robust" if the classification is identical in at least 80% of configurations.
- H4 replication: same direction for every pre-specified effect, and the ERCB
  estimate inside the NEPTUNE 95% CI (or the between-cohort difference CI includes 0).
- Asymmetric rule: "mostly damage" is claimed only if it survives correction and
  robustness; absence of evidence is never reported as evidence of absence.
- Method validation: each method is first run on simulated mixtures with known
  depletion and damage and is used on real data only if it recovers the known
  fraction within 10 percentage points. Primary abundance estimator is chosen from
  simulation (lowest mean absolute error), not from real data.
- Power analysis by simulation, before analysis, for the smallest groups.

## 7. Limitations stated in advance
- No outcome or kidney function data; no age or sex adjustment.
- Abundance is estimated, not measured.
- Group-level pattern, not per-patient decomposition.
- ERCB MCD n = 5: replication of H2 is judged by direction and size.
- Mouse-derived score; platform and centre are confounded in ERCB.
- Literature check status: scoop check not complete (citation lookups pending at
  the time of freezing); repeat before write-up.

## 8. Deviations log
| Date | What changed | Why | Effect |
|---|---|---|---|
| 2026-11-06 | Library-size rule changed from "exclude beyond 3 MAD" to "keep all; 8 outliers in sensitivity" | 3 MAD too strict for deep, tight data; decided from depth distribution only, before any score or group result | none expected |
