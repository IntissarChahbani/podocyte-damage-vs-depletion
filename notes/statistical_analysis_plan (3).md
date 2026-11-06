# Statistical analysis plan (DRAFT, freeze together with PREREGISTRATION.md)

Items in [square brackets] are decisions still to be made. Fix them BEFORE any
analysis of the real data, then commit and tag. Anything changed afterwards goes
into the "Deviations" section at the bottom, with a date and a reason.

## 0. Principles
- Every primary result is reported with three things: an effect size, a 95%
  confidence (or credible) interval, and a p-value where one applies.
- Tests and thresholds are fixed in advance. Everything not listed here is
  labelled "exploratory".
- All analyses are run within a cohort first. Cohorts are never pooled before
  the replication step.
- Random seed fixed at [seed]; software versions recorded in environment.yml.

## 1. Data and units
- Unit of analysis: one glomerular biopsy sample.
- Cohorts: NEPTUNE (GSE254957, main); ERCB site H7 (GPL22945, replication);
  ERCB GPL24120 (secondary, descriptive); NEPTUNE tubulointerstitium
  (GSE219185, negative control).
- Group comparisons are made within one platform only.
- Pre-analysis checks: (a) duplicate or repeat biopsies from one person (if
  detectable); (b) library size / quality outliers. Exclusion rule, fixed now:
  PRIMARY analysis keeps all samples; the only exclusion criteria are technical
  failure (under 10 million reads, or not a glomerular biopsy). NEPTUNE
  GSE254957 has 282 samples with 22.5-67.5 million reads (median 36 million);
  none meets either criterion. SENSITIVITY analysis repeats the main analyses
  without the 8 samples whose log10 library size is more than 3 scaled MADs
  from the median (Sample_298, 244, 264, 336, 110, 270, 326, 179).
  Note (deviation log): a first rule (exclude beyond 3 MADs) was replaced by
  this one after viewing the depth distribution (not any score or group result)
  because 3 MADs is too strict for a tight, deep dataset.
- ERCB H7 (GPL22945, log2 scale, 71 samples): primary groups LD (n=18, controls),
  FSGS (10), MCD (5). TN (3) is a sensitivity control only. FSGS&MCD (6) is
  excluded from primary analyses (ambiguous label). RPGN (22) and DN (7) are
  exploratory. With MCD n=5, replication of H2 is judged by direction and size
  agreement with NEPTUNE, not by significance on its own.
- Age and sex are not available in the public metadata, so no adjustment is
  possible. This is stated as a limitation.

## 2. Score and abundance estimates
- Damage score (PDS), defined exactly as in Padvitski et al., Cell Genomics 2026
  (STAR Methods, "Damage score algorithm"): AUCell scores are computed
  separately for the genes that go UP (DS_up) and DOWN (DS_down) with damage,
  then combined as
      DS = (n_up / N) * DS_up  -  (n_down / N) * DS_down
  with n_up, n_down = numbers of up and down genes in the signature and N = 42
  for the full signature. HIGHER DS = MORE damage. Scores are scaled and
  centred within each cohort (the authors scale per study).
- Our adaptation (a documented deviation): for the common 34-gene score used in
  ERCB (32 down, 2 up) we use N = n_up + n_down = 34. Orthologs come from
  Ensembl BioMart (the authors used MGI); all 42 mapped one-to-one.
- Implementation: use the authors' code (Zenodo, function PDS_calc; DOIs
  10.5281/zenodo.19916603 and 10.5281/zenodo.20087796) for the score, or
  reimplement and validate against it on the same input. [Decide after
  inspecting the code.] The paper notes that their AUCell score is negative
  for low-ranking genes, which differs from standard AUCell, so do not assume
  that off-the-shelf AUCell in another package gives identical values.
- Direction check: because podocyte-identity genes (the "down" genes) rank
  lower when podocytes are fewer OR sicker, both depletion and damage push the
  bulk DS upwards. This is the confounding the project is about.
- FGFBP1 (a "down" gene) is detected (count > 0) in only 109 of 282 NEPTUNE
  samples; the other 41 genes in at least 240. It stays in the primary score;
  a sensitivity analysis removes it (part of the H5 grid).
- Gene sets: FULL = all usable one-to-one orthologs (primary). REDUCED = genes
  showing the expected direction in the other cohort (sensitivity only; list
  fixed before testing; selection and testing never use the same cohort).
- Abundance estimators (at least two): (A) marker-based, using podocyte markers
  that are NOT among the PDS genes; (B) reference-based non-negative least
  squares using a kidney single-nucleus reference, with PDS genes removed.
- Primary estimator: chosen from the SIMULATION benchmark by the lowest WORST-CASE
  mean absolute error of the podocyte fraction across the healthy-podocyte and
  damaged-podocyte scenarios. Chosen from simulated, not real, data.
  First benchmark run (2026-11-06, invented parameters: markers reduced to 70%
  in damaged podocytes): healthy MAE 0.007 (A marker) vs 0.023 (B NNLS); damaged
  bias -0.042 (A) vs -0.024 (B). Both under-count podocytes when they are damaged,
  which inflates the apparent depletion share. That fraction-error rule
  was replaced after the mediation simulations (errors in the final depletion
  share matter more): A (marker-based) is the primary estimator, B (NNLS) the
  sensitivity estimator. Decided from simulated, not real, results.
- Because under-counting in damaged tissue pushes the depletion share upward, the
  corrected share is read as an upper bound on depletion (see decision rule in H3).
- Abundance fractions are analysed on the logit scale (fractions are bounded
  and compositional).

## 3. Primary analyses

### H1: the score is HIGHER in disease than in healthy donors (higher = more damage)
- Omnibus test across groups within a cohort: [ANOVA if residual checks pass,
  otherwise Kruskal-Wallis]. Rule: residual normality (QQ plot and Shapiro-Wilk)
  and equal-variance check (Levene); the other test is reported as a check.
- Post-hoc: each disease group versus donors, Holm-corrected.
- Effect size: standardised mean difference (or Cliff's delta for rank tests)
  with a 10,000-resample bootstrap percentile 95% CI.

### H2: estimated podocyte abundance is lower in FSGS than in MCD
- Two-sided test (Welch t-test or Mann-Whitney, by the same rule as H1);
  predicted direction FSGS < MCD stated in advance.
- Effect size and bootstrap CI as above, on the logit scale.

### H3: a damage component remains after adjusting for abundance
- Model 1: score ~ group. Model 2: score ~ group + abundance.
- Mediation (statsmodels): total effect c, direct effect c', indirect effect
  a*b; proportion mediated = a*b / c. 95% bootstrap CI (10,000 resamples).
- Decision rule (fixed now, chosen on simulated data only): the proportion
  mediated (depletion share) is estimated with SIMEX correction and a stratified
  bootstrap 95% CI. "A damage component is supported" if the upper limit of that CI
  is below 1.0 and the total score difference between groups has a CI excluding 0.
  Otherwise "depletion alone is not rejected". No "mostly depletion" or "mostly
  damage" claim is made. The corrected share is an UPPER bound on depletion,
  because damaged podocytes lose marker signal and look like lost podocytes.
  Simulation (30 per group, 20 datasets per scenario): false-alarm rate with
  depletion only 5% (marker estimator), power 45% in a mixed scenario and 100% when
  damage dominates. The false-alarm estimate is uncertain (20 datasets).
- Measurement-error correction: abundance is an estimate, and adjusting for a
  noisy estimate under-corrects. We repeat the analysis with simulation-
  extrapolation (SIMEX), using the error variance measured in the simulation
  benchmark. Corrected and uncorrected results are both reported. Primary
  conclusion: the SIMEX-corrected one (noise SD taken from the benchmark; the
  simulation showed the uncorrected share is attenuated, e.g. 0.47 vs truth 1.02
  for the NNLS estimator).
- Extra (exploratory): logistic regression, group (MCD vs FSGS) ~ abundance +
  damage component, reported as odds ratios with 95% CI, compared with a model
  containing abundance alone (likelihood-ratio test). Described as
  descriptive, not as clinical prediction.

### H4: replication
- Same analysis in ERCB site H7 for the shared groups (FSGS, MCD, donors).
- Criteria (fixed now): same direction of every pre-specified effect, and the
  ERCB estimate inside the NEPTUNE 95% CI, or the difference between the two
  estimates with CI including zero. [Adjust if needed.]
- Permutation null: group labels shuffled within each cohort 10,000 times;
  the observed statistic is compared with the null distribution.
- Bayesian hierarchical model (PyMC or bambi): group effect with a cohort-level
  random effect (partial pooling). Priors: weakly informative, [Normal(0, 1) on
  standardised effects], stated and checked with prior-predictive simulation.
  Reported: posterior median and 95% credible interval per cohort and pooled.

### H5: robustness
- Grid of configurations: 2 estimators x FULL / REDUCED gene set x PDS genes
  removed from the markers (yes/no). The H3 decision (damage supported /
  depletion alone not rejected) is reported for every cell. Criterion: the
  conclusion is "robust" if the decision is the same in at least 80% of
  configurations; otherwise it is "not robust"
  and reported as such.

### H6: negative control (tubulointerstitium)
- Equivalence test (TOST): the group effect on the score and on abundance must
  lie within +/- [0.3] standard deviations. Margin fixed now.
- "Not significant" is never presented as evidence of no effect.

## 4. Multiple testing
- Primary family: H1 post-hoc comparisons, H2, H3 (indirect effect), H4. Holm
  correction within the family.
- Everything else (logistic model, extra groups, GPL24120 cohort, REDUCED gene
  set) is exploratory and labelled as such; no formal significance claims.
- ERCB GPL24120: group effects with centre as a random effect (mixed model),
  descriptive only.

## 5. Power and detectable effect
- Before analysis: simulation-based power analysis for the smallest group
  (8 donors in NEPTUNE, 18 in ERCB H7). Reported: the minimum detectable
  standardised effect at 80% power. Results go in the README.

## 6. Validation on simulated data
- Every method (score, estimators, mediation, SIMEX, Bayesian model) is first run on
  artificial mixtures with known depletion and known damage. A method is used on
  real data only if it recovers the known fraction in the healthy scenario with
  absolute bias <= 0.02 and MAE <= 0.03. The worst-case bias (damaged scenario)
  is reported with every abundance-based conclusion. (An earlier tolerance of 10
  percentage points was dropped as too lax for fractions between 0.03 and 0.25.)

## 7. Limitations stated in advance
- No kidney function or outcome data; no age or sex adjustment.
- Abundance is estimated, not measured; podocytes are rare in references.
- Group comparisons support a pattern, not a per-patient decomposition.
- ERCB platform and centre are confounded; GPL24120 has only 5 controls.
- Mouse-derived score: the number of usable human orthologs is reported.

## 8. Deviations log
(Date | what changed | why | effect on conclusions)
