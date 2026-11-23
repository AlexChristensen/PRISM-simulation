# PRISM-simulation

**Recovering Hierarchical Structures in Psychological Measurement**

Simulation code, results, and figures for the PsyArXiv preprint introducing **PRISM** (*Pairwise Resolution Iteration via Subgraph Modularity*), a community detection procedure for recovering every level of a hierarchical structure without a tunable resolution parameter.

> Christensen, A. P., Goldberg, S. M., Watts, A. L., & Garrido, L. E. (2026). *Recovering hierarchical structures in psychological measurement*. PsyArXiv preprint (under review). The full manuscript is included in this repository as [`PP_PRISM.pdf`](PP_PRISM.pdf).

---

## Table of contents

- [Overview](#overview)
- [Summary of the manuscript](#summary-of-the-manuscript)
  - [1. Introduction](#1-introduction)
  - [2. Present research](#2-present-research)
  - [3. Method](#3-method)
  - [4. Results](#4-results)
  - [5. Discussion](#5-discussion)
- [Key findings at a glance](#key-findings-at-a-glance)
- [Repository structure](#repository-structure)
- [Reproducing the simulation](#reproducing-the-simulation)
- [Software](#software)
- [Citation](#citation)

---

## Overview

Psychological constructs are often organized hierarchically (e.g., the Big Five in personality, HiTOP in psychopathology), so accurately identifying the *narrowest* dimensions is essential to recovering the broader levels above them. Existing network approaches are limited either by the **resolution limit of modularity** (the standard Louvain algorithm merges small communities) or by the **variability introduced by stopping Louvain early** (lower order Louvain).

PRISM addresses both problems:

1. Run standard Louvain on the full network to get an initial partition.
2. For each community, build induced subgraphs: the community on its own, and its union with every other community (one at a time).
3. Re-run Louvain within each subgraph. A community is updated only if it shows a **modularity gain in every subgraph it participates in**. The new partition is the one most frequently observed across those subgraphs.
4. Iterate until no partition changes.

Because an induced subgraph shrinks the effective network size (the `2m` term in the modularity null model), small communities become detectable again, with no resolution parameter to tune. If no substructure is found, PRISM returns the standard Louvain solution unchanged, so it works for non-hierarchical structures too.

The simulation compares PRISM against hierarchical exploratory graph analysis (hierEGA) with lower order Louvain and against exploratory factor analysis (EFA) paired with several traditional dimensionality methods, across **6,400 conditions × 100 replicates = 640,000 samples** (160 structural conditions × 5 sample sizes × 2 skew levels × 2 overlap × 2 error; verified from the saved results, 5,120,000 rows = 640,000 samples × 8 methods).

The PRISM algorithm itself is implemented as `community.prism()` in the [`{EGAnet}`](https://github.com/hfgolino/EGAnet) R package (v2.5.0). This repository wraps it (see `source/methods.R`) and contains everything needed to reproduce the simulation study.

---

## Summary of the manuscript

### 1. Introduction

- **Motivation.** Hierarchical frameworks are everywhere in psychology, yet simulation studies of exploratory methods for hierarchical structure are rare. Only one prior study (Jiménez et al., 2025) had evaluated them, and it used normally distributed continuous data with few overlapping dimensions.
- **Contemporary methods (§1.1).**
  - *Bass-ackward* (Goldberg, 2006) is a top-down factor-analytic approach. Its depth is fixed in advance by a traditional dimensionality estimate, so it cannot recover levels beyond that depth. The extended version (Forbes, 2023) adds redundancy removal but does not increase depth.
  - *hierEGA* is a bottom-up network approach: estimate an EBICglasso network, detect communities with Louvain, compute network scores, and repeat at the next level.
- **Resolution limit (§1.2).** Two communities are merged by modularity maximization whenever `W_AB > k_A·k_B / 2m`, so the tolerable between-community weight shrinks as the network grows. Lower order Louvain sidesteps this by stopping after one pass, at the cost of order-dependent variability. Existing fixes (e.g., a resolution parameter ν) add researcher degrees of freedom.
- **PRISM (§1.3).** Builds on the hierarchical multiresolution idea of Granell et al. (2012) but replaces the tunable parameter with a unanimous-gain rule over induced subgraphs, which guards against oversplitting.
- **Demonstrations.** On two classic synthetic networks, PRISM recovers the true structure where the alternatives fail:
  - *Fortunato & Barthélemy ring network* (10 true communities): standard Louvain finds 6; lower order Louvain and PRISM find 10.
  - *Arenas et al. hierarchical network* (16 first-level communities): standard Louvain finds 4–6; lower order Louvain returns 15–22 communities depending on the seed; **PRISM recovered all 16 in 100 of 100 seeds**.

### 2. Present research

The study extends Jiménez et al. (2025) in four ways: (1) 5-point ordinal data with moderate and extreme skew, (2) substantial cross-loadings (0.40–0.50) creating overlapping dimensions, (3) evaluation of *structure* recovery and not only the number of dimensions, and (4) population error to emulate model misspecification. A single second order factor condition tests whether PRISM can also recover a non-hierarchical (correlated factors) structure.

Methods compared:

| Level | Methods |
|---|---|
| First **and** second order | PRISM, hierEGA with lower order Louvain, parallel analysis with principal axis factoring (PApaf), parallel analysis with principal components (PApca) |
| First order only (bass-ackward depth setting) | Goldberg's maximum loading criterion, Velicer's MAP, VSS1, VSS2 |

### 3. Method

- **Design (§3.1).** Data are generated from second order factor models and varied on: second order factors (1, 2, 4, 6), first order factors per second order factor (3, 5), indicators per first order factor (5, 10), first/second order loadings (0.50, 0.70), second order correlations (0.10, 0.30, 0.50), sample size (500, 1,000, 2,500, 5,000, 10,000), skew (moderate, extreme), overlapping dimensions (absent, present), and population error (absent, present). The manuscript reports 3,200 conditions × 100 replicates = 320,000 samples, but the saved results and simulation code contain 6,400 conditions × 100 replicates = 640,000 samples (the manuscript's count appears to omit the 5-vs-10 indicators factor).
- **Data generation (§3.2).**
  - *Population model:* `R_R = Λ_F Φ_F Λ_F′ + Ψ_F` with `Φ_F = Λ_S Φ_S Λ_S′ + Ψ_S`.
  - *Overlapping dimensions:* cyclical cross-loadings (0.40–0.50) on 20% of items per factor.
  - *Population error:* Tucker–Koopman–Linn minor factors (`{fungible}`, q = 50), calibrated to **CFI = 0.97**, chosen through a preliminary simulation as the most misfit that still allowed recovery of the structure.
  - *Convergence criteria:* positive-definite matrices, communalities ≤ 0.90, target misfit, residual correlations ≤ 0.10, and recoverable population structure.
  - *Skew:* moderate skew drawn from U(1, 2). Extreme skew drawn from a gamma distribution (shape 1.70, rate 1) fit to item skew in the SCL-90-R and HiTOP-SR (mean 1.70, range −0.05 to 12.55).
- **Estimation (§3.3).** All methods use polychoric correlations (with empty-cell correction). hierEGA uses EBICglasso networks with consensus clustering (1,000 Louvain runs). Second order structure is estimated from composite correlations of network scores (PRISM, Louvain) or factor scores (parallel analysis). First-order-only methods (max loading, MAP, VSS) are implemented in-house following `{psych}`.
- **Evaluation metrics (§3.4).** Exact dimensionality **accuracy**, **mean bias error** (MBE; positive = overfactoring), **mean absolute error** (MAE), and the **Omega index** (ω; Collins & Dent, 1988), a generalization of the adjusted Rand index that handles overlapping dimensions. ω > 0.65 / 0.80 / 0.90 is read as moderate / good / excellent recovery.
- **Transparency (§3.5).** R 4.5.1, `{latentFactoR}` 0.0.9, `{EGAnet}` 2.5.0, `{psych}` 2.6.5, `{fungible}` 2.4.8. All code is in this repository.

### 4. Results

- **Overall first order recovery (Table 1).** PRISM matched or beat every method on nearly every metric. The one exception: with overlap and moderate skew, PApaf had slightly higher accuracy when population error was absent (87.9% vs. 85.0%), though ω was comparable and PRISM had smaller MAE. Moving from moderate to extreme skew cost PRISM less than 9% accuracy, while overlapping dimensions cost about 8–14%.
- **Traditional depth-setting heuristics performed poorly.** Max loading, MAP, VSS1, and VSS2 never exceeded 75% accuracy and underfactored by more than two factors on average (MBE −1.2 to −8.0).
- **PApaf** was strong under moderate skew but very sensitive to extreme skew (accuracy dropped about 43–47%, with large overfactoring). **PApca** was the best traditional alternative and less affected by overlap, but needed very large samples.
- **Figure 3 (first order).** With a single second order factor (a correlated factors model), PRISM had near-zero bias (MBE = 0.051), whereas lower order Louvain tended to oversplit (MBE = 0.249). In hierarchical conditions PRISM reached about 90% accuracy or better once N ≥ 2,500. PApca matched or slightly exceeded PRISM only at N = 10,000.
- **Figure 4 (second order).** PRISM generally reached at least 80% accuracy at every sample size. PApca degraded as the number of second order factors grew, and PApaf never exceeded 75%.
- **Figure 5 (sensitivity).** PRISM was the least affected method across every manipulated variable (first order accuracy range 0.3–21.5 points, second order 0.1–9.2). Sample size was the largest driver for the other methods, and skew was the largest for PApaf.
- **Figure 6 (hardest condition).** All methods struggled when *overlapping dimensions*, *few variables per factor (5)*, *low first order loadings (0.50)*, and *high second order loadings (0.70)* co-occurred. Each factor alone, and every pair, was recoverable by PRISM and PApca at N = 10,000.

### 5. Discussion

- **PRISM as a general procedure (§5.1).** PRISM recovers both levels accurately, is comparatively robust to skew, and works at smaller sample sizes than competing methods. At N = 2,500 it often did as well as other methods did at much larger N. It rarely oversplits when no substructure exists, which supports its use for hierarchical and non-hierarchical structures alike.
- **Bass-ackward depth (§5.2).** Most traditional dimensionality methods underfactor, especially under extreme skew and overlap, which are both common in psychopathology data. Because bass-ackward cannot recover levels deeper than its initial depth, and its refinements cannot fix underfactoring, the authors advise caution when interpreting existing bass-ackward hierarchies. PApca was the most viable traditional option.
- **Applied implications (§5.3).** PRISM does not require a commitment to network *theory*. Network methods are used here purely as a tool for structure recovery. Applied recursively, PRISM lets the depth and organization of a hierarchy emerge bottom-up from one procedure.
- **Limitations (§5.4).**
  - Data came only from two-level second order factor models. Deeper hierarchies, unequal indicators, and other generating models (e.g., hierarchical stochastic block models) remain to be tested.
  - Only polychoric correlations were evaluated, not Spearman.
  - The max loading search range contained the true number of factors in 97.1% of replicates.
- **Conclusion (§5.5).** PRISM had the most accurate recovery of first and second order structures across a diverse set of conditions, including sample sizes where traditional methods were inaccurate and biased.

---

## Key findings at a glance

First order recovery by method, with population error absent. Values are exact dimensionality accuracy / Omega index (ω). Full results, including error-present conditions, MBE, and MAE, are in Table 1 of the manuscript.

| Method | No overlap, moderate skew | No overlap, extreme skew | Overlap, moderate skew | Overlap, extreme skew |
|---|---|---|---|---|
| **PRISM** | **0.988 / 0.996** | **0.903 / 0.984** | 0.850 / 0.952 | **0.822 / 0.944** |
| Lower order Louvain | 0.923 / 0.992 | 0.861 / 0.978 | 0.739 / 0.942 | 0.708 / 0.931 |
| Parallel analysis (PCA) | 0.842 / 0.950 | 0.756 / 0.932 | 0.788 / 0.938 | 0.750 / 0.923 |
| Parallel analysis (PAF) | 0.859 / 0.973 | 0.394 / 0.937 | **0.879 / 0.973** | 0.450 / 0.943 |
| Max loading | 0.741 / 0.878 | 0.323 / 0.854 | 0.683 / 0.833 | 0.330 / 0.815 |
| MAP | 0.668 / 0.789 | 0.620 / 0.760 | 0.697 / 0.839 | 0.645 / 0.809 |
| VSS1 | 0.439 / 0.659 | 0.364 / 0.645 | 0.216 / 0.437 | 0.160 / 0.428 |
| VSS2 | 0.424 / 0.759 | 0.350 / 0.737 | 0.366 / 0.680 | 0.315 / 0.667 |

> Bold marks the best accuracy and/or ω in each column. The one place PRISM is not best is overlap with moderate skew, where PApaf is slightly ahead on both accuracy and ω (the manuscript calls the ω difference trivial and notes PRISM's smaller MAE).

**Takeaways**

- PRISM is the most accurate and least biased method in nearly every condition, at both levels of the hierarchy.
- It is the most robust to extreme skew and to small samples (N ≤ 1,000).
- Overlapping dimensions are the main weakness for the network methods, and the combination of overlap, few indicators, weak first order loadings, and strong second order loadings is hard for every method.
- PApca is the strongest traditional alternative but only catches up at N = 10,000.
- Max loading, MAP, VSS1, and VSS2 consistently underfactor and are poor choices for setting bass-ackward depth.

Figures from the manuscript are in [`simulation/figures/`](simulation/figures) (`figure3.pdf` to `figure6.pdf`).

---

## Repository structure

```
PRISM-simulation/
├── README.md
├── PP_PRISM.pdf                        # Preprint manuscript
├── .gitignore
└── simulation/
    ├── 00_simulate_error.R             # Preliminary study: for target CFI 0.95–0.99, add population error to
    │                                   #   each condition and check the structure stays recoverable
    │                                   #   (this is how CFI = 0.97 was chosen)
    ├── 01_simulate_conditions.R        # Generate 100 population-model replicates per structural condition
    │                                   #   (no overlap / overlap × no error / error)
    ├── 02_data_generation.R            # Draw 5-point ordinal samples from each population model at every N,
    │                                   #   for both moderate and extreme skew
    ├── 03_simulation.R                 # Apply all 8 methods to every sample; save results and timing
    ├── 04_tables-and-figures.R         # Table 1 and Figures 3–6; also checks the max-loading search range
    ├── 05_timing-analysis.R            # Collate per-run computation times into timing_df.RData
    ├── PRISM_talk.R                    # Computes item skew of the generated extreme-skew data → skew_df.RData
    │
    ├── source/                         # Functions sourced by the numbered scripts
    │   ├── generate_error.R            #   Used by 00: population error at each target CFI, recoverability check
    │   ├── generate_condition.R        #   Used by 01: builds the four population structures per replicate
    │   │                               #   (valid loadings, overlap, TKL error at CFI = 0.97, recoverability)
    │   ├── generate_no_error.R         #   Earlier version of the generator (no population error, 0.30–0.50
    │   │                               #   cross-loadings); not sourced by any current script
    │   ├── generate_data.R             #   Used by 02: multivariate normal draws, then 5-category skew
    │   ├── methods.R                   #   Used by 00, 01, 03: polychoric correlations, parallel analysis
    │   │                               #   (PAF/PCA), MAP, VSS1/2, max loading, EFA, network loadings,
    │   │                               #   lower order Louvain, and PRISM (wraps EGAnet::community.prism)
    │   ├── obtain_results.R            #   Used by 03: run and time all methods on one sample
    │   └── results.R                   #   Used by 00, 01, 03: membership assignment, Omega index,
    │                                   #   per-method first and second order dimensions
    │
    ├── cross_loadings/
    │   └── 00_cross-loadings.R         # Supplementary simulation: network-loading-based factor correlations
    │                                   #   with simple-structure vs. full loadings (see below)
    │
    ├── figures/
    │   ├── figure3.pdf                 # First order accuracy and structural recovery
    │   ├── figure4.pdf                 # Second order accuracy and structural recovery
    │   ├── figure5.pdf                 # Worst vs. best accuracy by condition
    │   └── figure6.pdf                 # Most challenging condition combination
    │
    ├── results.RData                   # Final results (data frame, 5,120,000 rows × 21 columns;
    │                                   #   one row per sample × method); loaded by 04
    ├── results_all.RData               # Same layout, from the earlier run without the empty-cell correction
    ├── timing.RData                    # Per-sample run times (640,000 rows) for the final run
    ├── timing_all.RData                # Run times from the earlier run without the empty-cell correction
    └── skew_df.RData                   # Item skew values (1,776,000 rows) from PRISM_talk.R
```

Directories created when running the pipeline (git-ignored, so not in the repository): `simulation/error/`, `simulation/conditions/`, `simulation/data/{moderate,extreme}/`, `simulation/results/`, and `simulation/timing/`.

**Notes on specific files**

- `results.RData` columns: `number` (structural condition, 1–160), `LowV`, `LowF`, `LowL`, `HighF`, `HighL`, `HighC`, `N`, `Error`, `Overlap`, `SKEW`, `Replicate`, `METHOD` (`prism`, `lower_louvain`, `pca`, `paf`, `map`, `vss1`, `vss2`, `max`), estimated `lower_dimensions` / `higher_dimensions`, and structure recovery plus accuracy and bias flags (`LOWER_CORRECT`, `LOWER_MBE`, `HIGHER_CORRECT`, `HIGHER_MBE`). The columns named `lower_ari` / `higher_ari` hold the **Omega index**, which equals the adjusted Rand index when no dimensions overlap.
- `results_all.RData` / `timing_all.RData` come from the repository commit "updated results with no correction" (before `polychoric()` gained the empty-cell correction of +0.50 for 2 × 2 tables). They differ from the final files in roughly a fifth of the rows. Only `results.RData` is used for the manuscript's tables and figures.
- `cross_loadings/00_cross-loadings.R` simulates correlated factors with cross-loadings (4/8/12 variables per factor, 2/4/6 factors, loadings 0.50–0.70, cross-loadings 0.05–0.15, N = 1,000), estimates EGA network loadings, and compares factor correlations computed from simple-structure loadings against the full loading matrix. Results are pasted as comments at the bottom of the script: simple-structure loadings had lower mean absolute error (paired *t* = −6.98, Cohen's *d* ≈ 0.78), which is consistent with the paper using simple-structure network loadings to build the composite correlations.

---

## Reproducing the simulation

Run the scripts in order from inside the `simulation/` directory (they use relative paths such as `./source/` and `./conditions/`):

| Step | Script | Produces |
|---|---|---|
| 0 | `00_simulate_error.R` | Preliminary population-error calibration (`error/`); prints recoverability, CFI, and residual summaries by target CFI |
| 1 | `01_simulate_conditions.R` | Four population structures per replicate (no overlap / overlap × no error / error), 100 replicates per condition (`conditions/`) |
| 2 | `02_data_generation.R` | Ordinal samples at every N for moderate and extreme skew (`data/moderate/`, `data/extreme/`) |
| 3 | `03_simulation.R` | Per-sample results and timing (`results/`, `timing/`), combined into `results.RData` and `timing.RData` |
| 4 | `04_tables-and-figures.R` | Table 1 and Figures 3–6 |
| 5 | `05_timing-analysis.R` | Aggregated timing (`timing_df.RData`) |

**Notes**

- Steps 0–3 are computationally heavy. The scripts parallelize with `parallel::mclapply` (often 50 cores), pin BLAS/OpenMP to a single thread, and call `sudo` to drop the page cache between runs. Edit the core counts and remove or adapt the `sudo` lines for your machine.
- To regenerate only the tables and figures, you do not need to rerun the simulation: `results.RData` is included, so you can run `04_tables-and-figures.R` directly.

---

## Software

R 4.5.1 with:

| Package | Version | Used for |
|---|---|---|
| `latentFactoR` | 0.0.9 | Data generation, skew tables |
| `EGAnet` | 2.5.0 | EGA networks, PRISM (`community.prism`), consensus Louvain, network loadings |
| `fungible` | 2.4.8 | Population error (`noisemaker`) |
| `psych` | 2.6.5 | Reference implementations for MAP, VSS, parallel analysis |
| `ggplot2`, `ggpubr`, `tidyverse` | 4.0.3, 0.6.3 | Tables and figures |
| `RhpcBLASctl` | | Thread control |

---

## Citation

If you use PRISM or this simulation, please cite the preprint:

```
Christensen, A. P., Goldberg, S. M., Watts, A. L., & Garrido, L. E. (2026).
Recovering hierarchical structures in psychological measurement. PsyArXiv preprint.
```

The manuscript has not been peer-reviewed and is currently under review.
