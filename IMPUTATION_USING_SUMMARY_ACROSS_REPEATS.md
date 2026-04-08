# Imputation using summary statistics across repeats

This note documents an alternative way of using the 200 repeated cross-fits
in `repeated_cv_proxy_crf_asinh.RData`. Instead of computing performance
metrics within each repeat and averaging, we collapse the 200 firm-year
proxies into two summary statistics — selection frequency `p_i` and mean
proxy `proxy_mean` — and build the prediction pipeline directly on those.

The note records the conceptual choices, the diagnostics that motivated each
step, and the final numbers. Each result is linked to the script that
produced it. Output objects live in
[`analysis/active/output/`](analysis/active/output/).

---

## 1. Setup and motivation

[`build_repeated_cv_proxy_crf.R`](analysis/active/build_repeated_cv_proxy_crf.R)
produces, for every firm-year `i` in the training panel, 200 held-out elastic
net proxy values `proxy_{i,1}, ..., proxy_{i,200}` — one per repeat of the
CRF-group K-fold cross-validation. Panel B of
[`figures_tables/table_main_with_mixed_crf.R`](figures_tables/table_main_with_mixed_crf.R)
uses these per-repeat proxies as ranking signals: in each repeat it
redistributes (CRF group, year) totals proportionally, then averages
performance metrics across the 200 repeats. The 200 repeats are used purely
as a variance-reduction device for the metric.

The alternative explored here is to collapse the 200 repeats into per
firm-year summaries:

- `p_i = mean(proxy_{i,r} > 0)` — the **selection frequency**, i.e. the share
  of repeats in which the elastic net assigned firm-year `i` a strictly
  positive proxy. Conceptually a stability score in the spirit of
  Meinshausen & Bühlmann's stability selection, but at the prediction level.
- `proxy_mean = mean(proxy_{i,r})` — the average held-out proxy across the
  200 repeats.

These two summaries are then the inputs to the prediction pipeline:
`p_i` drives the extensive margin (a tuned threshold gates which firm-years
receive any predicted emissions), and `proxy_mean` drives the within-cell
ranking among the classified emitters.

---

## 2. Is the proposal even alive? — distribution of `p_i`

Before doing anything else, we asked whether the 200 repeats are
sufficiently informative as a *distribution*, or whether they collapse to
near-binary {0, 1} for almost every firm-year (in which case `p_i` would
carry no information beyond the single-fit answer).

Script: [`diagnostic_selection_frequency.R`](analysis/active/diagnostic_selection_frequency.R)

Result on the full training panel of 26,608 firm-years:

| `p_i` bin | share |
|---|---|
| `= 0` | 38.5% |
| `(0, 0.05]` | 15.9% |
| `(0.05, 0.20]` | 10.5% |
| `(0.20, 0.40]` | 5.9% |
| `(0.40, 0.60]` | 3.1% |
| `(0.60, 0.80]` | 5.7% |
| `(0.80, 0.95]` | 7.3% |
| `(0.95, 1)` | 1.9% |
| `= 1` | 11.2% |

**32.3% of firm-years sit in the interior** `0.05 < p_i < 0.95`, well above
the 15% rule-of-thumb cutoff for the proposal to be worth pursuing. Splitting
by true emit status (`y > 0`):

| | n | mean p_i | share `p_i = 0` | share `p_i = 1` | share interior |
|---|---|---|---|---|---|
| emit = 0 | 23,513 | 0.196 | 43.5% | 1.6% | **35.6%** |
| emit = 1 | 3,095 | 0.964 | 0.8% | **84.6%** | 7.3% |

True emitters are essentially solved at the extensive margin: 85% have
`p_i = 1` and the mean is 0.96. The interior mass lives almost entirely in
the true non-emitters: 36% sit in an ambiguous middle band where the EN
sometimes selects them and sometimes does not. **This is the false-positive
region of any single-fit thresholded proxy**, and it is where a smarter
extensive-margin filter could plausibly add value.

Output: [`diag_selfreq_hist.png`](analysis/active/output/diag_selfreq_hist.png),
[`diag_selfreq_hist_by_emit.png`](analysis/active/output/diag_selfreq_hist_by_emit.png),
[`diag_selfreq_summary.rds`](analysis/active/output/diag_selfreq_summary.rds).

---

## 3. Tuning a threshold on `p_i`

### Design

Three sectors in the training panel contain both EU ETS firms (with verified
positive emissions) and non-ETS firms (assumed zero by construction):

- **paper** — NACE 2d 17 and 18, merged into "17/18"
- **refining** — NACE 2d 19
- **metals** — NACE 2d 24 and 25

These are exactly the three sectors used in panel B of
`table_main_with_mixed_crf.R`. We tune a threshold on `p_i` by
**leave-one-sector-out**:

1. For each held-out sector `g`:
   1. On each of the other two sectors `g'`, find the threshold maximising
      Youden's J = TPR − FPR via a pooled ROC over all firm-years in that
      sector.
   2. Average the two trained thresholds → `tau_g`.
   3. Apply `tau_g` to the held-out sector and compute pooled TPR / FPR / J.

Two schemes were compared:

- **(a) raw `p_i`** — global threshold on the [0, 1] scale.
- **(b) within-cell percentile of `p_i`** — convert `p_i` to a mid-rank
  percentile within (sector × year) cells; threshold the percentile. This
  controls for cross-sector heterogeneity in EN behaviour at the cost of
  discarding magnitude information.

Both were also evaluated for `log_revenue` as a size baseline.

Script: [`threshold_p_mixed_sectors.R`](analysis/active/threshold_p_mixed_sectors.R)

### Result

Held-out J averaged across the three LOSO folds:

| predictor | mean J | mean TPR | mean FPR |
|---|---|---|---|
| **p_i, scheme (a)** | **0.786** | 0.960 | 0.174 |
| log_revenue, scheme (a) | 0.728 | 0.940 | 0.212 |
| log_revenue, scheme (b) | 0.560 | 0.731 | 0.171 |
| p_i, scheme (b) | 0.491 | 0.642 | 0.151 |

**Scheme (a) wins decisively**, and raw `p_i` beats raw `log_revenue` by
0.06 J on average. The threshold `tau_g` on raw `p_i` lands in [0.97, 0.99]
for all three held-out sectors — a sharp rule of "predict emit iff `p_i` is
essentially 1".

**Why scheme (b) collapses**: the three mixed sectors have emit-rate base
rates differing by ~10× (refining is 34% emitters, paper is 0.6%, metals is
14%). Within-cell percentile thresholds do not transport across base rates.
Concretely, the threshold averaged from paper and metals (`τ ≈ 0.855` on
the percentile scale) is so high that *no* refining firm-year clears it,
giving held-out refining J = 0 in scheme (b). Magnitude information in
`p_i` is doing real work.

Output: [`threshold_p_mixed_sectors_results.rds`](analysis/active/output/threshold_p_mixed_sectors_results.rds).

---

## 4. Levels prediction: pipeline and metrics

### Pipeline

For each held-out sector `g`:

1. `tau_g` is the LOSO-tuned threshold from the previous step.
2. Classify firm-year `i` as a predicted emitter iff `p_i ≥ tau_g`.
3. Within each (`primary_crf_group`, year) cell, the calibration target is
   the **true** sum of `y` over the cell — exactly equivalent to panel B's
   `E_target = E_total − E_train`. This is "oracle" with respect to the cell
   total in both pipelines.
4. Among classified emitters in the cell, distribute `E_target` according to
   a chosen redistribution shape applied to `proxy_mean`.

Cells are `primary_crf_group × year` to match panel B exactly. The eval
subset is firm-years with `nace2d ∈ {17/18, 19, 24, 25}`, again matching
panel B. The threshold is sector-specific (one `tau_g` per held-out sector,
applied to firm-years based on which mixed sector they belong to), but the
redistribution loop runs on the full mixed subset in one pass.

Script: [`threshold_p_levels_eval_all.R`](analysis/active/threshold_p_levels_eval_all.R)
(supersedes the earlier
[`threshold_p_levels_eval.R`](analysis/active/threshold_p_levels_eval.R) and
[`threshold_p_levels_eval_gpa.R`](analysis/active/threshold_p_levels_eval_gpa.R)).

### Sample size

| sector | n | n_emit | n_nonemit |
|---|---|---|---|
| paper | 21,910 | 128 | 21,782 |
| refining | 179 | 61 | 118 |
| metals | 1,860 | 267 | 1,593 |
| **total** | **23,949** | **456** | **23,493** |

LOSO-tuned thresholds (raw p_i, scheme a):

- Held out paper: τ from refining 1.000 + metals 0.980 → **0.9900**
- Held out refining: τ from paper 0.965 + metals 0.980 → **0.9725**
- Held out metals: τ from paper 0.965 + refining 1.000 → **0.9825**

---

## 5. Choosing the redistribution shape: L-moment goodness of fit

### Concern

Comparing redistribution shapes only on mixed-sector eval results is
ad hoc — three sectors is not enough surface area to do model selection
without overfitting. We therefore picked a shape based on **goodness of fit
to a reference distribution that uses no mixed-sector data**.

### Reference distribution

`y > 0` emitter firm-years from CRF groups *other than* paper / refining /
metals (energy, chemicals, food, minerals, mfg_other, commercial,
agriculture, transport — 7 groups, 2,642 firm-years), with `log y` demeaned
by year and then by `primary_crf_group`.

Script: [`diag_lmoment_reference.R`](analysis/active/diag_lmoment_reference.R)

### L-moments

Sample L-moments of the reference distribution:

| moment | value |
|---|---|
| `l_1` | 0 (by construction, double-demeaned) |
| `l_2` | 0.835 |
| `t_3` (L-skewness) | **0.0079** ≈ 0 |
| `t_4` (L-kurtosis) | **0.198** |

The reference is essentially symmetric in the body (`t_3 ≈ 0`) but
substantially leptokurtic.

### Hosking goodness-of-fit Z

For each candidate 3-parameter family, fit by L-moments, compute the
implied `t_4`, and form `Z = (t_4_model − t_4_sample) / SE_boot(t_4_sample)`.
Bootstrap (B = 1000) gives `SE_boot(t_4) = 0.00713`. Distributions with
`|Z| ≤ 1.64` are acceptable at 90%.

| dist | t_4 model | t_4 sample | Δt_4 | Z | acceptable |
|---|---|---|---|---|---|
| **GLO** | 0.167 | 0.198 | −0.031 | **−4.35** | no (best of 5) |
| GNO | 0.123 | 0.198 | −0.075 | −10.5 | no |
| PE3 | 0.123 | 0.198 | −0.075 | −10.5 | no |
| GEV | 0.108 | 0.198 | −0.090 | −12.6 | no |
| GPA | 0.002 | 0.198 | −0.196 | **−27.5** | no (worst) |

**GLO is the best of the standard 3-parameter families on pure
goodness-of-fit grounds.** GPA — what panel B uses — is by far the worst:
its theoretical `t_4` at `t_3 ≈ 0` is essentially zero, so it cannot capture
the heavy tails of the demeaned log-emission distribution at all. All five
families are formally rejected at 90%, meaning no standard 3-parameter
family is a great fit for these data — GLO is simply the closest.

Selecting GLO this way is fully sector-agnostic: the reference distribution
contains zero mixed-sector observations, and the selection criterion
(`|Z|`) cannot be computed from mixed-sector data.

Output: [`diag_lmr_diagram.png`](analysis/active/output/diag_lmr_diagram.png),
[`diag_lmr_qq.png`](analysis/active/output/diag_lmr_qq.png),
[`diag_lmr_results.rds`](analysis/active/output/diag_lmr_results.rds).

---

## 6. Final pipeline performance

### Pooled across the three mixed sectors

| metric | panel B EN (prop) | this pipeline (GLO) |
|---|---|---|
| RMSE (kt) | 113.7 | 123.3 |
| MAPD | 0.687 | 0.858 |
| Pearson | **0.801** | 0.768 |
| Spearman | 0.329 | **0.666** |
| TPR | 0.995 | 0.956 |
| FPR | 0.199 | **0.022** |

Panel B numbers from the "Elastic Net" row of "Panel B: Mixed sectors"
in
[`inferring_emissions_output/table_main_with_mixed_crf_cv.tex`](../inferring_emissions_output/table_main_with_mixed_crf_cv.tex).

**Trade summary.** Compared to panel B, this pipeline:

- Cuts the false-positive rate by ~10× (0.199 → 0.022).
- Doubles within-sector-year Spearman (0.329 → 0.666).
- Loses ~0.04 TPR (0.995 → 0.956).
- Loses ~0.03 Pearson and ~10 kt RMSE on absolute levels.

The extensive-margin and ranking gains are large; the level-accuracy loss
is small and (as the diagnostics in section 7 show) attributable to
residual model misspecification of the redistribution shape, not to
anything fixable by tightening the classifier.

### Per held-out sector (classified mask, GLO redistribution)

| sector | RMSE (kt) | median APD | Pearson | Spearman |
|---|---|---|---|---|
| paper | 7.25 | 0.532 | 0.718 | 0.615 |
| refining | 776.9 | 0.608 | 0.770 | 0.838 |
| metals | 370.1 | 2.001 | 0.788 | 0.687 |

Output: [`threshold_p_levels_eval_all_results.rds`](analysis/active/output/threshold_p_levels_eval_all_results.rds).

---

## 7. Diagnostics that supported the design choices

### 7.1 The "median APD = 1.000" puzzle

In the initial run of
[`threshold_p_levels_eval.R`](analysis/active/threshold_p_levels_eval.R), the
pooled median APD came out at *exactly* 1.000 across both proportional and
GPA redistributions. Diagnosis:

Script: [`diag_apd.R`](analysis/active/diag_apd.R).

```
emitters with yhat == 0 (false negatives, APD = 1 exactly) :  20
emitters with APD < 1                                       : 212
emitters with APD > 1                                       : 224
total                                                       : 456
```

The 20 false negatives form a mass point at exactly APD = 1 spanning sorted
positions 213–232. The median index falls at positions 228–229, both inside
the mass — so `median(APD) = (1 + 1) / 2 = 1.000` exactly. A knife-edge
artifact, not a deep fact about the distribution. The substantive concern
underneath it is real, though: 49% of true emitters (244 of 456) have
APD ≥ 1, almost all of them in metals.

Per-sector breakdown:

| sector | n_emit | APD < 1 | APD = 1 (FN) | APD > 1 | median APD |
|---|---|---|---|---|---|
| paper | 128 | 94 | 11 | 23 | 0.65 |
| refining | 61 | 44 | 0 | 17 | 0.48 |
| **metals** | **267** | **74** | **9** | **184** | **3.25** |

The pooled median is dragged toward 1 entirely by metals.

### 7.2 Is the metals failure due to too many false positives, or to a wrong shape?

Two competing mechanisms could produce metals' levels-accuracy failure:

- **Mechanism (1)**: too many false positives among classified emitters
  diluting the cell mass. With 305 FP and 258 true emitters in metals
  classified-positives, this is plausible *a priori*.
- **Mechanism (2)**: the shape of `proxy_mean` doesn't match the shape of
  `y` in metals. Even with perfect classification, proportional
  redistribution by `proxy_mean` would misallocate magnitudes.

These are separable. Two diagnostics:

**Diagnostic A — does `proxy_mean` rank true emitters above false positives
within metals cells?** If yes, mechanism (1) can be fixed by a tighter
threshold and mechanism (2) is the real bottleneck.

Script: [`diag_metals_ranker.R`](analysis/active/diag_metals_ranker.R).

| | value |
|---|---|
| Mean within-year AUC of `proxy_mean` (true emit > FP) | **0.921** |
| Median within-year AUC | 0.923 |
| Pooled metals AUC | 0.913 |
| True-emitter median percentile from top of cell | 0.75 (top quartile) |
| False-positive median percentile from top of cell | 0.28 (bottom third) |

**The ranker is excellent.** `proxy_mean` correctly orders true emitters
above false positives 92% of the time. Mechanism (1) cannot be the main
story.

**Diagnostic B — oracle-classifier ablation.** Replace the LOSO-tuned `p_i`
threshold with the true emit label (`y > 0`) as the classification mask, run
the same redistribution pipeline, compare metrics.

Script: [`threshold_p_levels_eval_all.R`](analysis/active/threshold_p_levels_eval_all.R)
(the `oracle` rows in the output).

Pooled:

| mask | shape | RMSE | Pearson |
|---|---|---|---|
| classified | GLO | 123.3 | 0.768 |
| oracle | GLO | 121.5 | 0.776 |

Perfect classification buys ~2 kt of RMSE and 0.008 of Pearson over the
LOSO-tuned threshold. **Mechanism (1) is essentially not real.** The
levels-accuracy gap to panel B is overwhelmingly mechanism (2): the
`proxy_mean → predicted level` mapping is the wrong shape for metals'
heavy-tailed emission distribution. This is what motivated section 5 (the
shape selection is the only lever left) and what disposes of "tighten the
threshold" as a way forward.

### 7.3 Why GLO and not GPA / GEV / GNO / PE3?

In the levels-eval table, all four shapes were run side by side
([`threshold_p_levels_eval_all.R`](analysis/active/threshold_p_levels_eval_all.R)).
Pooled (classified mask):

| shape | RMSE | MAPD | Pearson | Spearman |
|---|---|---|---|---|
| prop | 147.6 | 1.000 | 0.658 | 0.665 |
| gpa | 145.4 | 1.000 | 0.652 | 0.667 |
| **glo** | **123.3** | **0.858** | **0.768** | 0.666 |
| gev | 132.9 | 0.918 | 0.722 | 0.666 |

GLO is best on every level metric, by a wide margin. But this on its own
would be ad hoc — GLO's pooled win is driven by metals (per-sector Pearson:
prop 0.52, GLO 0.79 in metals, with GLO actually *worse* than prop in paper
and refining), and three sectors is not enough to justify a shape choice on
post hoc levels-eval results alone.

The L-moment goodness-of-fit analysis in section 5 provides the
sector-agnostic justification: GLO is also the best 3-parameter family on
the non-mixed-sector reference distribution (Hosking |Z| = 4.35 vs the next
best at 10.5 and the worst at 27.5). The two diagnostics independently
agree, which is the basis for committing to GLO without overfitting to
metals.

---

## 8. Pipeline summary

The final pipeline:

1. **Inputs**: `proxy_matrix` from
   [`build_repeated_cv_proxy_crf.R`](analysis/active/build_repeated_cv_proxy_crf.R)
   (200 repeats × CRF-group K-fold CV).

2. **Per firm-year summaries**:
   - `p_i = mean(proxy_{i,r} > 0)` (selection frequency)
   - `proxy_mean = mean(proxy_{i,r})` (mean held-out proxy)

3. **Threshold tuning** (extensive margin):
   - Restrict to the three mixed sectors.
   - LOSO over {paper, refining, metals}: for each held-out sector, train
     Youden's-J-optimal `tau` on each of the other two sectors and average.
   - Apply the held-out `tau` to firm-years in the held-out sector.

4. **Redistribution shape selection**:
   - Build reference distribution from non-mixed-sector emitter firm-years
     (`y > 0`, primary CRF group ∉ {paper, refining, metals}), `log y`
     demeaned by year and primary CRF group.
   - Fit GPA / GLO / GEV / GNO / PE3 by L-moments.
   - Pick the family with the smallest Hosking |Z| → **GLO**.

5. **Levels prediction** (per held-out mixed sector):
   - For each (`primary_crf_group`, year) cell, calibration target = sum of
     true `y` in the cell.
   - Among firm-years with `p_i ≥ tau_g`, distribute the target via GLO
     weighted by within-cell ranks of `proxy_mean`.
   - Firm-years not classified as emitters get yhat = 0.

6. **Evaluation**: pool predictions across the three held-out sectors and
   compute metrics in [`utils/calc_metrics.R`](utils/calc_metrics.R).

Final pooled numbers vs panel B's "EN (prop)" mixed row:

| | RMSE (kt) | MAPD | Pearson | Spearman | TPR | FPR |
|---|---|---|---|---|---|---|
| panel B EN (prop) | **113.7** | **0.687** | **0.801** | 0.329 | 0.995 | 0.199 |
| this pipeline (GLO) | 123.3 | 0.858 | 0.768 | **0.666** | 0.956 | **0.022** |

---

## 9. Files

**Scripts** (all under [`analysis/active/`](analysis/active/)):

- [`diagnostic_selection_frequency.R`](analysis/active/diagnostic_selection_frequency.R) — section 2
- [`threshold_p_mixed_sectors.R`](analysis/active/threshold_p_mixed_sectors.R) — section 3
- [`threshold_p_levels_eval.R`](analysis/active/threshold_p_levels_eval.R) — first version of section 4 (proportional only, sector × year cells; superseded)
- [`threshold_p_levels_eval_gpa.R`](analysis/active/threshold_p_levels_eval_gpa.R) — GPA-only one-off (superseded)
- [`threshold_p_levels_eval_all.R`](analysis/active/threshold_p_levels_eval_all.R) — final levels eval (4 shapes × 2 masks, primary_crf_group × year cells), section 4 + 6 + 7.2 + 7.3
- [`diag_apd.R`](analysis/active/diag_apd.R) — section 7.1
- [`diag_metals_ranker.R`](analysis/active/diag_metals_ranker.R) — section 7.2 (diagnostic A)
- [`diag_lmoment_reference.R`](analysis/active/diag_lmoment_reference.R) — section 5

**Outputs** (all under [`analysis/active/output/`](analysis/active/output/)):

- `diag_selfreq_hist.png`, `diag_selfreq_hist_by_emit.png`, `diag_selfreq_summary.rds`
- `threshold_p_mixed_sectors_results.rds`
- `threshold_p_levels_eval_results.rds`, `threshold_p_levels_eval_gpa_results.rds`
- `threshold_p_levels_eval_all_results.rds`
- `diag_lmr_diagram.png`, `diag_lmr_qq.png`, `diag_lmr_results.rds`

**Reference points outside this work**:

- [`figures_tables/table_main_with_mixed_crf.R`](figures_tables/table_main_with_mixed_crf.R) — generates panel B
- [`inferring_emissions_output/table_main_with_mixed_crf_cv.tex`](../inferring_emissions_output/table_main_with_mixed_crf_cv.tex) — panel B output (located in `Documents/inferring_emissions_output/`)
- [`utils/calc_metrics.R`](utils/calc_metrics.R) — metric definitions
