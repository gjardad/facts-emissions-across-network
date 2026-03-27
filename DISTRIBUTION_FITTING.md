# Distribution Fitting: From Proxy to Firm-Level Emissions

This document explains the full pipeline for converting EN proxy values (in asinh space) to firm-level emission predictions in levels, and the reasoning behind each design choice.

## 1. Why the EN is trained in asinh space

The elastic net predicts asinh(y), not y directly. Two reasons:

**Scale compression.** Emissions range from 0 to ~18 million tonnes. Minimizing squared error in levels would be dominated by fitting the few largest emitters — a coefficient of 0.001 matters enormously for a firm emitting 1M tonnes but is negligible for a firm emitting 100 tonnes. The asinh transformation compresses the scale (asinh(100) ≈ 5.3, asinh(1M) ≈ 14.5), so the EN treats prediction errors at small and large emitters as roughly comparable. This lets it learn from the full cross-section.

**Handling zeros.** Log(y) would achieve the same compression but is undefined at y = 0. Our training sample includes confirmed non-emitters (y = 0 in sectors 17/18, 19, 24). asinh(0) = 0, so the transformation handles zeros cleanly.

## 2. The retransformation problem

The EN produces a prediction in asinh space — the "proxy." To get firm-level emissions in tonnes, we need to convert back to levels. The naive approach is to invert the transformation: apply sinh to the proxy.

### 2.1 Why sinh is convex and what Jensen's inequality implies

A function f is convex if f''(x) > 0. For sinh: f'(x) = cosh(x), f''(x) = sinh(x) > 0 for all x > 0. Geometrically, the curve bends upward.

Jensen's inequality states that for any convex f and random variable X:

> f(E[X]) ≤ E[f(X)]

In our context: let X = proxy + ε, where proxy = E[asinh(y) | features] and ε is the prediction error. We compute sinh(proxy) = sinh(E[X]). What we want is E[y] = E[sinh(X)]. By Jensen:

> sinh(E[X]) ≤ E[sinh(X)] = E[y]

So sinh(proxy) — what we compute — is systematically below E[sinh(X)] — what we want. We underestimate.

### 2.2 Why underestimation is worse for large emitters

The second-order Taylor expansion of the Jensen gap is:

> E[sinh(X)] − sinh(E[X]) ≈ ½ sinh''(E[X]) × Var(ε) = ½ sinh(E[X]) × Var(ε)

Since sinh grows exponentially, this gap is exponentially larger at high proxy values. At proxy = 5 (small emitter), sinh(5) ≈ 74. At proxy = 15 (large emitter), sinh(15) ≈ 1,634,508. Even if Var(ε) is the same, the gap is 22,000× larger for the large emitter. The underestimation is driven overwhelmingly by where you sit on the curve, not by how noisy the prediction is.

### 2.3 Why this leads to extreme concentration of shares

When we use sinh values as proportional weights within a (sector, year) cell:

> share_i = sinh(proxy_i) / Σⱼ sinh(proxy_j)

the ratio sinh(a)/sinh(b) for a > b grows exponentially in (a − b). For large arguments, sinh(x) ≈ eˣ/2, so sinh(a)/sinh(b) ≈ e^(a−b). If two firms differ by 5 units in asinh space, their share ratio is ≈ e⁵ ≈ 148:1, even if the actual emission ratio is 5:1.

This is not specific to asinh — any transformation where the inverse maps additive differences to multiplicative ratios will do this. The issue is that the EN's prediction errors in asinh space, which may be moderate (±2 units), become extreme ratios in levels (e² ≈ 7× to e⁻² ≈ 0.13×).

Empirically, the top firm in a typical cell absorbs 50–80% of the cell's total emissions under sinh-calibrated allocation, vastly exceeding its true share. This inflates within-sector dispersion: p90/p10 bias of +412.

### 2.4 Why Duan's smearing estimator doesn't help

**What Duan's smearing does.** Instead of sinh(prediction), compute:

> Ê[y|X] = (1/n) Σᵢ sinh(prediction + eᵢ)

averaging sinh over the training residuals eᵢ. This "smears" the prediction over the residual distribution, correctly accounting for the nonlinearity. Under homoscedasticity, this gives a constant multiplicative correction S that cancels in proportional allocation (shares are unchanged).

**Why heteroscedasticity matters.** If residual variance varies with emission size, the smearing correction is observation-specific. You need to use residuals from similar firms.

**Why our pattern is the wrong direction.** Our residual variance *decreases* with emission size (from the Jensen residual diagnostic). This means:
- Small emitters: wide residuals → large smearing correction → share increases
- Large emitters: tight residuals → small smearing correction → share decreases

But the problem we're trying to fix is that large emitters' shares are already too high (from sinh concentration). Duan would push in the wrong direction — increasing small emitters' shares further while barely touching the large emitters.

The retransformation bias itself is driven primarily by the convexity of sinh at high values, not by residual variance. Even well-predicted large emitters are massively underestimated in levels because sinh amplifies even small asinh-space errors enormously (imputed/actual ≈ 0.21 for large emitters).

## 3. The fix: separate ranking from magnitude

### 3.1 Core insight

The EN proxy is good at *ranking* firms (who emits more than whom within a sector) but bad at assigning *magnitudes* (how much more), because sinh distorts magnitudes. We decouple the two:

- **Ranking:** proxy determines the ordering within each cell
- **Magnitude:** a reference distribution estimated from training data determines the within-cell shape
- **Level:** calibration to the known sector-year total E_target determines the overall scale

### 3.2 Why we trust ranks despite Pearson > Spearman

The Pearson correlation in levels (0.85) exceeds the Spearman rank correlation (0.64). This seems to contradict the claim that the proxy is "better at ranking than at levels." The resolution:

Pearson is variance-weighted. Observations far from the mean contribute more to both numerator and denominator. In our data, the top 5 emitters in a sector account for ~80% of the cross-sectional variance. If the model correctly identifies these 5 firms and assigns them roughly the right magnitudes, Pearson is high — even if the ordering of the remaining 95% of firms is garbled.

Spearman replaces values with ranks, giving every firm equal weight. It measures ordering across the full distribution. It's lower because the model is mediocre at ordering the many small emitters.

For distributional claims (carbon productivity dispersion, emission concentration), we need the full distribution to be right, not just the top. The ranking information is usable across the full distribution; the level information is only reliable for the top few firms that drive the Pearson. We keep the ranking and replace the levels.

## 4. The reference distribution

### 4.1 Construction

For each CV fold k, we estimate the within-sector shape from training emitters (fold_k ≠ k, y > 0):

1. **Year-demean:** tilde_j = log(y_j) − μ_t, where μ_t = mean of log(y) across all emitters in year t. Removes time trends.
2. **Sector-demean:** d_j = tilde_j − μ_s, where μ_s = mean of tilde within sector s (pooling across years). Removes sector levels.
3. **Pool** d_j across all training sectors. Gives the distribution of within-sector deviations in log space, net of year effects.

We pool across sectors because held-out sectors have no training data (mimics deployment). Year-demeaning before sector-demeaning avoids fragmenting into thin sector-year cells.

### 4.2 L-moment diagnostic: sectors don't share a common shape

We checked whether sectors share a common within-sector distributional shape using L-moment ratio diagrams — a standard diagnostic from regional frequency analysis in hydrology.

**L-moments** are linear combinations of order statistics: λ₁ (mean), λ₂ (L-scale), with ratios τ₃ = λ₃/λ₂ (L-skewness) and τ₄ = λ₄/λ₂ (L-kurtosis) characterizing the shape. Unlike conventional moments, L-moments are nearly unbiased in small samples and robust to outliers.

Different parametric families (Normal, Log-normal, GEV, GPA, GLO) trace distinct curves on the (τ₃, τ₄) plane. If sectors clustered along one curve, a common parametric family would be justified.

**Result:** Sectors scatter widely. L-skewness ranges from −0.50 to +0.44 (IQR = 0.27), L-kurtosis from −0.19 to +0.34 (IQR = 0.23). No single family fits all sectors. The pooled distribution (τ₃ ≈ 0, τ₄ ≈ 0.21) is a compromise that doesn't represent any individual sector well.

**Implication:** The "common shape" assumption is strained. We proceed anyway — the pooled shape is better than the sinh-distorted shape, even if imperfect for individual sectors. This is a limitation we acknowledge. See `figure_lmoment_ratio_diagram.pdf` and `table_lmoment_sectors.tex`.

### 4.3 Fitting the GPA via L-moments

Rather than using the raw empirical quantile function (quantile mapping), we fit a **Generalized Pareto distribution (GPA)** to the pooled deviations. The GPA has CDF:

> F(x) = 1 − [1 − k(x − ξ)/α]^(1/k)

with three parameters: location ξ, scale α, shape k.

**Fitting via L-moments** means computing sample L-moments (λ̂₁, λ̂₂, λ̂₃) from the pooled deviations, then solving for (ξ, α, k) such that the GPA's theoretical L-moments match the sample. This is method-of-moments estimation using L-moments instead of conventional moments. The `lmom` package provides `pelgpa()` for this inversion and `quagpa()` for the quantile function.

**Why GPA over empirical quantiles:** The GPA gives a smooth parametric quantile function, less sensitive to the specific composition of training sectors in each fold. Empirically, GPA slightly outperforms empirical quantile mapping on distributional metrics (p90/p10 RMSE: 86 vs 90).

## 5. Redistribution within each cell

For a held-out (sector, year) cell with n firms predicted as emitters:

1. **Rank** firms by proxy: highest proxy = rank n, lowest = rank 1.

2. **Map rank to quantile position:** p_i = (rank_i − 0.5) / n. This is the Hazen plotting position — it maps rank 1 to 0.5/n (just above 0) and rank n to (n−0.5)/n (just below 1), placing each rank at the midpoint of its probability bin.

3. **Look up GPA quantile:** w_i = Q_GPA(p_i). The GPA quantile function maps probability p to a deviation in log space — how far a firm at that quantile is from the sector mean. The top-ranked firm gets a large positive deviation (well above sector mean), the bottom-ranked firm gets a negative deviation (well below).

4. **Convert to shares and calibrate:**

> yhat_i = E_target × exp(w_i) / Σⱼ exp(w_j)

The w_i are log-space deviations. exp(w_i) converts to multiplicative factors relative to the sector mean. Normalizing by the sum gives shares; multiplying by E_target gives tonnes. This is equivalent to log(yhat_i) = c + w_i, where c absorbs the sector-year level.

## 6. Cross-validated threshold for the extensive margin

### 6.1 The problem

The extensive margin (proxy > 0 = emitter, proxy = 0 = non-emitter) produces false positives: non-emitters with positive proxy values. Under sinh-calibrated allocation, false positives got scraps (the top firm absorbed most of E_target), so their severity was artificially low. Under Pareto redistribution, false positives get meaningful emission levels, exposing the underlying misclassification.

### 6.2 The fix

We apply a percentile threshold to the proxy within each (sector, year) cell. Firms below the threshold are classified as non-emitters (zero emissions) regardless of their proxy value.

**Tuning:** For each of the three mixed sectors (17/18, 19, 24) — the only sectors where we observe both emitters and non-emitters — find the percentile p* of proxy > 0 values that maximizes Youden's J = TPR − FPR.

**Cross-validation (LOSO):** For held-out sector h, the threshold is the average of p* from the other two mixed sectors. This avoids using the held-out sector's data for tuning.

**Application:** Within each (sector, year) cell of a mixed sector, zero out firms below the p*-th percentile of proxy > 0 values. For pure-emitter sectors (EU ETS only), no threshold is applied — the proxy > 0 rule is unchanged.

**Why percentiles, not absolute values:** The proxy distribution differs across sectors, so an absolute threshold of τ = 0.5 means very different things in sector 17/18 vs 19. A percentile threshold is scale-invariant and more likely to transport across sectors.

### 6.3 Results

CV percentile thresholds (EN proxy):

| Held-out sector | Trained on | Avg p* | OOS TPR | OOS FPR |
|---|---|:-:|:-:|:-:|
| Paper & printing (17/18) | 19, 24 | 0.544 | 0.981 | 0.067 |
| Petroleum refining (19) | 17/18, 24 | 0.690 | 0.825 | 0.068 |
| Iron & steel (24) | 17/18, 19 | 0.614 | 0.992 | 0.242 |

The threshold works well for 17/18 and 19 (FPR drops to ~0.07 with high TPR). Sector 24 retains higher FPR (0.242) — the threshold trained on 17/18 and 19 isn't aggressive enough for 24's composition. This transportability limitation is inherent: the emitter-to-non-emitter ratio varies across sectors.

### 6.4 Transportability concern

Thresholds tuned on the three mixed sectors must be applied at deployment to sectors with unknown emitter compositions. The LOSO results show the threshold is reasonably stable (p* ≈ 0.54–0.69) but not perfectly transportable — sector 19 (34% emitter share) needs a different threshold than sector 17/18 (0.6% emitter share). This is a limitation we acknowledge.

## 7. Evaluation

### 7.1 Main results (all models, Pareto + CV threshold)

20 repeats, K=5 sector-level CV. All models use Pareto redistribution with CV threshold.

| Ranking signal | RMSE (kt) | nRMSE | MAPD | Pearson | Spearman | FPR | TPR | p50 | p99 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Revenue | 218.3 | 1.000 | 0.839 | 0.768 | 0.560 | 0.366 | 0.985 | 0.000 | 0.329 |
| **Elastic Net** | **185.2** | **0.848** | **0.777** | **0.838** | **0.767** | **0.079** | **0.960** | **0.000** | **0.265** |
| NACE | 242.7 | 1.112 | 0.825 | 0.709 | 0.577 | 0.301 | 0.981 | 0.000 | 0.549 |
| Gated Rev | 204.3 | 0.936 | 0.828 | 0.799 | 0.742 | 0.092 | 0.955 | 0.000 | 0.344 |

The EN proxy is the best ranking signal: lowest nRMSE (0.848), highest Pearson (0.838), lowest FPR (0.079) among non-revenue models.

### 7.2 Redistribution method comparison (EN proxy only)

| Method | p90/p10 bias | p90/p10 RMSE | Gini RMSE | Pearson | Spearman | RMSE (levels) |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| Sinh-calibrated | +412 | 850 | 0.24 | 0.849 | 0.512 | 625,229 |
| Quantile mapping | −41 | 90 | 0.19 | 0.849 | 0.694 | 510,110 |
| Pareto (GPA) | −26 | 86 | 0.20 | 0.850 | 0.700 | 514,079 |

Both alternatives reduce distributional bias by ~10× relative to sinh. Pareto has a slight edge and is our preferred method. Quantile mapping results are available for robustness.

### 7.3 Spearman improvement is real, not a bug

Spearman improves from 0.51 (sinh) to 0.70 (Pareto) to 0.77 (Pareto + threshold). This reflects correction of cross-cell magnitude distortion. Under sinh, the top firm in each cell gets 50–80% of emissions, inflating predictions for top-ranked firms in small cells beyond genuinely large emitters in other cells. Pareto assigns realistic within-cell shares, fixing cross-cell comparisons. A pure within-cell-rank Spearman (stripping all magnitudes) is 0.53, confirming that sinh magnitudes actively hurt the global ranking.

## 8. Deployment considerations

### 8.1 How deployment works

At deployment, for each (sector, year) cell:

1. **ETS firms:** Emissions are observed directly — no prediction needed.
2. **Non-ETS firms:** We compute E_non-ETS = E_total (from NIR) − E_ETS (observed). This budget is distributed among non-ETS firms using the same pipeline: rank by EN proxy, apply percentile threshold (for mixed sectors), assign Pareto-shaped emissions calibrated to E_non-ETS.

### 8.2 Shape mismatch between training and deployment populations

The GPA reference distribution was estimated from within-sector deviations of **ETS emitters** — firms above the regulation threshold. At deployment, we apply this shape to **non-ETS firms** — firms below the threshold. The calibration constant ensures the *level* is correct (predictions sum to E_non-ETS). The tension is about the *shape*.

The within-group dispersion of ETS emitters reflects variation among large emitters. The within-group dispersion of non-ETS emitters could be different:

- **If non-ETS emitters are more homogeneous** (all "small" by definition — below the ETS threshold), the ETS-estimated spread is too wide. We'd overstate the dispersion among non-ETS firms: the top-ranked non-ETS firm gets too much, the bottom-ranked gets too little.
- **If non-ETS emitters are more heterogeneous** (spanning from a bakery with negligible combustion to a medium-sized plant just below the threshold), the ETS-estimated spread is too narrow.

We cannot sign this bias without observing non-ETS emissions — which is the whole point of the exercise. What we can say is that this is a **second-order concern** relative to the first-order improvement: the GPA shape is estimated from data (even if a different population) rather than imposed by the sinh transformation (which has no empirical basis at all). The level is correct by construction via calibration. The shape is where the residual uncertainty lives.

### 8.3 Threshold transportability at deployment

The CV percentile threshold (~0.55–0.69) was tuned on three mixed sectors (17/18, 19, 24). At deployment, it must be applied to sectors with unknown emitter compositions. The LOSO diagnostic shows the threshold is reasonably stable across the three mixed sectors but not perfectly transportable — sector 19 (34% emitter share) needs a different threshold than sector 17/18 (0.6% emitter share).

At deployment, the emitter share of non-ETS firms is unknown by definition. A conservative choice would be to use the pooled threshold (~0.70 percentile, tuned on all three mixed sectors together), accepting that it may be too aggressive for some sectors and too lenient for others.

## 9. Scripts

| Script | Purpose | Location |
|---|---|---|
| `diagnostic_lmoment_sectors.R` | L-moment ratio diagnostic | `analysis/active/` |
| `diagnostic_redistribution_comparison.R` | Sinh vs quantile mapping vs Pareto comparison | `analysis/active/` |
| `diagnostic_pareto_threshold.R` | LOSO threshold tuning on mixed sectors | `analysis/active/` |
| `table_pareto_all_models.R` | Main results + zero-emitter tables, Pareto, no threshold | `figures_tables/` |
| `table_pareto_all_models_cvthresh.R` | Main results + zero-emitter tables, Pareto + CV threshold | `figures_tables/` |
| `table_pareto_redistribution.R` | Mixed table: sinh vs Pareto for EN, plus all models | `figures_tables/` |
| `figure_lmoment_diagnostic.R` | Publication-quality L-moment ratio diagram + table | `figures_tables/` |

## 9. Data

| File | Contents |
|---|---|
| `repeated_cv_proxy_sector_asinh.RData` | proxy_matrix (26608 × 200), panel, syt |
| `firm_year_panel_with_proxies.RData` | Revenue, NACE proxy (proxy_tabachova) |
| `training_sample.RData` | Full training panel with emissions |
