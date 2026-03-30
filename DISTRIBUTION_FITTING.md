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

### 5.1 Overview

The core idea: given a known total E_target for a group of firms, distribute it among ranked firms so that the relative spacing between firms in log-emission space matches the shape of a fitted GPA distribution.

The procedure has three ingredients:
- A **ranking** of firms (from the EN proxy)
- A **distributional shape** (from the fitted GPA)
- A **total** to calibrate to (E_target, known from NIR minus ETS)

The ranking tells us *who* emits more; the GPA tells us *how much more*; the total tells us *how much in aggregate*.

### 5.2 Formal definitions

**Setup.** Consider a group g (a sector-year cell in CV, or a CRF-group × year cell at deployment) with:
- n firms classified as emitters (proxy > 0, above the CV threshold if applicable)
- A known total E_target to distribute among them
- Proxy values {v₁, ..., vₙ} from the EN, where vᵢ = sinh(proxyᵢ) > 0

**Step 1: Rank.** Assign integer ranks r₁, ..., rₙ by ascending proxy value:

> rᵢ = rank(vᵢ),    rᵢ ∈ {1, ..., n}

The firm with the smallest proxy gets rank 1 (bottom of the distribution), the firm with the largest proxy gets rank n (top).

**Step 2: Plotting positions.** Convert ranks to probabilities using the Hazen plotting position:

> pᵢ = (rᵢ − 0.5) / N

where N is the total number of emitters in the distribution. In the CV evaluation setting, N = n (only the firms in the held-out cell). At deployment, N = n_deploy + n_ETS, where deployment firms occupy ranks 1, ..., n_deploy and ETS firms implicitly occupy ranks n_deploy + 1, ..., N (see section 5.5).

The Hazen position maps rank 1 to (0.5)/N ≈ 0 and rank n to (n − 0.5)/N ≈ n/N, placing each rank at the midpoint of its probability bin. This avoids boundary artifacts (p = 0 or p = 1 would produce infinite quantiles under many distributions).

**Step 3: GPA quantiles.** Evaluate the GPA quantile function at each plotting position:

> wᵢ = Q_GPA(pᵢ; ξ, α, k)

where Q_GPA is the inverse CDF of the Generalized Pareto distribution:

> Q_GPA(p; ξ, α, k) = ξ + (α/k) [1 − (1 − p)^k]

Each wᵢ is a deviation in log-emission space: how far above or below the sector mean a firm at quantile pᵢ should be. Firms near the top (high pᵢ) get large positive wᵢ; firms near the bottom get negative wᵢ.

**Step 4: Exponentiate to get unnormalized weights.**

> ωᵢ = exp(wᵢ − max(w))

Subtracting max(w) before exponentiating is a numerical stability trick — it doesn't affect the shares. The weights ωᵢ are proportional to the emission level each firm should have relative to the others: if wᵢ − wⱼ = 1, firm i should emit exp(1) ≈ 2.72 times more than firm j.

**Step 5: Normalize and calibrate.**

> ŷᵢ = E_target × ωᵢ / Σⱼ ωⱼ

This ensures Σᵢ ŷᵢ = E_target exactly. The share of firm i is sᵢ = ωᵢ / Σⱼ ωⱼ.

**What the GPA controls.** The GPA parameters (ξ, α, k) determine only the *relative spacing* between quantiles — how spread out firms are in log-emission space. The absolute level is entirely determined by E_target and the number of firms n. Two groups with the same GPA but different E_target or n will have different emission levels but the same Gini coefficient, the same p90/p10 ratio, and the same variance of log emissions. This is the key property: the GPA governs *dispersion*, while calibration governs *level*.

### 5.3 Worked example (CV evaluation setting)

Suppose we have a held-out sector-year cell with 6 firms predicted as emitters and E_target = 1,000 tonnes. The fitted GPA has parameters ξ = −1.05, α = 1.62, k = −0.10 (representative values from our data).

**Step 1.** The firms have proxy values (in levels, after sinh): {12, 45, 78, 150, 310, 820}. Ranks: {1, 2, 3, 4, 5, 6}.

**Step 2.** Plotting positions (N = 6, same as n in CV):

| Firm | Proxy | Rank rᵢ | pᵢ = (rᵢ − 0.5)/6 |
|---|---|---|---|
| A | 12 | 1 | 0.083 |
| B | 45 | 2 | 0.250 |
| C | 78 | 3 | 0.417 |
| D | 150 | 4 | 0.583 |
| E | 310 | 5 | 0.750 |
| F | 820 | 6 | 0.917 |

**Step 3.** GPA quantiles wᵢ = Q_GPA(pᵢ; −1.05, 1.62, −0.10):

Using Q_GPA(p) = ξ + (α/k)[1 − (1−p)^k] = −1.05 + (1.62/−0.10)[1 − (1−p)^(−0.10)]:

| Firm | pᵢ | wᵢ = Q_GPA(pᵢ) |
|---|---|---|
| A | 0.083 | −0.908 |
| B | 0.250 | −0.577 |
| C | 0.417 | −0.153 |
| D | 0.583 | 0.432 |
| E | 0.750 | 1.359 |
| F | 0.917 | 3.520 |

Note how the spacing increases toward the top — the gap between firms E and F (2.16 units in log space) is much larger than between A and B (0.33 units). This is the heavy right tail of the GPA: the top firm is disproportionately larger than the rest.

**Step 4.** Unnormalized weights (subtract max = 3.520):

| Firm | wᵢ − 3.520 | ωᵢ = exp(wᵢ − 3.520) |
|---|---|---|
| A | −4.428 | 0.0119 |
| B | −4.097 | 0.0166 |
| C | −3.673 | 0.0254 |
| D | −3.088 | 0.0456 |
| E | −2.161 | 0.1152 |
| F | 0.000 | 1.0000 |
| **Total** | | **1.2148** |

**Step 5.** Shares and calibrated emissions (E_target = 1,000):

| Firm | Share sᵢ | ŷᵢ (tonnes) |
|---|---|---|
| A | 1.0% | 9.8 |
| B | 1.4% | 13.7 |
| C | 2.1% | 20.9 |
| D | 3.8% | 37.5 |
| E | 9.5% | 94.8 |
| F | 82.3% | 823.2 |
| **Total** | **100%** | **1,000.0** |

**Interpretation.** The top-ranked firm (F) gets 82% of the total — this reflects the heavy tail of the within-sector emission distribution. The ratio between the top and bottom firms is 823.2 / 9.8 ≈ 84:1 in emission levels, corresponding to wF − wA = 3.520 − (−0.908) = 4.43 in log space (exp(4.43) ≈ 84). This ratio is entirely determined by the GPA shape, not by the proxy values themselves — the proxy only determines who is ranked where.

**The role of the proxy.** Notice that firm F has a proxy 68× larger than firm A (820 vs 12), but gets only 84× more emissions. Under sinh-calibrated allocation, the ratio would be sinh(proxy_F)/sinh(proxy_A) — potentially thousands-to-one depending on the asinh-space values. The Pareto approach discards the magnitude information in the proxy and uses only the ranking, replacing magnitude with the empirically-estimated GPA shape.

### 5.4 What "matching the Pareto" means precisely

When we say the within-sector distribution "matches" the fitted GPA, we mean:

1. **The quantile function is preserved.** If you take the imputed emissions {ŷ₁, ..., ŷₙ}, compute log(ŷᵢ), demean them within the cell, sort them, and plot the i-th smallest against the plotting position pᵢ = (i − 0.5)/n, the resulting curve traces the GPA quantile function Q_GPA(p). This is by construction.

2. **Dispersion statistics are determined by the GPA.** The Gini coefficient, p90/p10 ratio, and variance of log emissions of {ŷ₁, ..., ŷₙ} are all functions of the GPA parameters (ξ, α, k) and the number of firms n, not of the proxy values or E_target. Changing E_target rescales all ŷᵢ proportionally, leaving dispersion statistics unchanged. Changing the ranking (swapping which firm gets rank 3 vs rank 4) changes who gets what, but not the distribution of emissions across quantiles.

3. **What is NOT matched.** The GPA determines the shape of the distribution within each cell but not:
   - The level (set by E_target)
   - Which firm gets which quantile position (set by the proxy ranking)
   - The number of emitters (set by the CV threshold)

### 5.5 Deployment: combined ranking with ETS firms

At deployment, each CRF-group × year cell may contain both ETS firms (with observed emissions) and deployment firms (to be imputed). The key difference from the CV setting is the **combined ranking**: deployment firms are placed *below* ETS firms in the within-cell distribution.

**Setup.** A CRF-group × year cell has:
- n_ETS firms with observed emissions {e₁, ..., e_{n_ETS}}
- n_deploy deployment firms classified as emitters (proxy > threshold)
- E_deploy = E_NIR − E_ETS to distribute among deployment firms
- N = n_ETS + n_deploy total emitters in the cell

**Ranking convention.** Deployment firms occupy ranks 1, ..., n_deploy (bottom of the distribution). ETS firms occupy ranks n_deploy + 1, ..., N (top of the distribution). Within each group, firms are ranked by their proxy (deployment) or observed emissions (ETS). The implicit assumption is that every deployment firm emits less than every ETS firm.

**Plotting positions for deployment firms.** With N total emitters:

> pᵢ = (rᵢ − 0.5) / N,    rᵢ ∈ {1, ..., n_deploy}

Since rᵢ ≤ n_deploy < N, all deployment plotting positions satisfy pᵢ < n_deploy/N < 1. The deployment firms occupy only the lower portion of the distribution. The highest-ranked deployment firm has p = (n_deploy − 0.5)/N, which is well below 1 when n_ETS is large relative to n_deploy. This compresses all deployment firms into the left tail of the GPA.

**Worked example.** CRF group "food", year 2010:
- 3 ETS firms: emissions = {50, 100, 150} tonnes. min_ets_emit = 50.
- 5 deployment firms with proxy > threshold. E_deploy = 200 tonnes.
- N = 3 + 5 = 8 total emitters.

Deployment firm plotting positions (ranks 1–5 out of N = 8):

| Deploy firm | Rank | pᵢ = (rank − 0.5)/8 |
|---|---|---|
| d1 (lowest proxy) | 1 | 0.0625 |
| d2 | 2 | 0.1875 |
| d3 | 3 | 0.3125 |
| d4 | 4 | 0.4375 |
| d5 (highest proxy) | 5 | 0.5625 |

Compare with the CV setting where N = n = 5: the same 5 firms would have positions {0.10, 0.30, 0.50, 0.70, 0.90}. With 3 ETS firms added to the combined ranking, the deployment firms are compressed into the range [0.06, 0.56] instead of [0.10, 0.90]. This has two effects:

1. **Lower quantiles.** Each deployment firm maps to a smaller GPA quantile, yielding a more negative wᵢ and a smaller share. The top deployment firm gets Q_GPA(0.56) instead of Q_GPA(0.90) — much less extreme.

2. **Less dispersion among deployment firms.** The spread of plotting positions is narrower (0.50 range vs 0.80 range), so the ratio between the top and bottom deployment firm is smaller.

Using the same GPA as the earlier example (ξ = −1.05, α = 1.62, k = −0.10):

| Deploy firm | pᵢ | wᵢ | ωᵢ (unnorm.) | Share | ŷᵢ (tonnes) |
|---|---|---|---|---|---|
| d1 | 0.0625 | −0.945 | 0.275 | 10.1% | 20.2 |
| d2 | 0.1875 | −0.710 | 0.348 | 12.7% | 25.5 |
| d3 | 0.3125 | −0.431 | 0.460 | 16.8% | 33.7 |
| d4 | 0.4375 | −0.091 | 0.646 | 23.7% | 47.4 |
| d5 | 0.5625 | 0.346 | 1.000 | 36.7% | 73.3 |
| **Total** | | | **2.728** | **100%** | **200.0** |

Note: the shares are much more even (10%–37%) compared to the CV example (1%–82%). This is because all deployment firms are compressed into the lower-middle portion of the GPA, where the quantile function is relatively flat.

**Upper-bound constraint check.** max(ŷᵢ) = 73.3. Is 73.3 ≥ min_ets_emit = 50? **Yes — violated.** The top deployment firm would emit more than the smallest ETS firm, breaking the assumption that deployment firms are smaller.

In this case, the code would:
1. Halve the CV threshold and retry (adding more firms → larger N → smaller plotting positions → smaller shares per firm)
2. If the threshold reaches zero and the constraint is still violated, **cap** the violating firms at min_ets_emit × (1 − 10⁻⁶) and redistribute the freed emissions proportionally to uncapped firms

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

## 9. Physics-based upper bound on non-ETS firm emissions

### 9.1 The constraint problem

At deployment, mixed CRF-group × year cells contain both ETS firms (observed) and deployment firms (imputed). The Pareto allocation places deployment firms below ETS firms in the combined ranking, which implies that the highest-emission deployment firm should emit less than the lowest-emission ETS firm. In practice, this constraint is violated in 86% of mixed sector-years (see section 5.5), because the Pareto shape concentrates too much weight on top-ranked deployment firms relative to the small ETS firms at the bottom of the ETS distribution.

The original constraint — deployment < min(ETS within CRF group) — is problematic for two reasons. First, the minimum ETS emitter within a CRF group may be an outlier (a plant that shut down mid-year) or a firm in a different NACE sector with structurally lower emissions. Second, in the CRF groups that aggregate many NACE sectors (e.g., mfg_other), the minimum ETS emitter can be as low as 24 tonnes, constraining deployment firms in unrelated sectors with much higher typical emissions.

### 9.2 Derivation of the 30 kt CO2/year cap

We replace the data-driven min(ETS) constraint with a physics-based upper bound derived from the EU ETS regulatory threshold.

**Regulatory threshold.** The EU ETS covers combustion installations with a rated thermal input exceeding 20 MW (Directive 2003/87/EC, Annex I). Non-ETS firms, by definition, operate combustion installations at or below this capacity.

**Fuel assumption: natural gas.** We assume the marginal non-ETS installation combusts natural gas. This is justified for Belgium by three lines of evidence:

1. **Aggregate fuel mix.** The Belgian NID (2025 submission) reports that gaseous fuels accounted for 55–65% of manufacturing energy consumption during 2005–2011, rising from 45% in 1990 to 70% by 2021 (NID Section 3.1.1.1, Figure 3.2.b). In the commercial and institutional sector, natural gas and gaseous fuels represent 77% of energy consumption (NID Section 3.1.1.2). These shares are for all installations, including large ETS-covered plants.

2. **Concentration of solid fuels in large ETS installations.** Coal and coke consumption in Belgium is concentrated in a small number of large facilities that are all covered by the EU ETS: blast furnaces and coke ovens in the iron and steel sector (NACE 24), cement kilns (NACE 23), and large power plants (NACE 35). The NID notes that coal was phased out of power generation in the Flemish region by 2017, and was already concentrated in a handful of large Electrabel plants during 2005–2011 (NID Section 3.2.6). None of the sub-20MW emission sources discussed in the NID use solid fuels.

3. **Infrastructure and economics.** Belgium has extensive natural gas pipeline infrastructure operated by Fluxys (high-pressure transmission) and regional distribution system operators. The NID documents that DSOs have reported gas offtakes per NACE code since 2005, confirming pipeline access at the firm level across industrial and commercial sectors (NID Section 3.2.5). Coal combustion at sub-20MW scale requires dedicated bulk storage, handling equipment, and specialized boilers whose capital costs are disproportionate to the thermal capacity. Combined with increasingly restrictive air quality regulations on particulates and SO₂ in Flanders, Brussels, and Wallonia, sub-20MW coal installations are not economically viable in Belgium's context.

Using the natural gas assumption is conservative: the IPCC 2006 CO₂ emission factor for natural gas (56.1 t CO₂/TJ) is the lowest among fossil fuels. A non-ETS installation burning gas oil, heavy fuel oil, or coal at the same capacity would emit 1.3–1.7× more, making 30 kt a conservative upper bound.

**Capacity factor assumption: 60%.** We assume the installation operates at 60% of its rated capacity on an annual basis (≈ 5,256 full-load hours per year). This is the midpoint of the 50–70% range typical for industrial boilers and small combined heat-and-power units, based on:

- The IPCC Good Practice Guidance (2000, Section 2.2) notes that industrial stationary combustion sources operate intermittently depending on process demand, with utilization rates varying by application: process heat boilers typically run 4,000–6,000 hours/year (46–68%), while peaking units may run as few as 500 hours.
- EU energy efficiency benchmarking studies report average utilization rates of 50–65% for industrial steam boilers in Western Europe (IEA, 2007; European Commission, Reference Document on Best Available Techniques for Large Combustion Plants, 2006).
- The 60% assumption implies the installation is used regularly but not as baseload — consistent with process heat demand that follows production schedules (nights, weekends, seasonal variation).

**Calculation.**

> Emissions = Capacity × CF × Hours/year × Seconds/hour × EF

> = 20 MW × 0.60 × 8,760 h/yr × 3.6 GJ/MWh × 56.1 t CO₂/TJ × 10⁻³ TJ/GJ

> = 20 × 0.60 × 8,760 × 3.6 × 56.1 / 10⁶ kt

> ≈ 21.2 kt CO₂/year

We round up to **30 kt CO₂/year** to provide a margin for installations that operate at higher capacity factors (up to ~80%), use small amounts of liquid fuel alongside gas, or operate multiple sub-20MW units at the same site.

**Reference table.** Annual CO₂ emissions (kt) at 20 MW rated thermal input by fuel and capacity factor:

|  | 40% | 50% | 60% | 70% | 80% |
|---|---|---|---|---|---|
| Natural gas (56.1 t/TJ) | 14.2 | 17.7 | **21.2** | 24.8 | 28.3 |
| Gas oil/diesel (74.1 t/TJ) | 18.7 | 23.4 | 28.0 | 32.7 | 37.4 |
| Heavy fuel oil (77.4 t/TJ) | 19.5 | 24.4 | 29.3 | 34.2 | 39.1 |
| Coal (95.0 t/TJ) | 24.0 | 30.0 | 36.0 | 41.9 | 47.9 |

The 30 kt cap corresponds to natural gas at ~80% capacity factor, gas oil at ~60%, or coal at ~50%. Any non-ETS firm exceeding 30 kt would need to either operate at a larger capacity than 20 MW (violating the ETS threshold), burn a higher-carbon fuel at high capacity (unlikely in Belgium's gas-dominated context), or both. A firm at exactly 20 MW burning natural gas at the assumed 60% capacity factor emits ~21.2 kt; the margin between 21.2 and 30 kt accommodates uncertainty in capacity factors and fuel mix.

### 9.3 Application in the Pareto allocation

The 30 kt cap replaces the data-driven min(ETS) constraint. In the iterative allocation loop (section 5.5), the upper-bound check becomes:

> if max(emissions_deployment) ≥ 30,000: cap and redistribute

This is applied uniformly across all CRF groups and years, rather than varying by the composition of the ETS firms in each cell. The cap does not depend on which ETS firms happen to be present, eliminating the pathological cases where a single small ETS firm (e.g., 24 tonnes in mfg_other) constrains thousands of deployment firms.

### 9.4 NACE 38 assignment to CRF energy

The energy CRF group (1.A.1) maps to NACE 35 (electricity, gas, steam supply) and **NACE 38** (waste collection, treatment, disposal, recovery). The inclusion of NACE 38 follows from the Belgian NID's treatment of waste incineration: since 2005 in the Flemish region and 2006 in Wallonia, all municipal waste incineration plants produce electricity and/or useful heat, and their combustion emissions are allocated to CRF category 1.A.1.a (public electricity and heat production) per IPCC guidelines (NID Section 3.2.6).

**This is a strong and problematic assumption.** Only a small fraction of NACE 38 firms are waste-to-energy incinerators — the IMJV data identifies roughly 5–9 such firms in Flanders (see `analysis/IMJV_README.md`), plus approximately 4 in Wallonia. The remaining ~200 NACE 38 deployment firms per year are regular waste collection, recycling, and treatment businesses whose combustion emissions (heating, vehicle fleets) belong to CRF 1.A.4 (commercial/institutional) or CRF 5C (waste), not to 1.A.1.a.

We adopt this mapping despite its imprecision because:

1. **The alternative is worse.** Without NACE 38, the energy CRF group has only ~37–113 NACE 35 deployment firms absorbing 4–7 million tonnes of E_deploy — an average of 80–185 kt per firm, far exceeding what any sub-20MW installation can plausibly emit. This makes the Pareto allocation infeasible in every year.

2. **With NACE 38, per-firm burden drops to 16–31 kt.** Adding ~200 NACE 38 firms brings deployment firm counts to 212–258, and per-firm allocation to within the physics-based cap. The Pareto shape then concentrates emissions on the top-ranked firms (which includes the actual waste-to-energy incinerators, since they have large B2B transaction volumes and high proxy values), while the many small waste firms at the bottom of the ranking receive near-zero emissions.

3. **The Pareto allocation self-corrects partially.** Under the GPA redistribution, firms are ranked by their proxy value. The large waste incinerators — which are the NACE 38 firms that genuinely belong in 1.A.1.a — will have high proxy values (large B2B transaction volumes with energy-related suppliers). Regular waste collection firms will have low proxy values and receive minimal emissions. The ranking mechanism thus approximates the correct assignment even though the CRF mapping treats all NACE 38 firms identically.

4. **The firm-level activity data needed to resolve this properly is not available.** Identifying which specific NACE 38 firms operate waste-to-energy installations would require matching IMJV enterprise numbers to anonymized VAT codes (possible only on RMD, and only for Flemish firms) or obtaining plant-level environmental permit data. Neither is currently feasible at scale.

This assumption should be revisited if plant-level activity data becomes available, or if the resulting emission distributions for NACE 38 firms appear implausible in downstream analysis.

**Effect on feasibility.** With the 30 kt cap and NACE 38 included in the energy CRF group, all mixed sector-years become feasible across all CRF groups for 2009–2021. The years 2005–2008 require the additional pre_ets backcast adjustment (section 8) to bring per-firm burden below 30 kt.

## 10. Scripts

| Script | Purpose | Location |
|---|---|---|
| `diagnostic_lmoment_sectors.R` | L-moment ratio diagnostic | `analysis/active/` |
| `diagnostic_redistribution_comparison.R` | Sinh vs quantile mapping vs Pareto comparison | `analysis/active/` |
| `diagnostic_pareto_threshold.R` | LOSO threshold tuning on mixed sectors | `analysis/active/` |
| `table_pareto_all_models.R` | Main results + zero-emitter tables, Pareto, no threshold | `figures_tables/` |
| `table_pareto_all_models_cvthresh.R` | Main results + zero-emitter tables, Pareto + CV threshold | `figures_tables/` |
| `table_pareto_redistribution.R` | Mixed table: sinh vs Pareto for EN, plus all models | `figures_tables/` |
| `figure_lmoment_diagnostic.R` | Publication-quality L-moment ratio diagram + table | `figures_tables/` |

## 11. Data

| File | Contents |
|---|---|
| `repeated_cv_proxy_sector_asinh.RData` | proxy_matrix (26608 × 200), panel, syt |
| `firm_year_panel_with_proxies.RData` | Revenue, NACE proxy (proxy_tabachova) |
| `training_sample.RData` | Full training panel with emissions |
