# Thoughts

---

## Uncertainty Propagation for Imputed Emissions

*Status: draft — needs further iteration before implementation*

### Context

The `facts-emissions-across-network` project computes distributional facts (dispersion, correlations, concentration indices) about GHG emissions across Belgium's production network. Observed emissions cover only EU ETS (~257 firms). For the vast majority of firms, emissions are imputed using the elastic-net B2B proxy from `inferring_emissions/`. The imputation introduces two independent sources of uncertainty that must be propagated to downstream statistics so that facts are reported with confidence intervals.

**Note**: IMJV data usage is currently uncertain. The plan proceeds without IMJV; if later included, it slots in as another source of observed (non-imputed) emissions.

### Method Overview

1. **Source 1 (coefficient instability)**: Run B = 50 subsample draws of EN estimation. Average the B proxy vectors into a single "best" proxy. Report cross-B standard deviation as a robustness metric (not folded into main CIs).
2. **Source 2 (prediction error)**: Apply K = 200 perturbations to the averaged proxy — extensive margin (proxy-value-dependent misclassification) + intensive margin (Dirichlet share perturbation). CIs come from these K realizations.

The two sources are independent and capture different things: Source 1 is model uncertainty (which coefficients?); Source 2 is prediction uncertainty conditional on the model (even with perfect coefficients, the proxy is noisy).

### Stage A: Repeated Subsampling & Averaging (Source 1) — runs on RMD

For b = 1, ..., B (B = 50):

1. **Subsample**: Leave-20%-out of training firms, stratified by sector and emitter status.
   - Why not bootstrap: duplicating firms distorts EN penalty path (cf. stability selection, Meinshausen & Buhlmann 2010).
2. **Estimate**: Elastic net on subsample (alpha=0.5, lambda via 10-fold inner CV grouped at firm level).
3. **Proxy**: Compute proxy for all deployment firms = coefficient-weighted sum of purchases from selected suppliers, floored at 0.

**Averaging**: `proxy_avg_i = (1/B) Σ_b proxy_i^(b)` for each firm-year i.

**Output**:
- `proxy_avg`: firm-year level averaged proxy, used for point estimates and Source 2 perturbation
- `proxy_sd`: cross-B standard deviation, used as robustness metric
- Individual `proxy_b` vectors (optional, for decomposition analysis)

**Classification from averaged proxy**: emitter if `proxy_avg > 0`.

**Calibration from averaged proxy**: For each sector-year (s,y):
- `E_sy_deploy = E_sy_NIR − observed ETS emissions in (s,y)`
- Distribute `E_sy_deploy` among predicted emitters proportionally to `sinh(proxy_avg_i)`

**Robustness reporting**: For each downstream statistic θ, also compute θ on each individual `proxy_b` and report the cross-B s.d. This shows sensitivity to coefficient instability.

### Stage B: Prediction Error Perturbation (Source 2) — runs on local 1

Apply K = 200 perturbations to the averaged proxy.

#### B.1 Extensive margin: proxy-value-dependent misclassification

**Key idea**: P(true emitter | proxy = π) should be increasing in π. Firms with high proxy are almost certainly true emitters; firms with proxy near 0 are more likely false positives.

**Calibration from CV data** (one-time):

1. From the M = 200 CV repeats in `inferring_emissions/`, collect all held-out firm observations with their (proxy value, true label) pairs.
2. Bin firms by proxy value (e.g., deciles among proxy > 0, plus the proxy = 0 bin).
3. Compute: `P(true emitter | proxy in bin)` = share of true emitters in each bin.
4. Smooth across bins via logistic regression of true label on log(1 + proxy): `p(π) = P(true emitter | proxy = π)`.
5. For proxy = 0 firms: `p(0) = 1 − FOR ≈ 0.006`.

**Note**: This uses training-sample prevalence (~12% emitters) as a basis for FDR/FOR. In deployment, the true emitter share is unknown and likely lower, so FDR would be even higher. This is the simplest defensible approach. The more accurate alternative (described in the next subsection) conditions on proxy value directly and avoids the prevalence problem, but the calibration still depends on training data.

**Perturbation for each draw k**:

1. Perturb logistic regression coefficients by drawing from their estimated s.e. → produces K slightly different `p_k(π)` functions.
2. For each firm i:
   - If `proxy_avg_i > 0`: draw `emitter_ik ~ Bernoulli(p_k(proxy_avg_i))`
   - If `proxy_avg_i = 0`: draw `emitter_ik ~ Bernoulli(FOR_k)`, where `FOR_k ~ N(0.006, σ²)` truncated to [0,1]
3. Firms flipped emitter → non-emitter: set emissions to 0
4. Firms flipped non-emitter → emitter: assign proxy from lower quartile of within-sector emitter proxy distribution

**Why not flat FDR/FOR rates**: Using a flat rate ignores the information in the proxy value. A firm with proxy = 1,000 is almost certainly a true emitter; a firm with proxy = 0.001 is plausibly a false positive. The logistic calibration exploits this signal and is more defensible to referees.

#### B.2 Intensive margin: Dirichlet share perturbation

Among predicted emitters in each sector-year cell (s,y) after extensive-margin perturbation:

1. Compute shares: `w_i = sinh(proxy_avg_i) / Σ_j sinh(proxy_avg_j)`
2. Draw `ρ_k ~ N(0.639, 0.054²)` truncated to [0,1] (target rank correlation from CV, mean and s.d. across M=200 repeats)
3. Look up `α = α(N_s, ρ_k)` from pre-computed table
4. Draw perturbed shares: `w* ~ Dirichlet(α · w)`
5. Compute perturbed emissions: `ê_i = E_sy_deploy × w*_i`

**Why Dirichlet**: Shares positive and sum to 1 by construction → NIR sector-year total preserved exactly. One parameter (α) calibrated to one summary statistic (rank correlation). For N_s = 1: no perturbation.

**Concentration parameter α calibration** (one-time):
- For each (N_s, ρ_target) on a grid, bisect α until expected Spearman(original ranks, Dirichlet-perturbed ranks) ≈ ρ_target (averaged over 5,000 simulation draws).
- Store as lookup table.

#### B.3 Recalibration

`E_sy_deploy` distributed among the perturbed emitter set proportionally to perturbed shares. Sector-year total preserved exactly.

### Stage C: Downstream Statistics — runs on local 1

For each of K = 200 realizations, merge observed ETS emissions (fixed) + perturbed imputed emissions and compute:

- **RQ1**: Within-sector dispersion of scope 1 (variance, Gini, percentile ratios)
- **RQ2**: Within-sector dispersion of upstream network-adjusted emissions
- **RQ3**: Within-sector correlation between scope 1 and upstream
- **RQ4**: Within-sector rank correlation between scope 1 and upstream
- **RQ5**: HHI of emission exposure across suppliers
- **RQ6**: Network depth of embodied emissions

Compute statistics on-the-fly within the perturbation loop (avoid storing 200 full emission vectors).

### Stage D: Confidence Intervals & Reporting

**From Source 2 (K = 200)**:
- Point estimate: statistic computed on unperturbed averaged proxy (or median across K)
- 90% CI: [5th, 95th percentile] across K realizations
- IQR for visual summaries

**From Source 1 (B = 50, reported separately)**:
- Cross-B s.d. of each statistic
- Fraction of B draws where qualitative finding holds (e.g., "upstream > scope 1 dispersion in 48/50 draws")

**In the paper**:
- Main tables: point estimates with 90% CIs from Source 2
- Robustness section/appendix: Source 1 variation
- Figures: shaded bands for IQR

### Implementation Steps (to be done)

1. Calibrate proxy-value-dependent error rates → `analysis/calibrate_extensive_margin.R`
2. Calibrate Dirichlet α lookup table → `analysis/calibrate_dirichlet_alpha.R`
3. Write Stage A subsampling script → runs on RMD
4. Transfer Stage A outputs to local 1
5. Write Stage B + C perturbation + statistics script → `analysis/uncertainty_propagation.R`
6. Write Stage D reporting script → `analysis/report_uncertainty.R`
