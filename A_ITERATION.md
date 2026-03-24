# Neumann Series Convergence and the Leontief A Matrix

## What the Neumann series computes

We compute upstream embodied emissions via the Neumann series approximation to the Leontief inverse:

```
m = (I - A)^{-1} eps = eps + A*eps + A^2*eps + A^3*eps + ...
```

where `A` is the domestic input-output matrix and `eps` is the emission intensity vector (emissions per unit cost). Each term `A^k * eps` captures emissions embodied k layers upstream in the domestic supply chain.

## How A is defined

The entry `A_{ij} = purchases_{ij} / cost_i`, where:

```
cost_i = wage_bill_i + domestic_B2B_inputs_i + total_imports_i + permit_cost_i
```

By construction, `domestic_B2B_inputs_i` is a strict subset of `cost_i`, so row sums of A are strictly less than 1. This guarantees the Neumann series converges.

## Convergence criterion

We use `NEUMANN_MAXIT = 50` and `NEUMANN_TOL = 1e-8`. After each iteration k, we check:

```
max(|term_k|) / (max(|m|) + 1e-15) < NEUMANN_TOL
```

If this is satisfied before iteration 50, the series is deemed converged. Otherwise, we stop at 50 iterations.

## Why some years do not formally converge at K=50

Convergence speed is governed by the spectral radius of A, which is bounded above by the max row sum. When the max row sum is close to 1, convergence is slow: the remaining signal decays as ~(max_row_sum)^K per iteration.

The firms driving high row sums are **domestic trade intermediaries** — wholesale (NACE 46xxx) and retail (NACE 47xxx) firms. For these firms, domestic B2B purchases constitute nearly 100% of their observable cost because they add very little value: thin margins, low wage bills relative to turnover, no imports, no manufacturing.

## What non-convergence means in practice

The Neumann series propagates emission intensities upstream through the domestic supply chain. Each term A^k * eps represents the emission content of inputs k layers back. When firms trade in cycles (e.g. a wholesaler and retailer that buy from each other), the emission signal circulates through these loops, decaying by the product of the row sums at each step.

Non-convergence does NOT mean we are missing an emission source. It means that for a small number of pass-through firms, we are slightly undercounting the propagation of emissions through circular domestic trade. All emission sources (ETS verified emissions and imputed deployment emissions) are fully accounted for in the first term. Subsequent terms redistribute those emissions through the network.

## Quantifying the gap: K=50 vs K=500 (full RMD data)

We extended the series to K=500 for comparison on the full data for selected years:

| Year | N firms | Max row sum | Aggregate gap | Median firm gap | 99th pctile gap | Max firm gap | Firms >1% gap |
|------|---------|-------------|---------------|-----------------|-----------------|--------------|---------------|
| 2005 | 120,389 | 0.9999399 | 0.019% | 0.000% | 0.003% | 64.4% | 525 |
| 2010 | 126,188 | 0.9999989 | 0.000% | 0.000% | 0.003% | 7.5% | 5 |
| 2015 | 129,249 | 0.9999975 | 0.009% | 0.003% | 0.024% | 16.7% | 421 |
| 2020 | 123,436 | 0.9999996 | 0.018% | 0.012% | 0.025% | 28.8% | 7 |
| 2021 | 123,435 | 0.9999996 | 0.001% | 0.000% | 0.003% | 16.1% | 12 |

For 2005–2021, the aggregate gap from stopping at K=50 is below 0.02% in all years. The median firm-level gap is negligible. A handful of individual firms (pass-through intermediaries with circular trade) can have gaps above 1%, but these are economically insignificant in absolute terms.

## Why this is not problematic (2005–2021)

1. **Aggregate accuracy is near-perfect.** Total upstream emissions at K=50 miss <0.02% relative to K=500 in all years.

2. **The gap is concentrated in a few outlier firms.** The 99th percentile firm-level gap is below 0.03%. Firms with >1% gap number in single digits to low hundreds out of ~120,000.

3. **The non-converging signal is circular re-attribution, not missing emissions.** The unresolved remainder at K=50 is emissions being re-counted as they loop through intermediary trade cycles — not an unobserved emission source.

---

## Dropping 2022 from the panel

Year 2022 is excluded from the analysis due to a data quality issue that invalidates the cost-based A matrix.

### The problem

In 2022, 51.4% of firms (62,364 out of 121,299) have A row sums above 0.99, compared to 1.0% in 2021. Twenty-one firms have row sums effectively equal to 1.0. The Neumann convergence gap in 2022 is 0.4% in aggregate — an order of magnitude worse than any other year.

### Root cause

The 21 firms with row sum ≈ 1.0 all have `wage_bill = 1 EUR` — a placeholder value, not a real wage bill. These firms have tens of millions in domestic B2B purchases but only €1 in recorded wages. More broadly, among the 62,364 firms with row sum > 0.99 in 2022:

- 100% are present in annual accounts with positive wage_bill (no missing data)
- But median `wage_bill / domestic_B2B` ratio is 0.4% (vs. typical values well above 1% in other years)
- All 21 firms with row sum ≈ 1.0 have wage_bill ∈ {1, 2} EUR, zero imports, zero ETS permits

This pattern — mass near-zero wage bills with no corresponding drop in B2B — does not appear in 2020 or 2021. It is consistent with incomplete annual accounts filing for 2022 (the most recent year in the data), where wage_bill fields were filled with placeholder values for firms that had not yet filed.

### Comparison with 2021

| Metric | 2021 | 2022 |
|--------|------|------|
| B2B firms | 123,435 | 121,299 |
| In annual accounts | 100% | 100% |
| wage_bill > 0 | 100% | 100% |
| Firms with row sum > 0.99 | 1,225 (1.0%) | 62,364 (51.4%) |
| Firms with row sum > 0.999 | — | 7,898 |
| Neumann K=50 vs K=500 aggregate gap | 0.001% | 0.401% |

### Justification for dropping

The wage_bill data quality issue in 2022 fundamentally undermines the cost-based A matrix: if wage_bill is artificially near-zero, the cost denominator is dominated by domestic B2B, row sums approach 1, and the Leontief inverse is inflated. This affects not just convergence speed but the economic interpretation of upstream emissions. We therefore exclude 2022 from the analysis panel. The sample covers 2005–2021 (17 years).
