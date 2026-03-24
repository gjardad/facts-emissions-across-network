# Neumann Series Convergence in `compute_b_loop.R`

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

## Convergence diagnostics (local downsampled data)

| Year | Max row sum | Median K | Max K | % converged |
|------|-------------|----------|-------|-------------|
| 2005 | 0.9982 | 27.5 | 36 | 100% |
| 2006 | 0.9912 | 50.0 | 50 | 10% |
| 2007 | 0.9999 | 42.0 | 42 | 100% |
| 2008 | 0.9935 | 20.0 | 28 | 100% |
| 2009 | 0.9984 | 20.0 | 20 | 100% |
| 2010 | 0.9896 | 18.0 | 20 | 100% |
| 2011 | 0.9927 | 20.0 | 20 | 100% |
| 2012 | 0.9933 | 35.0 | 43 | 100% |
| 2013 | 0.9917 | 31.0 | 31 | 100% |
| 2014 | 0.9959 | 50.0 | 50 | 0% |
| 2015 | 0.9794 | 50.0 | 50 | 0% |
| 2016 | 0.9787 | 50.0 | 50 | 0% |
| 2017 | 0.9824 | 50.0 | 50 | 0% |
| 2018 | 0.9945 | 50.0 | 50 | 0% |
| 2019 | 0.9984 | 50.0 | 50 | 0% |
| 2020 | 0.9892 | 50.0 | 50 | 0% |
| 2021 | 0.9949 | 50.0 | 50 | 0% |
| 2022 | 0.9999 | 50.0 | 50 | 0% |

For 2014-2022, the convergence criterion `rel < 1e-8` is not met within 50 iterations.

## Why some years do not converge at K=50

Convergence speed is governed by the spectral radius of A, which is bounded above by the max row sum. When the max row sum is close to 1 (e.g. 0.9999 for 2022), convergence is very slow: the remaining signal decays as ~(max_row_sum)^K per iteration.

The firms driving high row sums are **domestic trade intermediaries** -- wholesale (NACE 46xxx) and retail (NACE 47xxx) firms. For these firms, domestic B2B purchases constitute nearly 100% of their observable cost because they add very little value: thin margins, low wage bills, no imports, no manufacturing.

Example from 2022: the firm with the highest row sum (0.9999) is a fuel retailer (NACE 47781) with domestic B2B purchases of EUR 556M and a wage bill of EUR 42K. Its row sum is 556M / (556M + 42K) = 0.9999.

Diagnostics for firms with row sum > 0.99 in 2022:
- **30 firms** out of 1,302 total (2.3%)
- **All 30** are present in annual accounts with positive wage bill -- no missing data
- Dominated by NACE 46 (wholesale) and NACE 47 (retail) sectors
- Their wage bills are 3-4 orders of magnitude smaller than their domestic B2B purchases

## What non-convergence means in practice

The Neumann series propagates emission intensities upstream through the domestic supply chain. Each term A^k * eps represents the emission content of inputs k layers back. When firms trade in cycles (e.g. a wholesaler and retailer that buy from each other), the emission signal circulates through these loops, decaying by the product of the row sums at each step.

For a firm with row sum 0.999 involved in a cycle with another firm at row sum 0.989, the signal decays as ~0.988^K per round trip. At K=50, roughly 55% of the circular component remains unresolved. But this circular component is a re-attribution of already-counted emissions through intermediary loops, not an unobserved emission source.

Non-convergence does NOT mean we are missing an emission source. It means that for a small number of pass-through firms, we are slightly undercounting the propagation of emissions through circular domestic trade.

## Quantifying the gap: K=50 vs K=500 (year 2022)

We extended the series to K=500 for comparison:

**Aggregate level:**
- Total upstream emissions at K=50: 156,516,143 tonnes
- Total upstream emissions at K=500: 156,519,085 tonnes
- Gap: 0.002%

**Firm-level distribution of relative gap (upstream_500 - upstream_50) / upstream_500:**
- Median: 0.0002%
- 90th percentile: 0.001%
- 99th percentile: 0.005%
- Maximum: 6.95% (one firm)

**High row-sum firms (row sum > 0.99):**
- Relative gap: 0% (their upstream emissions come from ETS firms that converge quickly)
- Their share of total upstream emissions: 6.9%

The worst-affected firm (6.95% gap) has an absolute gap of 2,510 tonnes. The mean gap across all firms with positive upstream emissions is 2 tonnes.

## Why this is not problematic

1. **Aggregate accuracy is near-perfect.** The total upstream emissions at K=50 miss 0.002% relative to K=500. Statistics computed at the sector-year level (Gini coefficients, percentile ratios, rank correlations) are robust to this level of approximation.

2. **The gap is concentrated in a single firm.** The 99th percentile firm-level gap is 0.005%. Only one firm has a gap above 1%.

3. **The high row-sum firms themselves have zero gap.** Their upstream emissions are well-resolved because the emission signal reaching them from ETS firms converges quickly. The slow convergence is in the circular residual, which is economically negligible.

4. **The non-converging signal is circular re-attribution, not missing emissions.** All emission sources (ETS verified emissions and imputed deployment emissions) are fully accounted for in the first term. Subsequent terms redistribute those emissions through the network. The unresolved remainder at K=50 is a small fraction of emissions being re-counted as they loop through intermediary trade cycles.
