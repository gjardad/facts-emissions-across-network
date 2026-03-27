# Retransformation Bias in the Emission Imputation Pipeline

## The problem

Our elastic net (EN) predicts in asinh space:

$$
\text{proxy}_i \approx \mathbb{E}[\text{asinh}(y_i) \mid X_i]
$$

where $y_i$ is firm $i$'s emissions. To recover emissions in levels, we apply the inverse transformation $\sinh(\cdot)$ and use the result as weights to redistribute sector-year totals:

$$
w_i = \frac{\sinh(\text{proxy}_i)}{\sum_{j \in s} \sinh(\text{proxy}_j)}, \quad \hat{y}_i = E_{s,t}^{\text{deploy}} \times w_i
$$

The diagnostic shows that this pipeline **systematically underestimates emissions for large emitters** and overestimates for small emitters. Within sectors, 26.7% of firms receive imputed emissions below 10% of their true emissions (mostly large emitters), while only 2.1% receive more than 10× their true emissions. The median imputed/actual ratio falls from 0.94 in the bottom emission quartile to 0.21 in the top quartile.

This inflates within-sector carbon productivity dispersion by ~18% on the 90-10 log gap and ~48% on the median p90/p10 ratio relative to actual ETS emissions.

## Jensen's inequality and the retransformation bias

### The core result

Let $Z = \text{asinh}(y)$ and suppose we observe $\hat{Z} = \mathbb{E}[Z \mid X]$ (the EN prediction). We want $\mathbb{E}[y \mid X]$, but we compute $\sinh(\hat{Z})$ instead. Since $\sinh(\cdot)$ is a **strictly convex** function (its second derivative $\sinh''(z) = \sinh(z) > 0$ for $z > 0$), Jensen's inequality gives:

$$
\sinh\big(\mathbb{E}[Z \mid X]\big) < \mathbb{E}\big[\sinh(Z) \mid X\big] = \mathbb{E}[y \mid X]
$$

That is, $\sinh(\text{proxy})$ is a **downward-biased** estimate of $\mathbb{E}[y \mid X]$. The bias is larger when:

1. **The conditional variance of $Z$ given $X$ is larger** — more prediction uncertainty means more room for Jensen's inequality to bite.
2. **The level of $Z$ is higher** — $\sinh$ is more convex at larger values (its curvature $\sinh''(z)/\sinh(z) = 1$ is constant in relative terms, but the absolute curvature $\sinh''(z) = \sinh(z)$ grows exponentially), so the bias in levels is larger for firms with high true emissions.

### Intuition

Suppose the EN predicts $\hat{Z}_i = 10$ for a large emitter, but the true $Z_i$ is equally likely to be 9 or 11 (prediction error $\varepsilon \sim \pm 1$). Then:

- $\sinh(10) = 11{,}013$
- $\frac{1}{2}[\sinh(9) + \sinh(11)] = \frac{1}{2}[4{,}052 + 29{,}937] = 16{,}995$

So $\sinh(\hat{Z})$ gives 11,013 but the true expected emission is 16,995 — an underestimate of 35%. The asymmetry comes from sinh being steeper above the prediction than below it. Overshooting by 1 unit in asinh space adds much more in levels than undershooting by 1 unit subtracts.

For a small emitter with $\hat{Z}_i = 2$ and the same $\varepsilon \sim \pm 1$:

- $\sinh(2) = 3.63$
- $\frac{1}{2}[\sinh(1) + \sinh(3)] = \frac{1}{2}[1.18 + 10.02] = 5.60$

The underestimate is 35% again in relative terms, but in absolute levels it's 2 tonnes vs. 6,000 tonnes. When these firms compete for shares of a fixed sector-year total, the large emitter's share gets compressed.

## The log-transformation analog in econometrics

This problem is well-known in the econometrics of log-linear models. If we estimate:

$$
\log(y_i) = X_i \beta + \varepsilon_i
$$

and predict $\hat{y}_i = \exp(X_i \hat{\beta})$, we get a biased estimate of $\mathbb{E}[y_i \mid X_i]$ for the same reason: $\exp(\cdot)$ is convex, so $\exp(\mathbb{E}[\log y]) < \mathbb{E}[\exp(\log y)] = \mathbb{E}[y]$.

Under normality $\varepsilon_i \sim N(0, \sigma^2)$, the correct conditional mean is:

$$
\mathbb{E}[y_i \mid X_i] = \exp\left(X_i \beta + \frac{\sigma^2}{2}\right)
$$

The correction factor $\exp(\sigma^2/2)$ is always greater than 1. For our asinh case, the analogous result (derived below) is $\sinh(\hat{Z}) \times \exp(\sigma^2/2)$.

This problem was first highlighted by Goldberger (1968) and solutions were proposed by Duan (1983), Manning (1998), and others. It is sometimes called the "retransformation problem" or "smearing problem."

### References in econometrics

- Goldberger, A.S. (1968). "The Interpretation and Estimation of Cobb-Douglas Functions." *Econometrica*, 36(3-4), 464-472.
- Duan, N. (1983). "Smearing Estimate: A Nonparametric Retransformation Method." *JASA*, 78(383), 605-610.
- Manning, W.G. (1998). "The Logged Dependent Variable, Heteroscedasticity, and the Retransformation Problem." *Journal of Health Economics*, 17(3), 283-295.
- Manning, W.G. and Mullahy, J. (2001). "Estimating Log Models: To Transform or Not to Transform?" *Journal of Health Economics*, 20(4), 461-494.

## The Duan smearing estimator

### Setup

Suppose the model is $Z_i = g(X_i; \beta) + \varepsilon_i$ where $Z = \text{asinh}(y)$ and $\varepsilon_i$ are prediction errors with unknown distribution. We want $\mathbb{E}[y_i \mid X_i] = \mathbb{E}[\sinh(Z_i) \mid X_i]$.

### Parametric correction (under normality)

If $\varepsilon_i \sim N(0, \sigma^2)$, we can derive the correction analytically using the identity $\sinh(a + b) = \sinh(a)\cosh(b) + \cosh(a)\sinh(b)$ and the moment generating function of the normal distribution:

$$
\mathbb{E}[\sinh(\hat{Z}_i + \varepsilon_i)] = \sinh(\hat{Z}_i) \cdot \mathbb{E}[\cosh(\varepsilon_i)] + \cosh(\hat{Z}_i) \cdot \underbrace{\mathbb{E}[\sinh(\varepsilon_i)]}_{= 0 \text{ by symmetry}}
$$

Since $\cosh(\varepsilon) = \frac{1}{2}(e^\varepsilon + e^{-\varepsilon})$ and $\varepsilon \sim N(0, \sigma^2)$:

$$
\mathbb{E}[\cosh(\varepsilon)] = \frac{1}{2}\left(e^{\sigma^2/2} + e^{\sigma^2/2}\right) = e^{\sigma^2/2}
$$

Therefore:

$$
\mathbb{E}[y_i \mid X_i] = \sinh(\hat{Z}_i) \cdot e^{\sigma^2/2}
$$

This is the **asinh analog of the log-retransformation correction**. The correction factor $e^{\sigma^2/2}$ is the same as in the log case — it depends only on the variance of the prediction error in the transformed space.

### Nonparametric correction (Duan smearing)

If we don't want to assume normality, Duan (1983) proposed using the empirical distribution of residuals directly:

$$
\hat{\mathbb{E}}[y_i \mid X_i] = \frac{1}{n} \sum_{k=1}^{n} \sinh(\hat{Z}_i + \hat{\varepsilon}_k)
$$

where $\hat{\varepsilon}_k$ are the OOS residuals (from cross-validation). This averages over the empirical residual distribution and is consistent regardless of the error distribution. We have 200 CV repeats × ~3,000 ETS firm-years of residuals available.

## Why it depends on homo- vs. heteroscedasticity

### Homoscedastic case: correction cancels in shares

If $\sigma^2$ is the same for all firms within a sector-year, then the corrected weight for firm $i$ is:

$$
w_i^{\text{corrected}} = \frac{\sinh(\hat{Z}_i) \cdot e^{\sigma^2/2}}{\sum_j \sinh(\hat{Z}_j) \cdot e^{\sigma^2/2}} = \frac{\sinh(\hat{Z}_i)}{\sum_j \sinh(\hat{Z}_j)} = w_i
$$

The correction factor cancels. **Under homoscedasticity, the Duan smearing estimator does not change the redistribution shares.** The sector-year total calibration already absorbs the level correction, and the relative shares remain unchanged.

This means: if the EN has constant prediction variance across all firms within a sector-year, the retransformation bias cannot be fixed by smearing. The compression of large emitters would persist.

### Heteroscedastic case: correction reshapes shares

If prediction variance varies across firms — e.g., $\sigma^2_i$ is larger for firms with higher true emissions or for firms with certain characteristics — then:

$$
w_i^{\text{corrected}} = \frac{\sinh(\hat{Z}_i) \cdot e^{\sigma^2_i/2}}{\sum_j \sinh(\hat{Z}_j) \cdot e^{\sigma^2_j/2}}
$$

Now the correction factors **do not cancel**. Firms with larger $\sigma^2_i$ get a proportionally larger upward correction. If large emitters have more prediction uncertainty (plausible: extreme values are harder to predict), then $\sigma^2_i$ is larger for large $\hat{Z}_i$, and the correction pushes their shares up — counteracting the compression we observe.

### Empirical test

We can test for heteroscedasticity directly from the CV residuals. For each firm-year in the training sample, we have 200 OOS residuals $\hat{\varepsilon}_{i,r} = Z_i - \hat{Z}_{i,r}$. Compute $\hat{\sigma}^2_i = \text{Var}_r(\hat{\varepsilon}_{i,r})$ across repeats. Then check whether $\hat{\sigma}^2_i$ correlates with $|Z_i|$ (emission size), sector, or other firm characteristics.

If heteroscedasticity is present, the firm-specific Duan correction is:

$$
\hat{\mathbb{E}}[y_i \mid X_i] = \frac{1}{R} \sum_{r=1}^{R} \sinh(\hat{Z}_i + \hat{\varepsilon}_{i,r})
$$

where $\hat{\varepsilon}_{i,r}$ are firm $i$'s own residuals across the $R$ CV repeats (or residuals from firms with similar characteristics if firm-specific residuals are too noisy).

## Alternative: GLM with log link (avoiding retransformation entirely)

### Motivation

The retransformation problem arises because we estimate in a transformed space and need to invert the transformation. An alternative is to model the conditional mean of $y$ directly in levels, using a link function that ensures non-negativity without requiring back-transformation.

### Specification

A Poisson or Gamma GLM with log link models:

$$
\log\left(\mathbb{E}[y_i \mid X_i]\right) = X_i \beta
$$

or equivalently:

$$
\mathbb{E}[y_i \mid X_i] = \exp(X_i \beta)
$$

The key difference from the asinh approach:

- **asinh approach**: models $\mathbb{E}[\text{asinh}(y_i) \mid X_i] = X_i \beta$, then inverts with sinh.
- **GLM approach**: models $\mathbb{E}[y_i \mid X_i]$ directly through a log link. No inversion needed.

The predicted $\hat{y}_i = \exp(X_i \hat{\beta})$ is an estimate of the **conditional mean** $\mathbb{E}[y_i \mid X_i]$ — not a point prediction that needs retransformation. This is because the GLM is specified in terms of the mean of $y$, not the mean of a transformation of $y$.

### Poisson quasi-likelihood

The most common implementation is **Poisson quasi-maximum likelihood (PQML)**, which maximizes:

$$
\sum_i \left[ y_i \cdot X_i \beta - \exp(X_i \beta) \right]
$$

This estimator is consistent for $\beta$ as long as the conditional mean is correctly specified ($\mathbb{E}[y | X] = \exp(X\beta)$), regardless of the true distribution of $y$. It does **not** require $y$ to be Poisson-distributed or even integer-valued. The Poisson family is used only for the score equations — the estimator is robust to distributional misspecification (Santos Silva and Tenreyro, 2006).

This is the same estimator used in the gravity model literature for trade flows, where it is standard precisely because trade data is highly skewed and log-linearization introduces retransformation bias.

### Implementation with glmnet

`glmnet` supports `family = "poisson"` with a log link natively:

```r
fit <- cv.glmnet(
  x = X_full,
  y = emissions,              # in levels, not asinh
  family = "poisson",         # log link, Poisson quasi-likelihood
  alpha = 0.5,                # elastic net mixing
  penalty.factor = pf,
  foldid = inner_foldid,
  standardize = TRUE
)

# Predicted conditional mean — no back-transformation needed
yhat <- predict(fit, newx = X_deploy, s = "lambda.min", type = "response")
```

The `type = "response"` option returns $\exp(X\hat{\beta})$ directly. These predictions can be used as redistribution weights without any correction.

### Tradeoffs relative to asinh + EN

| | asinh + Gaussian EN | Poisson PQML + EN |
|---|---|---|
| **Loss function** | Penalizes absolute errors in asinh space — approximately equal weight to proportional errors across the emission distribution | Penalizes relative errors in levels — naturally downweights large emitters' absolute errors |
| **Retransformation** | Required (sinh). Introduces bias unless corrected | Not required. $\exp(X\hat{\beta})$ estimates $\mathbb{E}[y \mid X]$ directly |
| **Zeros** | Handled naturally (asinh(0) = 0) | Must be handled carefully — log(0) is undefined. Standard approach: include zeros in the Poisson likelihood (Poisson allows $y = 0$ natively) |
| **Sparsity** | EN zeroes out irrelevant supplier coefficients in asinh space | EN zeroes out suppliers in log-mean space |
| **Prediction at tails** | sinh back-transformation amplifies errors for large emitters | exp back-transformation can overshoot (exp is also convex), but no retransformation bias since we model $\mathbb{E}[y]$ directly |
| **Existing code** | Current pipeline | Requires changing `run_subsampled_en.R` and re-estimating |

### References

- Santos Silva, J.M.C. and Tenreyro, S. (2006). "The Log of Gravity." *Review of Economics and Statistics*, 88(4), 641-658.
- Wooldridge, J.M. (2010). *Econometric Analysis of Cross Section and Panel Data*, 2nd edition. MIT Press. Chapter 18.
- Gourieroux, C., Monfort, A., and Trognon, A. (1984). "Pseudo Maximum Likelihood Methods: Theory." *Econometrica*, 52(3), 681-700.

## Summary of options

| Approach | Changes required | Principled? | Addresses the bias? |
|---|---|---|---|
| **Duan smearing (heteroscedastic)** | Post-processing only; uses existing CV residuals | Yes — standard econometric correction | Only if prediction variance varies across firms |
| **GLM with log link** | Re-estimate EN with `family = "poisson"` | Yes — avoids retransformation entirely | Yes, by construction |
| **Distributional calibration** | Add variance-matching constraint to calibration | Partially — requires assumption about deployment firms | Mechanically, but relies on ETS-deployment similarity |

The **first diagnostic step** is to check whether the CV residuals are heteroscedastic. If $\text{Var}(\varepsilon_i)$ correlates with $|\hat{Z}_i|$ or with firm size, the Duan correction with firm-specific variance will reshape the within-sector shares and reduce the compression of large emitters — without changing the EN model at all.
