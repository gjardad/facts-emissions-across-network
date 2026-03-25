# Domestic Carbon Leakage via the Production Network

## Research Question

Is there evidence of **domestic carbon leakage** in Belgium in response to EU ETS carbon pricing? Specifically: does output reallocate from EU ETS-regulated firms to their untreated competitors?

The first step is to define the **set of competitors** for any given EU ETS firm using B2B transaction data.

## Step 1: Defining Competitor Sets

### The Problem

We need to identify, for each EU ETS-treated firm, which untreated firms compete with it — i.e., could absorb its output if carbon costs drive customers to switch suppliers.

### Approach 1: Co-Suppliers within Same NACE 4-digit Sector (Recommended Starting Point)

Two firms are competitors if they:
1. Share the same NACE 4-digit sector code, **and**
2. Co-supply at least $X$ common customers in the B2B data.

This is transparent, easy to explain to referees, and hard to argue against. The NACE filter ensures firms produce similar goods; the co-supplier condition ensures they actually serve the same market.

Optionally, add a **regional** filter (same province or arrondissement) given the evidence from Magerman et al. that markets are geographically local in most sectors.

### Approach 2: Co-Supplier Counts (Without Sector Restriction)

Define firm $i$'s competitor set as all firms sharing at least $X$ customers (or above some percentile of the co-supplier distribution). This is more flexible — it can capture cross-sector substitution — but harder to defend as genuinely identifying substitutes vs. complements (see below).

### Approach 3: Firm Embeddings (Magerman, Carvalho & Hansen)

A more sophisticated ML approach that learns firm similarity from the full network topology. See [Appendix: Firm Embeddings](#appendix-firm-embeddings-magerman-carvalho--hansen) for the full model.

**Assessment:** Likely overkill for this application. The embeddings are a powerful general-purpose tool, but they introduce researcher degrees of freedom (embedding dimension, clustering method, similarity threshold) and are computationally expensive (50K epochs on 7.3M observations). The simpler approaches above capture the same core idea — firms that serve the same buyers — in a more transparent way. The embeddings are best suited as a **robustness check** rather than the primary market definition.

## Step 2: Distinguishing Substitutes from Complements

A key challenge: two firms that co-supply the same buyer could be **substitutes** (both sell steel to a manufacturer) or **complements** (one sells steel, the other sells packaging). Only substitutes matter for leakage.

### Strategy A: Same Narrow Sector + Co-Suppliers

If two firms are in the same NACE 4-digit sector *and* share many customers, they are much more likely to be substitutes. This is the simplest and most defensible filter — it's the reason Approach 1 is recommended.

### Strategy B: Time-Series Substitution Patterns

Substitutes should exhibit **negative** correlation in their sales to a shared customer over time — when one gains, the other loses. Complements should show **positive** correlation (the customer buys both or neither).

With B2B data from 2002-2022, this can be tested at the buyer-supplier-pair level. For each shared customer $j$ and co-supplier pair $(i, i')$, compute:

$$\text{Corr}\left(\Delta \ln m_{ijt}, \, \Delta \ln m_{i'jt}\right)$$

where $m_{ijt}$ is the sales value from supplier $i$ to customer $j$ in year $t$. Negative correlation $\Rightarrow$ substitutes; positive $\Rightarrow$ complements.

This could be a useful empirical contribution on its own and provides a data-driven way to validate the competitor definitions.

### Strategy C: Product Overlap via Prodcom

For manufacturing firms covered by the Prodcom survey, check whether two co-suppliers produce similar PC8 products. If they do, they are substitutes. Limited to the Prodcom sample (~200k firms, NACE 07-33).

## Timing Considerations

Competitor sets should be defined from **pre-treatment** data (before EU ETS, so 2002-2004) to avoid endogeneity — if the network has already adjusted to the policy, we'd be defining competitors from an outcome of the treatment. B2B data is available from 2002.

## Data Requirements

All data is available on RMD:
- **B2B** (`B2B_ANO.dta`): buyer-supplier-year transactions, 2002-2022
- **Annual Accounts** (`Annual_Accounts_MASTER_ANO.dta`): firm characteristics including NACE codes
- **EUTL Belgium**: EU ETS installation-to-firm mapping

---

## Appendix: Firm Embeddings (Magerman, Carvalho & Hansen)

**Paper:** "Finding Markets via Firm Embeddings: A machine learning approach to relevant markets"
**Source:** INET Oxford seminar, May 2024 (work in progress)

### Motivation

The goal is to identify relevant markets using firm-to-firm transaction data. The paper develops **firm embeddings**: low-dimensional vector representations of firms, learned purely from the topology of the production network. The core intuition is that **two suppliers have similar embeddings if they tend to co-supply the same customers** — analogous to word2vec.

### Model Setup

- Set of producers: $i = 1, \ldots, N$.
- Sales relationship from supplier $i$ to customer $j$: $a_{ij} \in \{0, 1\}$.
- Adjacency matrix of the production network: $a_{ij} \in G = (N \times N)$.

**Important:** the algorithm receives only the binary network — no firm characteristics, sector, or trade values.

### Probabilistic Model

The conditional probability that firm $i$ supplies customer $j$, given all other suppliers to $j$, is Bernoulli:

$$i \in S_j \mid S_j^{-i} \sim \text{Bernoulli}(p_i), \quad \text{where } p_i \equiv \Pr\left(i \in S_j \mid S_j^{-i}\right)$$

The probability of the full supplier set $S_j$ for a given customer $j$:

$$P(S_j) = \prod_{i \in S_j} \Pr\left(i \in S_j \mid S_j^{-i}\right)$$

Aggregating over all customers $j \in C$:

$$P(G) = \prod_{j \in C} \prod_{i \in S_j} \Pr\left(i \in S_j \mid S_j^{-i}\right)$$

### Embedding Parameterization

Each firm $i$ is assigned two embedding vectors:

- $\rho_i \in \mathbb{R}^K$ — the firm's **supplier** embedding
- $\alpha_i \in \mathbb{R}^K$ — the firm's **context** (co-supplier) embedding

The conditional probability is parameterized via a sigmoid link function:

$$\Pr\left(i \in S_j \mid S_j^{-i}\right) = \sigma\left(\frac{1}{\#S_j^{-i}} \, \rho_i^T \sum_{s \in S_j^{-i}} \alpha_s\right)$$

The probability of firm $g$ **not** supplying customer $j$:

$$\Pr(g \notin S_j \mid S_j) = 1 - \sigma\left(\frac{1}{\#S_j} \, \rho_g^T \sum_{s \in S_j} \alpha_s\right)$$

### Objective Function

$$\mathcal{L}(\rho, \alpha) = \underbrace{\sum_{j \in C} \sum_{i \in S_j} \log \sigma\left(\frac{1}{\#S_j^{-i}} \, \rho_i^T \sum_{s \in S_j^{-i}} \alpha_s\right)}_{\text{conditional probability of supplying } j} + \gamma \underbrace{\sum_{j \in C} \sum_{k \in NS} \log\left(1 - \sigma\left(\frac{1}{\#S_j} \, \rho_k^T \sum_{s \in S_j^{-i}} \alpha_s\right)\right)}_{\text{negative sampling}} + \underbrace{\log p(\rho) + \log p(\alpha)}_{\text{priors on embedding vectors}}$$

**Estimation:**

- Minimize the negative of $\mathcal{L}$ to obtain $\rho_i$ and $\alpha_i$ for all $i \in N$.
- Hyperparameters: $K = 50$ (embedding dimension), epochs $= 50{,}000$, $NS = 100$ (negative samples).
- Gaussian prior constrains embeddings to be non-negative.

### Validation: Embeddings and Geography

For all firms $i$ and $j$ in the same NACE 4-digit sector:

$$\cos \text{sim}_{ij} = \beta_0 + \beta_1 \ln \text{distance}_{ij} + \varepsilon_{ij} \quad \forall \, i, j \in \text{NACE4}_K$$

| | $\ln \text{distance}_{i,j}$ |
|---|---|
| % of sectors with $\beta_1 < 0$ | 82.6% |
| % of sectors with $\beta_1 > 0$ | 0.9% |
| % of sectors with $\beta_1 = 0$ at 95% CI | 16.5% |
| N | 430 |

In 82.6% of sectors, similarity between suppliers decreases with distance — markets are geographically local.

### Data (Belgium, 2014)

| Source | Content |
|---|---|
| **NBB B2B Transactions** (Dhyne, Magerman & Rubinova, 2015) | Universe of firm-to-firm transactions. 14M links. Sales value $m_{ijt}$. |
| **Annual Accounts** | Sales, inputs, employment, labor cost, capital, ... |
| **Crossroads Bank of Enterprises** | Postal code, main NACE code. |
| **Prodcom** | Manufacturing production: PC8, year, value, quantity, unit. |

**Sample construction:** Binarize transactions to 0/1. Drop firms with only one supplier. Drop transactions < 1% of customer purchases. Final sample: 7.3M transactions, 840k customers, 480k suppliers. Mean suppliers per customer: 8.7 (median: 7).
