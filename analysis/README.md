# Climate TRACE to EUTL Matching Procedure

## Problem

Climate TRACE (CT) and the EU Transaction Log (EUTL) both track emissions from stationary installations in Belgium, but use different identifiers, naming conventions, and coordinate systems. Matching CT sources to EUTL installations is necessary to (a) validate CT emission estimates against EUTL verified emissions, and (b) identify CT sources that fall outside EU ETS coverage.

Matching is non-trivial because:
- Names differ in format (e.g., "Lhoist Industrie Sa - Site De On" vs "Usine de On")
- Coordinates may refer to different points within the same industrial site
- Dense industrial zones (e.g., Antwerp port) contain multiple facilities within 1 km
- CT subsectors and EUTL activity types use different classification systems

## Procedure overview

The matching proceeds in six steps:

1. **Candidate generation** — geographic filter
2. **Name similarity** — two complementary metrics
3. **Sector compatibility** — cross-classification filter
4. **Tier assignment** — four-tier decision rule
5. **One-to-one assignment** — greedy deduplication
6. **Diagnostics** — quality checks and ambiguity flags

## Step 1: Geographic distance

For a CT source at coordinates (phi_1, lambda_1) and an EUTL installation at (phi_2, lambda_2), the Haversine distance in km is:

    d = 2R * arcsin(sqrt(sin^2((phi_2 - phi_1)/2) + cos(phi_1) * cos(phi_2) * sin^2((lambda_2 - lambda_1)/2)))

where R = 6371 km (Earth's mean radius), and all angles are in radians.

All CT-EUTL pairs with d <= 10 km are retained as candidates. Pairs beyond 10 km are discarded.

## Step 2: Name similarity

Two complementary metrics are computed for each candidate pair. The combined similarity is the maximum of the two.

### Metric A: Bigram Jaccard index

Given two facility names a and b:

1. Strip to lowercase alphanumeric characters
2. Extract the set of consecutive 2-character substrings (bigrams):
   B(s) = {s[k..k+1] : k = 1, ..., len(s)-1}
3. Compute the Jaccard index:

       S_bigram(a, b) = |B(a) ∩ B(b)| / |B(a) ∪ B(b)|

This metric is robust to small insertions and character-level noise but fails when semantically equivalent names share few character sequences (e.g., "Usine de On" vs "Lhoist Industrie Sa - Site De On").

### Metric B: Token overlap score

1. Tokenize both names by splitting on whitespace and punctuation
2. Lowercase all tokens
3. Remove stopwords: {nv, sa, bv, bvba, sprl, srl, cvba, site, vestiging, usine, plant, de, le, la, les, het, van, du, des, den, der, een, the, and, et, en}
4. Keep only tokens with length >= 3
5. Compute:

       S_token(a, b) = |T(a) ∩ T(b)| / min(|T(a)|, |T(b)|)

   where T(s) is the set of unique non-stopword tokens of s.

The denominator uses the minimum of the two set sizes to handle asymmetric name lengths (e.g., a short EUTL name that is a substring of a longer CT name).

### Combined similarity

    S_name(a, b) = max(S_bigram(a, b), S_token(a, b))

## Step 3: Sector compatibility

A binary compatibility function C(i, j) = 1 if the CT subsector of source i is compatible with the EUTL activity type of installation j. The mapping is:

| CT subsector                 | Compatible EUTL activity_id(s)           |
|------------------------------|------------------------------------------|
| cement                       | 29 (cement clinker)                      |
| iron-and-steel               | 24 (pig iron/steel), 25 (ferrous metals) |
| glass                        | 31 (glass)                               |
| oil-and-gas-refining         | 21 (mineral oil refining)                |
| petrochemical-steam-cracking | 42 (bulk chemicals), 43 (hydrogen)       |
| pulp-and-paper               | 35 (pulp), 36 (paper/cardboard)          |
| lime                         | 30 (lime)                                |
| chemicals                    | 38, 39, 41, 42, 43                       |
| other-metals                 | 26, 27, 28 (aluminium, non-ferrous)      |
| electricity-generation       | 1, 20, 21, 42, 43 (see note below)      |
| food-beverage-tobacco        | 1 (combustion > 20 MW), 20 (combustion)  |
| textiles-leather-apparel     | 1 (combustion > 20 MW), 20 (combustion)  |

CT subsectors not in the table default to {1, 20} (generic combustion).

**Note on electricity-generation compatibility.** CT classifies combined heat and power (CHP) units as "electricity-generation" even when they are physically located within, and operationally part of, a refinery or chemical complex. In the EUTL, these CHP units are often registered under the co-located facility's entry (e.g., a refinery with activity_id 21, or a chemicals plant with activity_id 42). This is particularly common in the Antwerp port industrial zone, where multiple CHP plants serve petrochemical and refining operations. To avoid false negatives in these dense industrial areas, the electricity-generation subsector is also compatible with refining (21), bulk chemicals (42), and hydrogen/synthesis gas (43) activity types.

## Step 4: Tier assignment

Each candidate pair (i, j) is assigned to the highest applicable tier:

| Tier | Distance       | Name similarity        | Sector compatible | Interpretation                        |
|------|----------------|------------------------|--------------------|---------------------------------------|
| 1    | d <= 1 km      | S_name >= 0.2          | Yes                | Strong: close, name confirms, sector matches |
| 2    | d <= 1 km      | (any)                  | Yes                | Geographic + sector: very close with right activity type |
| 3    | 1 < d <= 5 km  | S_name >= 0.2          | Yes                | Medium range: name and sector both confirm |
| 4    | 5 < d <= 10 km | S_name >= 0.4          | (any)              | Extended: only if name strongly matches |

Pairs that do not qualify for any tier are discarded.

Tier 2 waives the name check for d <= 1 km when sector compatibility holds. This accommodates cases where the same facility has completely different names in the two databases (e.g., a company name in EUTL vs a facility description in CT) but the geographic and sectoral evidence is strong.

Tier 4 allows matches beyond 5 km when the name match is strong (>= 0.4), regardless of sector. This handles cases where EUTL coordinates are imprecise for large industrial complexes.

## Step 5: One-to-one greedy assignment

To prevent multiple CT sources from matching to the same EUTL installation (or vice versa):

1. Sort all qualifying candidate pairs by tier (ascending), then by distance (ascending) within tier
2. Iterate through sorted pairs:
   - If neither the CT source nor the EUTL installation has been claimed, accept the match and mark both as claimed
   - Otherwise, skip the pair

This greedy algorithm prioritizes higher-confidence matches (lower tier number, shorter distance).

## Step 6: Diagnostics

The script reports:
- Match counts by tier
- Distance and name similarity distributions
- Matched and unmatched CT sources by sector/subsector
- **Ambiguous cases**: CT sources with multiple EUTL candidates within 2 km (flagged for manual review)

## Output

`{PROC_DATA}/ct_eutl_match.RData` contains:
- `ct` — CT source-year panel (CO2 emissions, aggregated to annual)
- `ct_eutl_match` — matching table with columns: source_id, ct_name, ct_subsector, installation_id, eutl_name, eutl_activity, dist_km, sim_bigram, sim_token, sim_name, sector_compat, tier, match_type
- `be_inst` — EUTL Belgian installations (non-aircraft, non-maritime, with coordinates)
