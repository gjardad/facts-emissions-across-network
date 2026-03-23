# NACE–CRF Group Crosswalk

## Purpose

This crosswalk maps NACE Rev. 2 2-digit sectors to **CRF emission groups** so that
NIR sector-year totals (`E_NIR`) can be compared against ETS sector-year totals
(`E_ETS`) to compute the residual deployment emissions:

```
E_deploy[group, year] = max(0, E_NIR[group, year] - E_ETS[group, year])
```

`E_deploy` is the NIR budget distributed among non-ETS deployment firms within the
group, proportionally to their imputed emission proxy (see `analysis/compute_b_loop.R`).

## Files

| File | Description |
|---|---|
| `nace_crf_crosswalk.csv` | NACE 2-digit → CRF group (one row per NACE sector) |
| `crf_group_definitions.csv` | CRF group → regex pattern matching CRF codes in BEL-CRT data |

The regex patterns in `crf_group_definitions.csv` are matched against the compact
CRF codes in `belcrt_sector_denominators.tsv` (e.g. `1A2a`, `2C1`). All CRF
subcodes matching the pattern are summed to obtain `E_NIR` for the group.

## Source

Annex I to the AEA Manual (2015 edition):
*Correspondence between CRF-NFR and NACE Rev. 2 to the Manual for Air Emissions
Accounts*. Columns used: CRF code (col 9), CRF label (col 10), NACE Rev. 2 (col 15),
notes (col 16).

For CHP treatment: NIR Belgium sub2025 (submission 2025), pages 94–95 and 141–142.

## Scope

**Included CRF main categories:**
- **1A1** Energy industries (combustion)
- **1A2** Manufacturing industries and construction (combustion)
- **1A4** Other sectors, excluding transport (combustion)
- **2A–2C, 2E** Industrial processes: mineral industry, chemical industry, metal
  industry, electronics industry

**Excluded:**
- 1A3 Transport — see §Exclusions
- 1A4b Residential — see §Exclusions
- 1A5 Other stationary and mobile — too small, country-specific allocation
- 2D Non-energy products and solvents — see §Exclusions
- 2F ODS substitutes, 2G Other product manufacture, 2H Other — all reported as
  country-specific (CS) in the Annex I with no reliable single-NACE assignment

## CRF Groups

### `energy` — NACE 35
CRF 1A1a: Public electricity and heat production, including main-activity CHP.

NACE 35 has deployment firms: Annex XII to the NIR reports ~80% ETS coverage of
1A1a. The remaining 20% comprises small CHP plants and district heating operators
below the 20 MW ETS registration threshold.

### `refining` — NACE 19
CRF 1A1b: Petroleum refining.

Included for accounting completeness and as a verification check. In practice
E_deploy[refining] ≈ 0 because (i) ETS covers ~100% of C19 combustion and
(ii) non-ETS firms in NACE 19 are in the training sample with assumed zero
emissions, so there are no NACE 19 deployment firms to receive imputed emissions.

### `metals` — NACE 24, 25
CRF 1A1c + 1A2a + 1A2b + 2C (all subcodes).

**Assumption A — 1A1c absorbed into `metals`:** CRF 1A1c (Manufacture of solid
fuels and other energy industries) spans NACE B, C19, C24, and D — sectors assigned
to four different groups — making a clean single-group assignment impossible. The
dominant Belgian activity under 1A1c is coke-oven production (large ETS installations
in C24). We absorb all of 1A1c into `metals`. The fraction attributable to B (mining)
or D (energy utilities) is small in Belgium.

**Assumption B — C25 assigned to `metals`:** The Annex I assigns both 1A2a (iron &
steel) and 1A2b (non-ferrous metals) to "C24, C25 — split country-specific". Assigning
C25 to `mfg_other` would create a systematic mismatch: E_NIR[metals] would include
C25 combustion (reported jointly with C24 in 1A2a/b) while E_ETS[metals] would miss
any ETS firms in C25. Assigning C25 to `metals` aligns both sides of the residual.
Process emissions (2C) map cleanly to C24.

### `chemicals` — NACE 20, 21
CRF 1A2c + 2B (all subcodes).

1A2c (Chemicals combustion) maps to C20, C21 in the Annex I. All chemical-industry
process emissions (2B: ammonia, nitric acid, adipic acid, caprolactam, carbide,
titanium dioxide, soda ash, petrochemicals, fluorochemicals, other) map to C20.
NACE 21 (pharmaceuticals) has negligible process emissions under 2B; its combustion
is captured by 1A2c.

### `paper` — NACE 17, 18
CRF 1A2d: Pulp, paper and print combustion only.

Process emissions from printing (2D.3h, solvent use in C18) are excluded: they are
small, reported as country-specific, and not from fossil fuel combustion.

### `food` — NACE 10, 11, 12
CRF 1A2e: Food processing, beverages and tobacco combustion only.

No significant process emissions in category 2 for these sectors.

### `minerals` — NACE 23
CRF 1A2f + 2A (all subcodes).

1A2f (Non-metallic minerals combustion) maps cleanly to C23. Category 2A covers all
mineral-industry process emissions: cement (2A1), lime (2A2), glass (2A3), ceramics
and other carbonate use (2A4) — all assigned to C23 in the Annex I. This is the
group with the largest process-to-combustion ratio: cement calcination CO₂ dominates
total C23 emissions.

### `mfg_other` — NACE 05–09, 13–16, 22, 26–33, 41–43
CRF 1A2g + 2E (all subcodes).

CRF 1A2g (Other manufacturing and construction) is a catch-all pool covering mining
(1A2g.iii → B), wood (1A2g.iv → C16), construction (1A2g.v → F), textiles
(1A2g.vi → C13–C15), machinery (1A2g.i → C28), transport equipment (1A2g.ii →
C29–C30), and other manufacturing (1A2g.viii → C22, C26, C27, C31–C33). All
subcodes are summed into a single NIR pool.

Electronics-industry process emissions (2E → C27) are included since C27 is part of
this group. Category 2E is typically small in Belgium.

Off-road mobile combustion within manufacturing (1A2g.vii, country-specific) is
included if reported.

Note: C25 (fabricated metal products) is assigned to `metals`, not `mfg_other` —
see Assumption B above.

### `commercial` — NACE 36–39, 45–47, 52–53, 55–56, 58–96
CRF 1A4a (all subcodes: stationary combustion 1A4a.i and off-road mobile 1A4a.ii).

1A4a covers commercial and institutional combustion across nearly all service-sector
NACE sectors. E_ETS[commercial] is near zero (no large combustion plants in retail,
finance, education, etc.), so E_deploy[commercial] ≈ E_NIR[1A4a].

**CHP note:** Autoproducer CHP plants in service-sector firms (e.g. hospital
cogeneration) are allocated to 1A4a in Belgium's NIR and to their firm's NACE sector
in the EUTL. Both sides of the E_NIR − E_ETS comparison therefore land in
`commercial`, preserving coherence (NIR Belgium sub2025, pp. 94–95, 141–142).

NACE 49 (land transport), 50 (water transport), and 51 (air transport) are excluded
from this group — see §Exclusions.

### `agriculture` — NACE 01, 02, 03
CRF 1A4c (all subcodes: stationary 1A4c.i, off-road mobile 1A4c.ii, fishing
1A4c.iii).

Both stationary and off-road mobile combustion are included: an agricultural firm
burning diesel in tractors (1A4c.ii) has genuine combustion emissions attributable
to the firm-year.

## Exclusions

### Transport (NACE 49, 50, 51)
CRF 1A3 covers road, rail, aviation, and navigation combustion. These emissions are
not attributable to individual firm-years in our framework: (i) road transport
emissions are distributed across millions of vehicle owners including households;
(ii) the EN proxy is trained on stationary industrial combustion, making imputed
values for transport firms uninformative. NACE 49–51 have no CRF group assignment
and receive no imputed emissions.

### Residential (1A4b)
Household heating and personal vehicles have no firm counterpart. NACE 97–98
(households as employers of domestic workers) are excluded.

### Non-energy products and solvents (2D)
Category 2D covers emissions from using hydrocarbons as raw materials or solvents
(asphalt paving → F, printing solvents → C18, lubricants, etc.). Excluded because:
(i) most 2D subcategories are country-specific with no reliable NACE assignment;
(ii) they are not fossil fuel combustion emissions, which is the mechanism our model
targets; (iii) they are small relative to combustion in all affected sectors.

## Known Limitations

1. **`mfg_other` pool is large and heterogeneous.** Sectors ranging from coal mining
   (05) to furniture manufacturing (31) and construction (41–43) share a single NIR
   pool. The proxy signal helps distribute within the pool, but sectors with very
   different combustion intensities compete for the same E_deploy budget.

2. **1A1c absorption introduces a small cross-sector bias.** If Belgian coke-oven
   combustion is a minority of 1A1c (i.e. most of it is actually in B or D), then
   E_NIR[metals] is overstated and E_deploy[metals] inflated. Given Belgium's
   industrial structure (Liège and Charleroi steel complexes historically dominant
   in coke production) this fraction is expected to be small.

3. **E_deploy floored at zero.** When ETS verified emissions exceed the NIR total
   within a group-year, E_deploy is set to zero. This can result from cross-walking
   approximations or year-specific anomalies in ETS reporting. The floor prevents
   negative imputed emissions but means some NIR budget goes undistributed in
   affected group-years. A warning count is logged in `analysis/compute_b_loop.R`.

4. **NIR–ETS reconciliation within group-years is approximate.** Belgium guarantees
   consistency between CRT totals and ETS data at the national level via supplementary
   data sources (integrated environmental reports). Residual discrepancies at the
   sector-year level — from CHP reclassification, installation-level reporting lags,
   or cross-sector spill-overs — are absorbed by the floor.
