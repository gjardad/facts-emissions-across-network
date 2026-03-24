# HS6 to EXIOBASE v3 Industry Crosswalk

## Purpose

Maps HS6 product codes to EXIOBASE v3's 163 industries. This crosswalk
enables linking Belgian firms' Customs imports/exports (reported in CN8 codes)
to foreign (sector x country) nodes in an extended production network that
combines domestic B2B data with EXIOBASE's multi-regional input-output tables.

**Downstream use:** For a Belgian firm importing CN8 product _p_ from country
_c_, truncate _p_ to 6 digits (= HS6), look up the EXIOBASE industry via this
crosswalk, and connect the firm to the (industry, country) node in EXIOBASE.

## Mapping Chain

```
CN8  --(truncate first 6 digits)-->  HS6  --(WITS)-->  ISIC Rev 3 (4-digit)  --(hardcoded)-->  EXIOBASE industry
```

## Classification Versions

| Classification | Version | Notes |
|---|---|---|
| CN | 8-digit, yearly vintages | EU extension of HS; first 6 digits = HS6 |
| HS | 2002, 2007, 2012 (stacked) | Later vintage takes precedence |
| ISIC | Rev 3 | ~1:1 with NACE Rev 1.1 at 2-digit level |
| NACE | Rev 1.1 | Basis for EXIOBASE v3 industry classification |

## Source Data

| File | Source | Description |
|---|---|---|
| `hs2_isic3_wits.csv` | World Bank WITS | 5,224 HS 2002 to ISIC Rev 3 mappings |
| `hs3_isic3_wits.csv` | World Bank WITS | 5,053 HS 2007 to ISIC Rev 3 mappings |
| `hs4_isic3_wits.csv` | World Bank WITS | 5,177 HS 2012 to ISIC Rev 3 mappings |
| `industries.txt` | EXIOBASE v3 (Zenodo) | 163 industries with NACE Rev 1.1-based codes |

## HS Version Stacking

The three WITS concordances (H2, H3, H4) are stacked to maximise HS6
coverage. When an HS6 code appears in multiple vintages, the **latest vintage's
ISIC assignment takes precedence**. This handles two problems:

1. **New codes** introduced in HS 2007/2012 that don't exist in HS 2002
   (583 codes recovered).
2. **Redefined codes** where the product scope changed across revisions, causing
   the correct ISIC assignment to differ (120 codes corrected, of which 32
   changed EXIOBASE industry).

Examples of redefined codes:
- HS 080620 (dried grapes): H2 maps to ISIC 0113 (fruit growing) -> `i01.d`;
  H4 maps to ISIC 1513 (fruit processing) -> `i15.i`. H4 is correct.
- HS 040120 (milk, 1-6% fat): H2 maps to ISIC 1520 (dairy manufacturing) ->
  `i15.f`; H4 maps to ISIC 0121 (cattle farming) -> `i01.i`. H4 is correct.

## Output

**`hs6_exiobase_crosswalk.csv`** — 5,807 rows

| Column | Description |
|---|---|
| `hs6` | HS6 code (zero-padded string) |
| `isicrev3_4d` | ISIC Rev 3 4-digit code (zero-padded string) |
| `nace_div` | NACE Rev 1.1 division (first 2 digits of ISIC code) |
| `exio_code` | EXIOBASE industry CodeNr (e.g., `i23.2`, `i27.a`) |
| `exio_name` | EXIOBASE industry name |
| `exio_num` | EXIOBASE industry number (1-163, for matrix indexing) |
| `mapping_quality` | Confidence flag (see below) |

## Mapping Quality Flags

| Flag | Count | Meaning |
|---|---|---|
| `division` | 2,969 | ISIC code in a NACE division with a single EXIOBASE industry |
| `exact` | 1,303 | Unambiguous 1:1 mapping to a specific EXIOBASE sub-industry |
| `residual` | 1,524 | Multiple EXIOBASE sub-industries possible; assigned to the most likely one |
| `unmapped` | 11 | No meaningful EXIOBASE goods industry (all are ISIC 9999) |

## Coverage

- **5,796 of 5,807 HS6 codes mapped** (99.8%)
- **55 of 163 EXIOBASE industries reached** from trade data
- Unreached industries are services (NACE 50+), waste/recycling (i20.w,
  i24.a.w, etc.), and electricity sub-industries by source (i40.11.a-k)

## Design Decisions

### One-to-one mapping
Each HS6 code maps to exactly one EXIOBASE industry. Where ISIC collapses
multiple EXIOBASE sub-industries (e.g., ISIC 2720 covers all non-ferrous
metals but EXIOBASE splits by metal type), we assign to the most plausible
sub-industry and flag as `residual`. This keeps the downstream network
construction simple (no weighted edges).

### Output at HS6 level
The crosswalk is at HS6, not CN8. CN8 to HS6 is a trivial truncation (first 6
digits). Producing CN8-level output would replicate the HS6 mapping across
hundreds of CN8 sub-codes with no additional information.

### Waste/recycling industries excluded
EXIOBASE industries for secondary material processing (i20.w, i24.a.w, i26.a.w,
i27.a.w, etc.) are not reachable from trade data and are excluded.

### Electricity
HS code 271600 (electrical energy) maps to ISIC 4010, which EXIOBASE splits
into 12 sub-industries by generation source. Since trade data cannot distinguish
the source, we assign to `i40.11.l` (Production of electricity nec).

## Known Limitations

1. **HS 2017+ not covered.** H5 (HS 2017) to ISIC Rev 3 is not available on
   WITS. Codes introduced in HS 2017 or later may go unmatched. Could be
   constructed by chaining H5->H4->ISIC3 if needed.

2. **ISIC collapses HS-level information.** Some HS6 codes within the same ISIC
   code could be mapped to different EXIOBASE sub-industries (e.g., HS 270900
   = crude oil vs HS 271121 = natural gas both map to ISIC 1110). The ISIC
   bottleneck prevents this finer assignment.

3. **Residual mappings are approximate.** 1,524 HS6 codes (26%) are assigned to
   the dominant EXIOBASE sub-industry within their ISIC code, which may
   introduce noise for specific products.

4. **11 unmapped HS6 codes.** All are ISIC 9999 (miscellaneous goods not
   elsewhere classified). These represent a negligible share of trade value.

## Build Script

`preprocess/build_hs6_exiobase_crosswalk.R` — run from repo root.

## References

- WITS Product Concordance: https://wits.worldbank.org/product_concordance.html
- EXIOBASE v3: https://zenodo.org/records/5589597
- Stadler et al. (2018), "EXIOBASE 3", Journal of Industrial Ecology
- Magerman (2022), "Correspondences of EU Product Classifications"
