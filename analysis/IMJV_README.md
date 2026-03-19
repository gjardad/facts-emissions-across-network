# IMJV to EUTL Matching Procedure

## What is IMJV?

IMJV (Integraal Milieujaarverslag) is Flanders' mandatory environmental reporting system. Firms in Flanders consuming more than 0.1 PJ of primary energy per year must report annual emissions. The data covers 178 unique firms with CO2 emissions over 2004–2023.

IMJV identifies firms by CBE/KBO numbers (Belgian enterprise identifiers) and firm names, not by the anonymized VAT codes used in NBB data.

## Why match IMJV to EUTL?

To assess the marginal value of IMJV — how many firms with real, reported emissions it adds beyond what EU ETS already covers. Unlike Climate TRACE (which assigns uniform top-down allocations to non-ETS facilities), IMJV provides actual facility-level reported emissions.

## Matching procedure

### Step 1: Algorithmic name-based matching

For each IMJV firm, we search all Belgian EUTL installations and compute:

**Name similarity** (same functions as CT-to-EUTL matching, see analysis/README.md):
- Bigram Jaccard index
- Token overlap score (with stopword removal)
- Combined: S_name = max(bigram Jaccard, token overlap)

**Substring check**: whether one firm name (stripped to lowercase alphanumeric) is contained within the other. Catches cases like TEEPAK inside "ViskoTeepak".

**City match**: whether the IMJV gemeente matches the EUTL city field (exact or substring).

**NACE match**: at 2-digit or 4-digit level.

Auto-accept tiers (in priority order):

| Tier | Name sim | City | NACE 2d | Label |
|------|----------|------|---------|-------|
| 1 | >= 0.4 | any | Yes | name_nace |
| 2 | >= 0.4 | any | any | name |
| 3 | >= 0.2 | Yes | Yes | name_city_nace |
| 4 | >= 0.2 | Yes | any | name_city |
| 5 | substring | Yes | Yes | substr_city_nace |

### Step 2: Flagging ambiguous cases

Unmatched IMJV firms that have an EUTL installation in the same city with the same NACE code (but zero name similarity) are flagged for manual review. These could be:
- Corporate name changes (the same physical facility under a new company name)
- Genuinely different firms that happen to share a city and sector

We deliberately do NOT auto-accept city+NACE matches because of false positive risk, especially in cities with multiple firms in the same sector (e.g., Antwerp chemicals).

### Step 3: Manual verification of flagged cases

Flagged cases were investigated via web search to determine whether the IMJV firm and EUTL installation are the same physical facility. Confirmed matches are hard-coded in the script.

## Hard-coded matches (verified corporate name changes)

The following IMJV-to-EUTL links were verified as corporate name changes, mergers, or acquisitions:

| IMJV name | EUTL name | Explanation | Source |
|-----------|-----------|-------------|--------|
| BP CHEMBEL | INEOS Aromatics Belgium | BP sold its global petrochemicals business to INEOS in 2020. Same plant in Geel. | [BP press release](https://www.bp.com/en/global/corporate/news-and-insights/press-releases/bp-agrees-to-sell-its-petrochemicals-business-to-ineos.html) |
| OUDEGEM PAPIER | VPK Paper | VPK Group was founded as a paper factory in Oudegem in 1936. Same site at Oude Baan 120, Dendermonde. | [VPK Group Wikipedia](https://en.wikipedia.org/wiki/VPK_Group) |
| CYTEC SURFACE SPECIALTIES | ALLNEX Belgium | Cytec acquired UCB's Surface Specialties in 2005. Solvay acquired Cytec in 2015. The coating resins business was spun off as Allnex. Same address: Anderlechtstraat 33, Drogenbos. | [Allnex announcement](https://allnex.com/en/info-hub/news/former-cytec-industries-coating-resins-business-be) |
| SOLVIN | INOVYN BELGIUM | Solvin was a 75/25 Solvay/BASF PVC joint venture. BASF sold its stake in 2015 and Solvay merged chlorovinyls with INEOS to form INOVYN. | [ICIS report](https://www.icis.com/explore/resources/news/2015/07/01/9899962/basf-sells-solvin-stake-to-solvay-as-inovyn-chlorvinyls-jv-starts-up/) |
| SOLVIC | BASF Antwerpen - 127c | Solvic was the VCM-PVC production unit at the BASF Antwerp Verbund site. Operations were shut down and absorbed into the BASF complex. | [BASF Antwerpen partners](https://www.basf.com/be/en/who-we-are/Group-Companies/BASF-Antwerpen/About-the-site/Partners-on-site) |
| LATEXCO | Novaya Belgium | Latexco (founded 1953, Tielt) merged with Artilat in 2023 to form Novaya. Same site. | [KennisWest](https://www.kenniswest.be/organisatie/latexco-novaya-nv-tielt/22938) |
| MISA ECO | Nesar | MISA ECO produced sulfuric acid in Sint-Kruis-Winkel (Gent). After bankruptcy, Nesar took over the site. Nesar also subsequently closed (2012). Same physical location. | [VRT NWS](https://www.vrt.be/vrtnws/nl/2012/05/04/53_banen_verlorenbijchemischbedrijfnesar-1-1291315/) |
| NORTH EUROPEAN SULFERIC ACID REGENERATION | Nesar | Same site as MISA ECO — predecessor company at the same location. | Same as above |
| ORRION CHEMICALS REGEN | Nesar | Intermediate operator between MISA ECO and Nesar at the same site. | Same as above |
| V.B.G. | VBG nv | Same company — "V.B.G." is the abbreviated form of "VBG nv". Same city (Wijnegem), same NACE (23.99). | User-confirmed |
| ASFALT PRODUCTIE LIMBURG | APL nv | Same company — "APL" is the acronym for "Asfalt Productie Limburg". Same city (Heusden-Zolder), same NACE (23.99). | User-confirmed |

**EDF LUMINUS / S.P.E. entries**: S.P.E. was renamed EDF Luminus in 2011 (after EDF acquired Centrica's stake in SPE-Luminus in 2009), then renamed Luminus in 2019. EDFL refers to Electrabel-Luminus joint operations. These entries are matched by CBE establishment suffix:

| CBB suffix | Location | EUTL installation | Source |
|------------|----------|-------------------|--------|
| 0137 | Harelbeke | EDFL - Centrale Harelbeke | [Luminus Wikipedia](https://en.wikipedia.org/wiki/Luminus_(company)) |
| 0238 | Gent | Electrabel - Langerbrugge | [EDF Luminus thermal](https://edfluminus.edf.com/en/edf-luminus/activities/produce-energy/thermal-power) |
| 0339 | Gent | Electrabel - Langerbrugge | Same source |
| 0743 | Izegem | EDFL - Izegem | Same source |

## Cases verified as NOT the same firm

The following flagged pairs were confirmed to be different firms:

| IMJV name | EUTL candidate | Reason |
|-----------|----------------|--------|
| COMINBEL | Inex | Both in Sint-Lievens-Houtem but different sectors: Cominbel = meat processing (NACE 10.41), Inex = dairy (NACE 10.51). |
| ADPO GHENT | Rousselot | ADPO is a chemical tank storage terminal. Rousselot is a gelatin manufacturer. Different companies. |
| HALTERMANN | 3M Belgium | Haltermann Carless operates blending sites in Ghent, not in Zwijndrecht where 3M is located. Different firms. |
| COGEBI | Terca Beerse | Different cities (Beersel vs Beerse) and different NACE (23.99 vs 23.32). |

## Output

`{PROC_DATA}/imjv_eutl_match.RData` contains:
- `imjv_eutl` — matching table with columns: cbb_number, imjv_name, imjv_gemeente, imjv_nace, mean_annual_co2_t, eutl_name, eutl_city, eutl_nace, name_sim, substr_match, city_match, nace_match_2d, nace_match_4d, matched, match_type
- `imjv_firms` — unique IMJV firms with summary statistics

Match types: `hardcoded`, `name_nace`, `name`, `name_city_nace`, `name_city`, `substr_city_nace`, `none`
