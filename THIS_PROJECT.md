## Project Overview: The aim of the project is establish facts about the distribution of emissions across the producton network in Belgium between 2005-2021.

### Research Questions:

1. How dispersed are  scope 1 emissions across firms within sectors?
2. How disperse are upstream network-adjusted emissions across firms within sectors? How does this dispersion compare with within-sector, across firms dispersion in scope 1 emissions?
3. Does upstream network-adjusted emissions correlate with scope 1 emissions at the firm level across firms within sectors?
4. To what extent do firm-level organizational boundaries — rather than underlying production technology — drive within-sector heterogeneity in scope 1 emissions? (Measured by within-sector rank correlation between scope 1 and network-adjusted emissions; low correlation indicates outsourcing decisions distort scope 1 rankings.)
5. How concentrated is each firm's emission exposure across its suppliers?
6. How many steps upstream does most embodied emission come from?

#### Note on RQ 4

Comparing within-sector rankings under scope 1 vs. network-adjusted emissions reveals the extent to which the **organizational boundary** (vertical integration vs. outsourcing) distorts emission rankings. If rank correlation is well below 1 in a sector, it means that make-or-buy decisions — not underlying production technology — are driving scope 1 heterogeneity.

**When it matters and when it doesn't**

The comparison is **uninformative** in sectors where the dirty node is internal to the sector regardless of firm structure. Example: cement vs. ready-mix concrete (both NACE 23). Cement producers operate kilns (high scope 1); concrete producers buy cement (low scope 1). But the concrete producers' network-adjusted emissions trace back to the same cement kilns — so both scope 1 and network-adjusted measures identify the cement producers as the high emitters. Network adjustment doesn't change the ranking of who matters.

The comparison is **informative** in sectors where the dirty node is **outside** the sector — i.e., firms purchase carbon-intensive inputs from suppliers in other NACE sectors. In these sectors, scope 1 is uniformly low but firms vary in how much they buy from high-emission upstream suppliers. Network-adjusted emissions reveal within-sector heterogeneity that scope 1 cannot.

**Practical examples**

- **Plastics manufacturing (NACE 22):** Own processes (molding, extrusion) are low-temperature and low-emission. The key input is polymer resin from petrochemical producers (NACE 20), where steam cracking of naphtha is extremely energy-intensive. Within NACE 22, firms vary in how much virgin polymer they purchase from Antwerp crackers vs. using recycled plastics. Scope 1 looks similar; network-adjusted emissions diverge.

- **Fabricated metals (NACE 25):** Firms buy steel and aluminum from primary metals producers (NACE 24) and perform cutting, welding, coating — modest energy use. Variation in how much primary metal a firm purchases from blast-furnace steelmakers shows up only in network-adjusted emissions.

- **Paper products (NACE 17):** Non-integrated firms buy pulp from energy-intensive pulp mills. Their scope 1 is low but their upstream exposure is high. Integrated firms (with on-site pulping) have high scope 1 but similar network-adjusted emissions.

#### Note on RQ 5

For each firm, compute a Herfindahl-type index over its suppliers weighted by supplier emission intensity. Do firms spread their sourcing across many low-emission suppliers, or concentrate purchases with a few high-emission ones? This has direct implications for transition risk — firms with concentrated dirty-supplier exposure are more vulnerable to carbon pricing pass-through. The B2B data provides the transaction-level resolution needed to compute this; aggregate IO tables cannot.

### Emissions Data

Firm-level emissions are constructed from four sources:

**1. EU ETS.** The EU ETS compliance data provides verified annual emissions for each regulated installation in Belgium, covering more than 70% of emissions from stationary sources. Installations are linked to firms via the EUTL account data, which maps installation IDs to firm identifiers (BvD IDs and anonymized VAT codes). Firm-level EU ETS emissions are obtained by summing verified emissions across all installations owned by a given firm in a given year.

**2. Climate TRACE — not used.** We investigated Climate TRACE (CT) as a potential source of facility-level emissions for non-ETS firms. CT covers 367 sources in Belgium (2021–2025) across 6 sectors. After matching CT sources to EUTL installations using geographic proximity, name similarity, and sector compatibility (see analysis/README.md for the full procedure), we found that CT does not add useful information for this project:

- *High overlap with EUTL.* Of the facility-level CT sources (manufacturing, power, fossil fuel operations), most match to existing EUTL installations. Manufacturing subsectors like cement, glass, petrochemicals, and refineries have near-100% match rates.
- *Non-overlapping sources lack facility-level variation.* The ~25 genuinely non-ETS CT sources (mostly food processing and textiles) are estimated using CT's "Data-informed Emissions Disaggregation" method, which divides country-level emissions from EDGAR/CEDS equally across all known facilities in a subsector. All 17 food-beverage-tobacco sources in Belgium receive identical emissions (60,454 t/yr each); all 5 textile sources receive 15,894 t/yr each. There is zero cross-facility variation.
- *Noisy where it overlaps.* On the 226 matched installation-years where both CT and EUTL have positive emissions, CT systematically overestimates (median CT/EUTL ratio: 1.25) with moderate rank correlation (Spearman rho: 0.66).

The imputation model from inferring_emissions is strictly preferable for non-ETS firms since it exploits firm-level observables (B2B transactions, revenue, sector) rather than dividing a country total by N.

**3. IMJV.** [TODO: describe what this data provides.]

**4. Imputed emissions.** For firms not covered by the sources above, emissions are imputed using the prediction model developed in the inferring_emissions project.

### NACE source convention for the CRF mapping

The `nace2d → crf_group` lookup is the same in all scripts, but the source of `nace2d` differs by firm type:

- **ETS firms (true ETS + pre-ETS backcast):** `nace5d` field on `firm_year_belgian_euets.RData`. This was sourced from annual accounts upstream (it is *not* `nace_id_from_eutl`, which is the installation-level NACE). Coverage is structurally ~94% of distinct ETS vats — annual accounts simply don't cover every ETS firm (unfiled accounts, foreign-domiciled firms with Belgian installations, etc.). ETS firms with NA `nace2d` are dropped from cell ranking and from the `E_ETS` aggregate, with a logged warning.
- **Non-ETS firms (deployment):** `nace5d` field on `deployment_panel.RData` (also annual-accounts-derived).

Crucially, the GLO allocation pipeline (`analysis/b_allocation_glo.R` and the upstream `analysis/fit_glo_reference.R`, `analysis/calibrate_pooled_youden_tau.R`) does **not** depend on `annual_accounts_selected_sample_key_variables.RData` for the NACE join. That file is downsampled on local 1; using the EUTL panel's own `nace5d` for ETS firms makes the allocation runnable end-to-end on local 1 with identical filenames and identical conventions to RMD.

This convention matches the `inferring_emissions/` project, which also uses annual-accounts NACE (not EUTL installation NACE) when assigning ETS firms to CRF groups.

### Sector Conventions

**NACE 17 (paper) and NACE 18 (printing) are ALWAYS treated as a single sector "17/18". NEVER treat them as separate sectors.** Both sectors share upstream pulp/paper supply chains and the distinction is uninformative for emission analysis. This applies everywhere: stratification, fixed effects, dispersion statistics, Youden threshold calibration, and any sector-level aggregation.

### Data Sources

All data at the firm-level contain an unique firm-level identifier which is an anonymized VAT code. The code makes it possible to merge data sets.

NBB data sets (located in DATA_DIR/raw/NBB):

	Customs data: firm-product-year-level data on imports and exports for the universe of goods imported into and exported from Belgium between 2000 and 2022. Goods are 	identified by CN 8-digit product codes and firms are identified by anonymized VAT codes.

	B2B: buyer-supplier-year-level data on transactions between any two given private firms in Belgium between 2002 and 2022. It only contains the total nominal euros any   	two given firms transact, not which products were transacted. Buyers and suppliers are identified by their anonymized VAT code. The raw data is B2B\_ANO.dta

	Annual accounts: a large set of firm-year-level characteristics. It includes, among others, revenue, wage bill, capital, and nace 5-digit sectors. The raw data is 	called Annual\_Accounts\_MASTER\_ANO.dta

	PRODCOM: firm-product-year data on quantities and prices for each 8-digit code product produced by any given firm that is part of the survey sample. It only covers 	manufacturing goods (NACE codes between 07 and 33), for a sample of around 200k firms.

	EUTL\_Belgium: for each installation covered by the EUTL, it provides unique firm identifiers (BvD id and the corresponding anonymized VAT).

EUTL data sets (located in DATA_DIR/raw/EUTL/Oct_2024):

	account: for each installation, it informs account unique id and corresponding BvD id.

	compliance: for each installation-year regulated by the EUETS, it contains amount of emission permits allocated to the installation and verified emissions on that year. 	Installation are identified by unique installation\_id.

	installation: for each installation regulated by the EUETS, it contains installation\_id, geographic location (lat/lon), NACE sector id (2 or 4 digit), among others.

NIR data (located in DATA_DIR/raw/NIR):

	BEL-CRTs: for each year between 1990 and 2023, it informs GHG emissions by CRF category. The data is split into multiple .xlsx, one for each year.

	Annex XII: data on total GHG emissions and share of it regulated by the EUETS by CRF category. It consists of two .xslx files, one for 2022 (published in 2024) and one 	for 2023 (published in 2025).

CRF-to-NACE crosswalk (located in DATA_DIR/raw/Correspondences_and_dictionaries): maps CRF categories to NACE Rev. 2 sectors. Combined with BEL-CRTs, this enables sector-level calibration targets for deployment to non-EUETS firms.

Processed data (located in DATA_DIR/processed):

	training_sample.RData: firm-year panel used for cross-validation. Contains EU ETS firms (with verified emissions) and non-ETS firms from NACE 17/18, 19, and 24 (with emissions set to 0). Available on local 1.

**Downsampled data on local 1.** The following processed files on local 1 are **downsampled** versions of the full data (which is available only on RMD). The same applies to the corresponding raw .dta files in DATA_DIR/raw/NBB/.

- `annual_accounts_selected_sample.RData` (and `_key_variables` and `_more_selected_sample` variants)
- `b2b_selected_sample.RData`
- `df_b2b.RData`
- `df_trade.RData`
- `df_national_accounts_with_5digits.RData`
- `firms_in_selected_sample.RData`

The full training sample (`training_sample.RData``) are NOT downsampled and are available in full on local 1.

### Hardware setup:

**I use three desktops in this project: local 1, local 2, and RMD (remote desktop).
Access to the full NBB data is restricted to RMD through a VPN connection. RMD doesn't have access to the web browser, but it is connected to GitHub. I can only use the VPN connection through local 2. Local 2 has regular access to the web browser. Local 1 is my personal desktop and it is where I have Claude code and cursor downloaded.
When copying files from RMD to local 1, I first need to copy them to local 2, then from local 2 to the cloud (Dropbox/Claude), then from the Claude to local 1.
In local 1 I have available a downsampled version of the full NBB data sets as well as the full training sample. I built the training sample in RMD and copied it to local 2.
Any script that only requires `training_sample.RData` (e.g., CV scripts, alternative specs, rho comparisons) can be run locally on local 1. RMD is only needed for scripts that access the raw NBB data (e.g., preprocessing, proxy construction).

Files transferred from RMD to local 1 (full-data, not downsampled):
- `repeated_cv_proxy_crf_asinh.RData` — 200 CV repeats for training firms (M=200 × N=26,608)
- `deployment_proxy_list.RData` — 200 EN repeats for deployment firms (M=200, ~345k firm-years per draw). Generated on RMD with full B2B data, transferred 2026-04-08. All downstream outputs on local 1 that derive from this file (uncertainty draws, GLO allocation, scope 1 dispersion) use the full-data proxies and are therefore production-quality, not downsampled.**

### Preliminary Results: Scope 1 Carbon Productivity Dispersion (RQ1)

We compute within-sector dispersion of scope 1 carbon productivity (revenue / scope 1 emissions) at three levels of aggregation: NACE 5-digit, NACE 2-digit, and CRF category. Point estimates come from the deterministic GLO allocation; prediction intervals from 200 perturbation draws. All inputs are full-data (deployment_panel.RData for revenue/NACE, allocation_glo_balanced for scope 1, uncertainty_draws_balanced for PIs). Results cover 2005-2021, all sectors (not restricted to manufacturing), all firm sizes.

**Carbon productivity p90/p10 (median across sector-years, point estimate):**
- NACE 5-digit: 23.0x (90-10 log diff: 3.24)
- NACE 2-digit: 72.7x (90-10 log diff: 4.32)
- CRF category: 147.0x (90-10 log diff: 5.03)

**Comparison with the literature.** The most comparable benchmark is De Lyon & Dechezlepretre (2025, OECD WP 2025/24), who report carbon productivity p90/p10 of 24x averaged across 9 countries, within 3-digit NACE manufacturing industries (20+ employees, 10+ firms per cell). This is an unweighted average of the log 90-10 difference across industry-country-year cells, exponentiated. Their energy productivity p90/p10 by country (Figure 3.1) ranges from 12x (Croatia, Sweden) to 39x (Indonesia), with the Netherlands at 18x and France at 14x. For France specifically (Table 3.1), they show energy productivity p90/p10 declining from 17.5x at 2-digit to 12.2x at 4-digit — the same pattern of finer industries yielding less dispersion that we observe. Our NACE 5-digit result of 23x aligns closely with their 24x cross-country carbon productivity aggregate, despite differences in scope (we include all sectors and firm sizes, they restrict to manufacturing with 20+ employees) and in the emission measure (we include process emissions in scope 1, they measure energy-related CO2 only).

Lyubich, Shapiro & Walker (2018) report a 90-10 log difference of 2.27 (p90/p10 of 9.7x) for CO2 productivity in US manufacturing at the 6-digit NAICS plant level. The lower dispersion is consistent with their much finer industry definition and plant-level (vs. firm-level) unit of analysis.

### Dropped year: 2022

Year 2022 is excluded from the analysis. The annual accounts data for 2022 contains placeholder wage_bill values (EUR 1–2) for a large number of firms that had not yet filed their accounts at the time the data was extracted. While these firms appear in the data with 100% coverage, the near-zero wage bills cause the cost-based denominator of the Leontief A matrix to be dominated by domestic B2B purchases, pushing row sums to 1.0 for 51% of firms (vs. 1% in 2021). This invalidates the Leontief inverse computation for 2022. The sample covers 2005–2021 (17 years). See `A_ITERATION.md` for full diagnostics.
