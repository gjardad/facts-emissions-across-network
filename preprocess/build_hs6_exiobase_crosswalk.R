###############################################################################
# preprocess/build_hs6_exiobase_crosswalk.R
#
# PURPOSE
#   Build a crosswalk from HS6 product codes to EXIOBASE v3 industries
#   (163 industries, based on NACE Rev 1.1).
#
#   Mapping chain:
#     CN8 -> HS6 (truncate first 6 digits) -> ISIC Rev 3 (WITS) -> EXIOBASE
#
#   The ISIC Rev 3 -> EXIOBASE step is hardcoded domain knowledge.
#   See preprocess/crosswalks/CROSSWALK_EXIOBASE.md for full rationale.
#
# INPUT
#   {CROSS_DIR}/hs2_isic3_wits.csv  (HS 2002 -> ISIC Rev 3, 5224 rows)
#   {CROSS_DIR}/hs3_isic3_wits.csv  (HS 2007 -> ISIC Rev 3, 5053 rows)
#   {CROSS_DIR}/hs4_isic3_wits.csv  (HS 2012 -> ISIC Rev 3, 5177 rows)
#       — Source: World Bank WITS (H2/H3/H4 to I3 concordances)
#   {RAW_DATA}/Exiobase/IOT_2005_ixi/industries.txt
#       — 163 EXIOBASE v3 industries with codes and names
#
# OUTPUT
#   {REPO_DIR}/preprocess/crosswalks/hs6_exiobase_crosswalk.csv
#
# CLASSIFICATION VERSIONS
#   HS 2002, 2007, 2012 (stacked; later vintage takes precedence)
#   ISIC Rev 3 (~= NACE Rev 1.1, the basis for EXIOBASE v3)
#
# STACKING LOGIC
#   We stack H2, H3, H4 concordances to maximise HS6 coverage:
#   - Codes only in H2: use H2's ISIC assignment
#   - Codes in H2 and H3/H4: use latest vintage (handles redefined codes)
#   - Codes only in H3/H4: picked up as new codes
#   This recovers ~583 codes not in H2 and corrects ~32 EXIOBASE assignments
#   for codes whose product scope changed across HS revisions.
#
# KNOWN LIMITATION
#   HS 2017 (H5) -> ISIC Rev 3 concordance not available on WITS. Codes
#   introduced in HS 2017 or later may still go unmatched. H5 to ISIC Rev 3
#   could be constructed by chaining H5->H4->ISIC3 if needed.
#
# RUNS ON: local 1 or RMD
###############################################################################

source("paths.R")

CROSS_DIR <- file.path(REPO_DIR, "preprocess", "crosswalks")
EXIO_DIR  <- file.path(RAW_DATA, "Exiobase", "IOT_2005_ixi")

cat("==============================================================\n")
cat("  BUILD HS6 -> EXIOBASE CROSSWALK\n")
cat("==============================================================\n\n")

# =============================================================================
# Section 1: Read and stack WITS concordances (H2, H3, H4)
# =============================================================================

read_wits <- function(file, label) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  names(df) <- c("hs6", "hs_desc", "isicrev3_4d", "isic_desc")
  df$hs6         <- sprintf("%06d", as.integer(df$hs6))
  df$isicrev3_4d <- sprintf("%04d", as.integer(df$isicrev3_4d))
  cat(sprintf("   %s: %d rows, %d unique HS6\n", label, nrow(df), length(unique(df$hs6))))
  df[, c("hs6", "isicrev3_4d")]
}

cat("-- Reading WITS concordances ...\n")
h2 <- read_wits(file.path(CROSS_DIR, "hs2_isic3_wits.csv"), "H2 (HS 2002)")
h3 <- read_wits(file.path(CROSS_DIR, "hs3_isic3_wits.csv"), "H3 (HS 2007)")
h4 <- read_wits(file.path(CROSS_DIR, "hs4_isic3_wits.csv"), "H4 (HS 2012)")

# Stack: later vintage takes precedence for codes that appear in multiple
# vintages (handles redefined codes). For codes only in H2, H2's assignment
# is kept. For codes in H3/H4 not in H2, they are added.
h2$vintage <- 2L
h3$vintage <- 3L
h4$vintage <- 4L
stacked <- rbind(h2, h3, h4)
# Keep highest vintage per HS6 code
stacked <- stacked[order(stacked$hs6, -stacked$vintage), ]
wits <- stacked[!duplicated(stacked$hs6), c("hs6", "isicrev3_4d")]

cat(sprintf("\n   Stacked: %d unique HS6 codes, %d unique ISIC Rev 3 codes\n",
            length(unique(wits$hs6)), length(unique(wits$isicrev3_4d))))
cat(sprintf("   (H2-only: %d, added from H3: %d, added from H4: %d, redefined: %d)\n\n",
            sum(!wits$hs6 %in% c(h3$hs6, h4$hs6)),
            sum(wits$hs6 %in% setdiff(h3$hs6, h2$hs6)),
            sum(wits$hs6 %in% setdiff(h4$hs6, union(h2$hs6, h3$hs6))),
            sum(wits$hs6 %in% intersect(h2$hs6, h4$hs6) &
                wits$isicrev3_4d != h2$isicrev3_4d[match(wits$hs6, h2$hs6)],
                na.rm = TRUE)))

cat("-- Reading EXIOBASE industries ...\n")
exio_ind <- read.delim(file.path(EXIO_DIR, "industries.txt"),
                       stringsAsFactors = FALSE)
names(exio_ind) <- c("exio_num", "exio_name", "exio_code", "exio_short")
cat("   ", nrow(exio_ind), "industries\n\n")


# =============================================================================
# Section 2: Hardcoded ISIC Rev 3 (4-digit) -> EXIOBASE industry mapping
# =============================================================================
# Each row maps one ISIC Rev 3 4-digit code to one EXIOBASE industry CodeNr.
# mapping_quality:
#   "exact"    — unambiguous 1:1 mapping
#   "division" — ISIC code in a division with a single EXIOBASE industry
#   "residual" — multiple EXIOBASE sub-industries possible; assigned to most
#                likely one (see inline comments)
#   "unmapped" — no meaningful EXIOBASE goods industry
#
# Organisation: by NACE Rev 1.1 / ISIC Rev 3 division (first 2 digits).
# =============================================================================

cat("-- Building ISIC Rev 3 -> EXIOBASE mapping ...\n")

isic_exio <- rbind(

  # ── Division 01: Agriculture ──────────────────────────────────────────────
  # EXIOBASE splits div 01 into 15 crop/animal sub-industries.
  # ISIC Rev 3 has only 5 codes. Many-to-one is unavoidable.
  data.frame(isicrev3_4d = "0111", exio_code = "i01.h", mapping_quality = "residual"),
    # 0111 = cereals, oil seeds, sugar, fibers, other crops -> crops nec (catch-all)
  data.frame(isicrev3_4d = "0112", exio_code = "i01.d", mapping_quality = "exact"),
    # 0112 = growing of vegetables, horticultural -> vegetables, fruit, nuts
  data.frame(isicrev3_4d = "0113", exio_code = "i01.d", mapping_quality = "residual"),
    # 0113 = growing of fruit, nuts, spices -> vegetables, fruit, nuts
  data.frame(isicrev3_4d = "0121", exio_code = "i01.i", mapping_quality = "residual"),
    # 0121 = cattle, sheep, goats, horses -> cattle farming (dominant category)
  data.frame(isicrev3_4d = "0122", exio_code = "i01.m", mapping_quality = "residual"),
    # 0122 = other animal farming -> animal products nec

  # ── Division 02: Forestry ─────────────────────────────────────────────────
  data.frame(isicrev3_4d = "0200", exio_code = "i02", mapping_quality = "division"),
  data.frame(isicrev3_4d = "0210", exio_code = "i02", mapping_quality = "division"),
    # 0210 appears in H4 concordance (variant coding for forestry/logging)

  # ── Division 05: Fishing ──────────────────────────────────────────────────
  data.frame(isicrev3_4d = "0050", exio_code = "i05", mapping_quality = "division"),
    # 0050 appears in H4 concordance (variant coding for fishing)
  data.frame(isicrev3_4d = "0500", exio_code = "i05", mapping_quality = "division"),

  # ── Division 10: Coal mining ──────────────────────────────────────────────
  # EXIOBASE i10 = "Mining of coal and lignite; extraction of peat"
  data.frame(isicrev3_4d = "1010", exio_code = "i10", mapping_quality = "exact"),
    # hard coal
  data.frame(isicrev3_4d = "1020", exio_code = "i10", mapping_quality = "exact"),
    # lignite
  data.frame(isicrev3_4d = "1030", exio_code = "i10", mapping_quality = "exact"),
    # peat

  # ── Division 11: Crude petroleum & natural gas ────────────────────────────
  # EXIOBASE splits: i11.a (crude oil), i11.b (natural gas), i11.c (other).
  # ISIC 1110 covers both crude oil and natural gas — cannot distinguish.
  data.frame(isicrev3_4d = "1110", exio_code = "i11.a", mapping_quality = "residual"),
    # crude petroleum AND natural gas -> crude petroleum (dominant by value)

  # ── Division 12: Uranium/thorium ──────────────────────────────────────────
  data.frame(isicrev3_4d = "1200", exio_code = "i12", mapping_quality = "division"),

  # ── Division 13: Metal ore mining ─────────────────────────────────────────
  # EXIOBASE splits: i13.1 (iron), i13.20.11-16 (copper, nickel, aluminium,
  # precious, lead/zinc/tin, other non-ferrous).
  data.frame(isicrev3_4d = "1310", exio_code = "i13.1", mapping_quality = "exact"),
    # iron ores
  data.frame(isicrev3_4d = "1320", exio_code = "i13.20.16", mapping_quality = "residual"),
    # non-ferrous metal ores -> other non-ferrous (ISIC 1320 covers all non-ferrous)

  # ── Division 14: Other mining ─────────────────────────────────────────────
  # EXIOBASE: i14.1 (stone), i14.2 (sand/clay), i14.3 (chemical/fertilizer
  # minerals, salt, other)
  data.frame(isicrev3_4d = "1410", exio_code = "i14.1", mapping_quality = "exact"),
    # quarrying of stone
  data.frame(isicrev3_4d = "1421", exio_code = "i14.3", mapping_quality = "exact"),
    # chemical and fertilizer minerals
  data.frame(isicrev3_4d = "1422", exio_code = "i14.3", mapping_quality = "exact"),
    # extraction of salt
  data.frame(isicrev3_4d = "1429", exio_code = "i14.3", mapping_quality = "residual"),
    # other mining nec -> chemical/fertilizer minerals (residual)

  # ── Division 15: Food, beverages ──────────────────────────────────────────
  # EXIOBASE: i15.a-d (meats by type), i15.e (veg oils), i15.f (dairy),
  # i15.g (rice), i15.h (sugar), i15.i (food nec), i15.j (beverages),
  # i15.k (fish)
  data.frame(isicrev3_4d = "1511", exio_code = "i15.d", mapping_quality = "residual"),
    # meat production/processing -> meat products nec (ISIC 1511 covers all meats)
  data.frame(isicrev3_4d = "1512", exio_code = "i15.k", mapping_quality = "exact"),
    # fish processing -> fish products
  data.frame(isicrev3_4d = "1513", exio_code = "i15.i", mapping_quality = "residual"),
    # fruit/vegetable processing -> food products nec
  data.frame(isicrev3_4d = "1514", exio_code = "i15.e", mapping_quality = "exact"),
    # vegetable/animal oils and fats
  data.frame(isicrev3_4d = "1520", exio_code = "i15.f", mapping_quality = "exact"),
    # dairy products
  data.frame(isicrev3_4d = "1531", exio_code = "i15.i", mapping_quality = "residual"),
    # grain mill products -> food nec
  data.frame(isicrev3_4d = "1532", exio_code = "i15.i", mapping_quality = "residual"),
    # starches -> food nec
  data.frame(isicrev3_4d = "1533", exio_code = "i15.i", mapping_quality = "residual"),
    # animal feeds -> food nec
  data.frame(isicrev3_4d = "1541", exio_code = "i15.i", mapping_quality = "residual"),
    # bakery products -> food nec
  data.frame(isicrev3_4d = "1542", exio_code = "i15.h", mapping_quality = "exact"),
    # sugar refining
  data.frame(isicrev3_4d = "1543", exio_code = "i15.i", mapping_quality = "residual"),
    # cocoa, chocolate -> food nec
  data.frame(isicrev3_4d = "1544", exio_code = "i15.i", mapping_quality = "residual"),
    # pasta, noodles -> food nec
  data.frame(isicrev3_4d = "1549", exio_code = "i15.i", mapping_quality = "exact"),
    # other food products nec
  data.frame(isicrev3_4d = "1551", exio_code = "i15.j", mapping_quality = "exact"),
    # distilling spirits -> beverages
  data.frame(isicrev3_4d = "1552", exio_code = "i15.j", mapping_quality = "exact"),
    # wines -> beverages
  data.frame(isicrev3_4d = "1553", exio_code = "i15.j", mapping_quality = "exact"),
    # malt liquors -> beverages
  data.frame(isicrev3_4d = "1554", exio_code = "i15.j", mapping_quality = "exact"),
    # soft drinks -> beverages

  # ── Division 16: Tobacco ──────────────────────────────────────────────────
  data.frame(isicrev3_4d = "1600", exio_code = "i16", mapping_quality = "division"),

  # ── Division 17: Textiles (single EXIOBASE industry) ──────────────────────
  data.frame(isicrev3_4d = "1711", exio_code = "i17", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1721", exio_code = "i17", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1722", exio_code = "i17", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1723", exio_code = "i17", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1729", exio_code = "i17", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1730", exio_code = "i17", mapping_quality = "division"),

  # ── Division 18: Wearing apparel (single EXIOBASE industry) ───────────────
  data.frame(isicrev3_4d = "1810", exio_code = "i18", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1820", exio_code = "i18", mapping_quality = "division"),

  # ── Division 19: Leather, footwear (single EXIOBASE industry) ─────────────
  data.frame(isicrev3_4d = "1911", exio_code = "i19", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1912", exio_code = "i19", mapping_quality = "division"),
  data.frame(isicrev3_4d = "1920", exio_code = "i19", mapping_quality = "division"),

  # ── Division 20: Wood products ────────────────────────────────────────────
  # EXIOBASE: i20 (wood products), i20.w (recycling — not in trade)
  data.frame(isicrev3_4d = "2010", exio_code = "i20", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2021", exio_code = "i20", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2022", exio_code = "i20", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2023", exio_code = "i20", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2029", exio_code = "i20", mapping_quality = "division"),

  # ── Division 21: Pulp and paper ───────────────────────────────────────────
  # EXIOBASE: i21.1 (pulp), i21.2 (paper), i21.w.1 (recycled paper — not in trade)
  data.frame(isicrev3_4d = "2101", exio_code = "i21.2", mapping_quality = "residual"),
    # pulp + paper + paperboard -> paper (ISIC 2101 covers both; paper dominates trade)
  data.frame(isicrev3_4d = "2102", exio_code = "i21.2", mapping_quality = "exact"),
    # corrugated paper containers
  data.frame(isicrev3_4d = "2109", exio_code = "i21.2", mapping_quality = "exact"),
    # other articles of paper

  # ── Division 22: Publishing, printing ─────────────────────────────────────
  data.frame(isicrev3_4d = "2211", exio_code = "i22", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2212", exio_code = "i22", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2213", exio_code = "i22", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2219", exio_code = "i22", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2221", exio_code = "i22", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2222", exio_code = "i22", mapping_quality = "division"),

  # ── Division 23: Coke, refined petroleum, nuclear fuel ────────────────────
  # EXIOBASE: i23.1 (coke), i23.2 (petroleum refinery), i23.3 (nuclear fuel)
  data.frame(isicrev3_4d = "2310", exio_code = "i23.1", mapping_quality = "exact"),
    # coke oven products
  data.frame(isicrev3_4d = "2320", exio_code = "i23.2", mapping_quality = "exact"),
    # refined petroleum products
  data.frame(isicrev3_4d = "2330", exio_code = "i23.3", mapping_quality = "exact"),
    # processing of nuclear fuel

  # ── Division 24: Chemicals ────────────────────────────────────────────────
  # EXIOBASE: i24.a (plastics basic), i24.b (N-fertiliser), i24.c (P-fertiliser),
  # i24.d (chemicals nec)
  data.frame(isicrev3_4d = "2411", exio_code = "i24.d", mapping_quality = "residual"),
    # basic chemicals -> chemicals nec (ISIC 2411 is very broad)
  data.frame(isicrev3_4d = "2412", exio_code = "i24.b", mapping_quality = "residual"),
    # fertilizers -> N-fertiliser (ISIC 2412 covers both N and P)
  data.frame(isicrev3_4d = "2413", exio_code = "i24.a", mapping_quality = "exact"),
    # plastics in primary forms -> plastics basic
  data.frame(isicrev3_4d = "2421", exio_code = "i24.d", mapping_quality = "exact"),
    # pesticides -> chemicals nec
  data.frame(isicrev3_4d = "2422", exio_code = "i24.d", mapping_quality = "exact"),
    # paints, varnishes -> chemicals nec
  data.frame(isicrev3_4d = "2423", exio_code = "i24.d", mapping_quality = "exact"),
    # pharmaceuticals -> chemicals nec
  data.frame(isicrev3_4d = "2424", exio_code = "i24.d", mapping_quality = "exact"),
    # soap, detergents -> chemicals nec
  data.frame(isicrev3_4d = "2429", exio_code = "i24.d", mapping_quality = "exact"),
    # other chemical products -> chemicals nec
  data.frame(isicrev3_4d = "2430", exio_code = "i24.d", mapping_quality = "residual"),
    # man-made fibres -> chemicals nec

  # ── Division 25: Rubber and plastics (single EXIOBASE industry) ───────────
  data.frame(isicrev3_4d = "2511", exio_code = "i25", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2519", exio_code = "i25", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2520", exio_code = "i25", mapping_quality = "division"),

  # ── Division 26: Non-metallic minerals ────────────────────────────────────
  # EXIOBASE: i26.a (glass), i26.b (ceramics), i26.c (bricks/tiles/clay),
  # i26.d (cement/lime/plaster), i26.e (other non-metallic minerals)
  data.frame(isicrev3_4d = "2610", exio_code = "i26.a", mapping_quality = "exact"),
    # glass and glass products
  data.frame(isicrev3_4d = "2691", exio_code = "i26.b", mapping_quality = "exact"),
    # non-structural non-refractory ceramic
  data.frame(isicrev3_4d = "2692", exio_code = "i26.b", mapping_quality = "exact"),
    # refractory ceramic
  data.frame(isicrev3_4d = "2693", exio_code = "i26.c", mapping_quality = "exact"),
    # structural clay products (bricks, tiles)
  data.frame(isicrev3_4d = "2694", exio_code = "i26.d", mapping_quality = "exact"),
    # cement, lime, plaster
  data.frame(isicrev3_4d = "2695", exio_code = "i26.e", mapping_quality = "exact"),
    # articles of concrete
  data.frame(isicrev3_4d = "2696", exio_code = "i26.e", mapping_quality = "exact"),
    # cutting, shaping of stone
  data.frame(isicrev3_4d = "2699", exio_code = "i26.e", mapping_quality = "exact"),
    # other non-metallic mineral products

  # ── Division 27: Basic metals ─────────────────────────────────────────────
  # EXIOBASE: i27.a (iron/steel), i27.41 (precious), i27.42 (aluminium),
  # i27.43 (lead/zinc/tin), i27.44 (copper), i27.45 (other non-ferrous),
  # i27.5 (casting). Plus recycling variants (not in trade).
  data.frame(isicrev3_4d = "2710", exio_code = "i27.a", mapping_quality = "exact"),
    # basic iron and steel
  data.frame(isicrev3_4d = "2720", exio_code = "i27.45", mapping_quality = "residual"),
    # basic precious + non-ferrous metals -> other non-ferrous (ISIC 2720 covers
    # all non-ferrous; EXIOBASE splits by metal type but ISIC does not)

  # ── Division 28: Fabricated metals (single EXIOBASE industry) ─────────────
  data.frame(isicrev3_4d = "2811", exio_code = "i28", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2812", exio_code = "i28", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2813", exio_code = "i28", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2891", exio_code = "i28", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2892", exio_code = "i28", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2893", exio_code = "i28", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2899", exio_code = "i28", mapping_quality = "division"),

  # ── Division 29: Machinery (single EXIOBASE industry) ─────────────────────
  data.frame(isicrev3_4d = "2911", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2912", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2913", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2914", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2915", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2919", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2921", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2922", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2923", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2924", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2925", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2926", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2927", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2929", exio_code = "i29", mapping_quality = "division"),
  data.frame(isicrev3_4d = "2930", exio_code = "i29", mapping_quality = "division"),

  # ── Division 30: Office machinery (single EXIOBASE industry) ──────────────
  data.frame(isicrev3_4d = "3000", exio_code = "i30", mapping_quality = "division"),

  # ── Division 31: Electrical machinery (single EXIOBASE industry) ──────────
  data.frame(isicrev3_4d = "3110", exio_code = "i31", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3120", exio_code = "i31", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3130", exio_code = "i31", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3140", exio_code = "i31", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3150", exio_code = "i31", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3190", exio_code = "i31", mapping_quality = "division"),

  # ── Division 32: Radio, TV, communication equipment (single EXIOBASE) ─────
  data.frame(isicrev3_4d = "3210", exio_code = "i32", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3220", exio_code = "i32", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3230", exio_code = "i32", mapping_quality = "division"),

  # ── Division 33: Instruments (single EXIOBASE industry) ───────────────────
  data.frame(isicrev3_4d = "3311", exio_code = "i33", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3312", exio_code = "i33", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3313", exio_code = "i33", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3320", exio_code = "i33", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3330", exio_code = "i33", mapping_quality = "division"),

  # ── Division 34: Motor vehicles (single EXIOBASE industry) ────────────────
  data.frame(isicrev3_4d = "3410", exio_code = "i34", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3420", exio_code = "i34", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3430", exio_code = "i34", mapping_quality = "division"),

  # ── Division 35: Other transport equipment (single EXIOBASE industry) ─────
  data.frame(isicrev3_4d = "3511", exio_code = "i35", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3512", exio_code = "i35", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3520", exio_code = "i35", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3530", exio_code = "i35", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3591", exio_code = "i35", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3592", exio_code = "i35", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3599", exio_code = "i35", mapping_quality = "division"),

  # ── Division 36: Furniture, other manufacturing (single EXIOBASE) ─────────
  data.frame(isicrev3_4d = "3610", exio_code = "i36", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3691", exio_code = "i36", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3692", exio_code = "i36", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3693", exio_code = "i36", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3694", exio_code = "i36", mapping_quality = "division"),
  data.frame(isicrev3_4d = "3699", exio_code = "i36", mapping_quality = "division"),

  # ── Division 40: Electricity, gas, steam ──────────────────────────────────
  # EXIOBASE splits electricity by source (12 types), plus gas distribution,
  # steam. HS codes for electricity (2716) cannot distinguish source.
  data.frame(isicrev3_4d = "4010", exio_code = "i40.11.l", mapping_quality = "residual"),
    # electricity production -> electricity nec (source unknown from trade data)
  data.frame(isicrev3_4d = "4020", exio_code = "i40.2", mapping_quality = "exact"),
    # gas manufacture/distribution

  # ── Service divisions (rarely in goods trade) ─────────────────────────────
  data.frame(isicrev3_4d = "7421", exio_code = "i74", mapping_quality = "residual"),
    # architectural/engineering activities -> other business activities
  data.frame(isicrev3_4d = "7494", exio_code = "i74", mapping_quality = "residual"),
    # photographic activities -> other business activities
  data.frame(isicrev3_4d = "9211", exio_code = "i92", mapping_quality = "residual"),
    # motion picture production -> recreational/cultural activities
  data.frame(isicrev3_4d = "9214", exio_code = "i92", mapping_quality = "residual"),
    # dramatic arts -> recreational/cultural activities
  data.frame(isicrev3_4d = "9302", exio_code = "i93", mapping_quality = "residual"),
    # hairdressing -> other service activities

  # ── Unclassified ──────────────────────────────────────────────────────────
  data.frame(isicrev3_4d = "9999", exio_code = NA_character_, mapping_quality = "unmapped"),
    # goods not elsewhere classified

  stringsAsFactors = FALSE
)

cat("   ", nrow(isic_exio), "ISIC -> EXIOBASE mappings\n\n")
stopifnot(!any(duplicated(isic_exio$isicrev3_4d)))


# =============================================================================
# Section 3: Merge chain
# =============================================================================

cat("-- Merging: HS6 -> ISIC -> EXIOBASE ...\n")

# Step 1: Keep unique HS6 -> ISIC pairs from WITS
# (some HS6 codes may appear multiple times in WITS with same ISIC)
hs_isic <- unique(wits[, c("hs6", "isicrev3_4d")])

# Step 2: Merge with hardcoded ISIC -> EXIOBASE
crosswalk <- merge(hs_isic, isic_exio, by = "isicrev3_4d", all.x = TRUE)

# Step 3: Add NACE division (first 2 digits of ISIC code = NACE Rev 1.1 division)
crosswalk$nace_div <- substr(crosswalk$isicrev3_4d, 1, 2)

# Step 4: Merge with EXIOBASE industries for names and numbers
crosswalk <- merge(crosswalk, exio_ind[, c("exio_num", "exio_name", "exio_code")],
                   by = "exio_code", all.x = TRUE)

# Step 5: Reorder columns
crosswalk <- crosswalk[, c("hs6", "isicrev3_4d", "nace_div",
                           "exio_code", "exio_name", "exio_num",
                           "mapping_quality")]
crosswalk <- crosswalk[order(crosswalk$hs6), ]
rownames(crosswalk) <- NULL


# =============================================================================
# Section 4: Validation
# =============================================================================

cat("\n-- Validation --\n")

# Check for ISIC codes in WITS that are missing from hardcoded mapping
isic_in_wits <- unique(hs_isic$isicrev3_4d)
missing_isic <- setdiff(isic_in_wits, isic_exio$isicrev3_4d)
if (length(missing_isic) > 0) {
  cat("   WARNING: ISIC codes in WITS not in hardcoded mapping:\n")
  cat("   ", paste(missing_isic, collapse = ", "), "\n")
} else {
  cat("   All", length(isic_in_wits), "ISIC codes from WITS are mapped\n")
}

# Check for NA exio_code (should only be 9999 = unmapped)
n_unmapped <- sum(is.na(crosswalk$exio_code))
n_mapped   <- sum(!is.na(crosswalk$exio_code))
cat("   Mapped HS6 codes:", n_mapped, "\n")
cat("   Unmapped HS6 codes:", n_unmapped, "\n")

# EXIOBASE industry coverage
exio_reached <- unique(crosswalk$exio_code[!is.na(crosswalk$exio_code)])
cat("   EXIOBASE industries reached:", length(exio_reached), "of", nrow(exio_ind), "\n")

# Mapping quality breakdown
cat("   Mapping quality breakdown:\n")
tbl <- table(crosswalk$mapping_quality, useNA = "ifany")
for (nm in names(tbl)) cat("     ", nm, ":", tbl[nm], "\n")

# Spot checks
cat("\n-- Spot checks --\n")
spot <- function(hs, expected_exio) {
  row <- crosswalk[crosswalk$hs6 == hs, ]
  actual <- if (nrow(row) == 1) row$exio_code else "NOT FOUND"
  status <- if (actual == expected_exio) "OK" else "MISMATCH"
  cat(sprintf("   HS %s -> %s (expected %s) [%s]\n", hs, actual, expected_exio, status))
}
spot("270900", "i11.a")   # crude petroleum
spot("720110", "i27.a")   # pig iron
spot("271011", "i23.2")   # light petroleum oils
spot("640110", "i19")     # waterproof footwear
spot("271600", "i40.11.l") # electrical energy


# =============================================================================
# Section 5: Write output
# =============================================================================

out_file <- file.path(CROSS_DIR, "hs6_exiobase_crosswalk.csv")
write.csv(crosswalk, out_file, row.names = FALSE)
cat("\n-- Output written to:", out_file, "\n")
cat("   ", nrow(crosswalk), "rows\n")

cat("\n==============================================================\n")
cat("  DONE\n")
cat("==============================================================\n")
