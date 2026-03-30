# TODO

1. Copy locally results for all years from b_loop_scope1_dispersion_pareto from RMD
2. ~~Investigate why mixed_capped is so high~~ → DONE: replaced min(ETS) constraint with 30kt physics-based cap; added NACE 38 to energy CRF group. See DISTRIBUTION_FITTING.md sections 9.1–9.4.
3. Compare results with Lyubich et al (2018) and De Lyon & Dechezlepretre (2025) benchmarks
4. ~~Write a detailed document explaining the Pareto allocation mechanism~~ → DONE: see DISTRIBUTION_FITTING.md sections 5.1–5.5.
5. ~~Investigate the energy/electricity sector~~ → DONE: identified CRF-NACE mapping mismatch (CHP, waste incinerators). Fixed by adding NACE 38 to energy. Documented in DISTRIBUTION_FITTING.md section 9.4.
6. Re-run run_subsampled_en.R on RMD with year-specific deployment exclusion (late ETS entrants fix)
7. Implement pre_ets backcast in b_loop_scope1_dispersion_pareto.R: backcast emissions for firms entering ETS after the current year using sector-year fixed effects and first-two-year anchor
8. Implement 30kt physics-based cap in b_loop_scope1_dispersion_pareto.R: replace min(ETS) constraint with fixed 30kt upper bound
9. Update nace_crf_crosswalk.csv on RMD (NACE 38 → energy)
10. Cross-validate the pre_ets backcasting procedure: hold out ETS firms observed from 2005, pretend they enter in 2008 or 2012, backcast their 2005-2007 or 2005-2011 emissions using the sector-year FE + first-two-year anchor method, and compare backcasted vs actual emissions
