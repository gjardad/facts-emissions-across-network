# TODO

1. Re-run run_subsampled_en.R on RMD with year-specific deployment exclusion (late ETS entrants fix)
2. Cross-validate the pre_ets backcasting procedure: hold out ETS firms observed from 2005, pretend they enter in 2008 or 2012, backcast their 2005-2007 or 2005-2011 emissions using the sector-year FE + first-two-year anchor method, and compare backcasted vs actual emissions
3. Thorough comparison of what aspects of the algorithm to assign emissions to deployment firms have been CVed in training data in chapter 1 (e.g., threshold, GPA shape, proxy weights — which are validated out-of-sample and which are assumed?)
4. Re-run b_loop_scope1_dispersion_pareto.R on RMD (with updated alloc_flags that save thr_effective, thr_initial, n_deploy_total, n_deploy_emit) and check de facto threshold across sectors and years
