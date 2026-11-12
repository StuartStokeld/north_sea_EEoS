# Simulation report (layer 4.7): deferred

Not run. By scope decision (12 Nov 2026), the H2/H3 simulation study (null false-positive rate for each phase
slope, with Bonferroni and permutation; power, bias and coverage under known FP_between / FP_within slopes, on
the real panel structure and adjacency) is deferred. It will run on the **corrected** model (absolute residual,
F-A5) once the fix branches are merged, so that the calibration applies to the model actually reported.

Inputs already in place for it:

- `audit/scripts/audit_refit_helpers.R` (`audit_build_h2h3_data()`, `audit_fit_primary_car()`), which reproduces
  the live CAR fit exactly (max |Δ| 9e−18).
- Primary CAR fit time is about 0.3 s, so 1,000 null and 1,000 power fits take about 10 minutes. The
  permutation calibration is the costly part: 99 permutations × the number of simulations.
