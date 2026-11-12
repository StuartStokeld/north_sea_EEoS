# Mutation testing report (layer 4.8, scoped to 4 mutants)

Script: `audit/scripts/mutation_test.R`. Results: `audit/outputs/mutation/mutation_results.csv`. Run 12 Nov 2026.

## Method

1. Copy `pipeline/` to a temporary directory and inject one bug into the copy. The repo is never touched.
2. Run the whole suite (`pipeline/tests/testthat`) against the copy, via `EEOS_PIPELINE_DIR`.
3. A mutant is **killed** when at least one test that passes on the reference copy fails on the mutant.

The current code already contains two of the mutants (LngtCode "5" read as mm; signed residual). For those, the
reference copy is the **corrected** code (fixed in the temporary copy only), and the question is whether the
suite would catch the bug being reintroduced.

## Results: kill rate 4 / 4

| ID | Mutant | Reference | Killed | Killing test(s) |
|---|---|---|---|---|
| M1 | Remove the C × HaulDur/60 conversion (`datras_hh_helpers.R:141`) | current code | yes (2) | `test-integrity.R` accounting (rebuilt N, E ≠ outputs); `test-metamorphic.R` C ↔ R equivalence |
| M2 | Treat LngtCode "5" as mm | ICES-correct midpoints | yes (1) | `test-units-hand.R` [F-A2] midpoint |
| M3 | Signed instead of absolute residual (H2 CAR formula) | absolute-residual formulas | yes (1) | `test-units-hand.R` [F-A5] script formulas |
| M4 | FP_within not centred (`h2h3_within_between_helpers.R`) | current code | yes (1) | `test-units-hand.R` within–between decomposition |

## Caveats

- **M2, M3 and M4 are each killed by a single test.** The pre-existing tests (`test-datras_*.R`,
  `test-h1_*.R`, `test-h2_panel.R`) load code with `source("../../R/...")`, not through `EEOS_PIPELINE_DIR`. So
  they always test the repo copy and can't kill mutants. Moving them to `source_pipeline()` would add
  redundancy.
- **M3 is caught by a static check of the script formulas.** The model-object test (response = |residual|)
  fails on both copies, because the outputs aren't regenerated in a mutation run.
- **A harness bug was found and fixed during the run.** An empty reference key made M1 and M4 compare against
  nothing, so the first run reported 2/4. After the fix, M1 and M4 were re-run.
- **Scope.** The other brief mutants (duration filter, phase boundary, log base, permutation +1, join
  duplication) weren't run by scope decision. The suite has direct tests for each:
  - phase boundaries, natural log, permutation p: `test-units-hand.R`
  - join cardinality and duration accounting: `test-integrity.R`
