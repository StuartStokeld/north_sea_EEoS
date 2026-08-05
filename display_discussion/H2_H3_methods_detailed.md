## 1. Overview

H2 and H3 are proposed to share a single model, rather than being tested as two
separate analyses on two separate pre-aggregated datasets. H2 tests whether fishing
pressure explains prediction error *across space*; H3 tests whether it explains
prediction error *across time*. Both questions can be asked of the same underlying
data: individual hauls, individual rectangles, individual years :without pre-
aggregating into fixed zones or fixed time periods first. The exploratory analysis
(summarised in Section 4) is what motivates this: every fixed-aggregation approach
tried loses real information, in a way a shared hierarchical model avoids.

## 2. Spatial unit of analysis

**Unit: the individual haul, nested within its individual ICES rectangle.** No
spatial merging (block or pressure-tier zones) is applied.

**Model component**: rectangle enters as a random intercept (`(1 | rectangle)`),
allowing each rectangle its own baseline level of prediction error, but with
**partial pooling** — rectangles with few hauls are pulled toward the overall
average more than rectangles with many hauls, rather than being treated as equally
reliable. Spatial autocorrelation between neighbouring rectangles (the reason the
original H2 analysis needed a spatial error model, SEM) is incorporated directly into
this random-effects structure — e.g. a spatial correlation function over rectangle
centroids (a conditional autoregressive or Gaussian-process structure), rather than a
separate post-hoc spatial correction step. This is intended to do the same job the
original SEM did, but without first requiring data to be aggregated to a fixed zone.

Fishing pressure enters as two separate terms, not one. Rectangle-year fishing
pressure (`log(fishing_hours + 1)`, its natural resolution in the Couce dataset) is
decomposed into a **between-rectangle** component and a **within-rectangle**
component:

- `FP_between` — each rectangle's own mean fishing pressure across the whole study
period (time-invariant, one value per rectangle). This is the genuinely *spatial*
term, and is what answers **H2**: does a rectangle's persistent fishing-pressure
level predict its residual.
- `FP_within` — that rectangle-year's deviation from its own rectangle's mean
(time-varying, mean zero within each rectangle by construction). This is the
genuinely *temporal* term, and is what answers **H3**: does a rectangle's own
year-to-year fluctuation in fishing pressure track its own residual.

**Why the decomposition is necessary.** A single blended `log(fishing_hours + 1)`
term — the original approach — mixes both sources of variation into one coefficient.
Roughly 59% of the variance in fishing pressure is between-rectangle (persistent
differences) and 41% is within-rectangle (year-to-year fluctuation), per the ICC
established during the exploratory phase (Section 5). A model with only the blended
term cannot say whether an observed fishing-pressure effect reflects "persistently
high-fishing rectangles differ from low-fishing ones" (H2) or "a rectangle's own
fishing pressure moving up or down tracks its own residual moving up or down" (H3) —
and in principle the two could even point in different directions. This is the
within-between (or "Mundlak") conflation problem, standard in panel-data modelling
(Bell & Jones, 2015, "Explaining Fixed Effects," *Political Science Research and
Methods*). The decomposition above resolves it: each hypothesis now has its own
coefficient, rather than one shared coefficient standing in for both.

**Feasibility outcome**: three versions of this spatial structure were tested — no
spatial term (plain random intercept), a continuous distance-decay correlation
(range/sill over rectangle centroids), and a CAR structure on a discrete
queen-adjacency matrix reused directly from the original H2 SEM's own weights object.
The continuous version was not usefully identified (estimated range exceeded the
data's maximum inter-centroid distance; CIs spanned orders of magnitude) — a
data-type mismatch, since continuous distance-decay functions suit point-referenced
data, not the areal (lattice) structure of ICES rectangles. The CAR version converged
cleanly and its autocorrelation parameter was well-identified, but sits near the edge
of its admissible range — which corroborates, rather than contradicts, the original
H2 SEM's own strong autocorrelation estimate (λ = 0.858), obtained independently via a
different spatial model family. Across all three versions, the coefficients that
matter for H2/H3 (fishing pressure, phase, their interaction) were essentially
unchanged; only nuisance parameters (intercept, biomass coefficient) shifted.
**Decision**: the plain `(1 | rectangle)` random intercept is adopted as the primary
model (simplest, no identifiability caveats to defend), with the CAR and continuous
versions reported as sensitivity analyses demonstrating the result does not depend on
the choice of spatial structure.

*Note: these three feasibility comparisons, and the initial results run, used the
single blended `log(fishing_hours + 1)` term, predating the within-between
decomposition above. Removing a fixed effect or splitting one term into two does not
introduce new identifiability problems, so the convergence and partial-pooling
findings are not expected to change — but the decomposed model's H2/H3 coefficients
have not yet themselves been through the full feasibility/sensitivity sequence. Given
how mechanical this split is (arithmetic decomposition of an existing term, not a new
model family), a full re-run of every sensitivity check is not treated as required
before proceeding — see Section 6.*

*Separately: these same three feasibility comparisons were run with `mean_ln_B_obs`
included, as that was the specification at the time. Biomass has since been removed
from the primary model (see above) — again not expected to change convergence or
partial-pooling behaviour, but these checks did not specifically test "with vs.
without biomass," which remains part of the deferred exploratory work in Section 6.*

## 3. Temporal unit of analysis

**Unit: the individual haul's own year**, categorised into a small number of
phases rather than aggregated into a single before/after split or tested at the
level of individual years.

**Phase structure (primary):** policy-anchored breakpoints at **1992 / 2002 / 2008**,
giving four phases — 1985–1991 (pre-reform baseline), 1992–2001 (1992 CFP reform),
2002–2007 (2002 CFP reform, in force Jan 2003), 2008–2015 (2008 LTMP / MSFD).
Implemented as `phase_v2` in the live primary model.

**Phase structure (sensitivity / design history):** a Bai–Perron structural-break
analysis on the whole-area mean fishing-hours series located **3 statistically
robust breaks (1989, 2001, 2008)**. A candidate 4th break (1997) had only marginal
statistical support and is not used. Those data-driven phases informed early design
and remain available for comparison; they are **not** the presented primary.

**Model component**: phase enters as a fixed-effect categorical variable. The key test
for H3 is the **interaction between `FP_within` and phase**: does the effect of a
rectangle's own year-to-year fishing-pressure fluctuation on its own residual change
across phases — this is the version that genuinely tests within-rectangle temporal
covariation, now that fishing pressure has been decomposed (Section 2). Because
rectangle already has its own random intercept, and `FP_within` is by construction a
deviation from that rectangle's own mean, this interaction cannot be explained by
persistent differences between rectangles — it isolates change within a rectangle
over time, conceptually similar to the paired before/after comparison discussed
earlier, but handled by the model's structure rather than by manually differencing
pre-aggregated zone means. `FP_between`'s own interaction with phase is fitted
alongside it as a secondary check — whether the *spatial* (H2) pattern itself differs
across phases — but the primary H3 test is the `FP_within × phase` term specifically.

**Robustness outcome**: the categorical phase structure was compared against a
continuous linear year term and a varying-coefficient GAM (allowing a fully smooth
fishing-pressure-by-year effect). All three specifications show a statistically
significant fishing-pressure × time interaction. The phase model's step-function
estimate of the fishing-pressure effect correlates at ~0.78 with the GAM's smooth
estimate (77% same-sign agreement across years), and both capture the same broad
post-2002 rise in the presented phase windows; the linear-only term misses a later
downturn that both the phase structure and the GAM pick up. **Decision**: a
categorical phase structure is retained as the primary temporal specification — it is
a reasonable simplification of the smoother underlying pattern, and with
policy-anchored `phase_v2` it stays directly interpretable against reform history in a
way a smooth term is not. *(This check, like the spatial feasibility checks, used the
blended fishing-pressure term and the earlier data-driven phases — it validates phase
vs. continuous time as the temporal structure, not the within/between split or the
final breakpoint choice specifically; not treated as needing a full re-run for the
same reason given in Section 2.)* The GAM is kept as a supplementary robustness
figure showing the effect isn't an artefact of discretising time.

## 4. Shared model — summary

A single model tests both hypotheses, using the within-between decomposition of
fishing pressure (Section 2):

- **H2 (spatial)**: significant `FP_between` main effect on residual, with
rectangle-level spatial autocorrelation accounted for.
- **H3 (temporal)**: significant `FP_within × phase_v2` interaction — whether a
rectangle's own year-to-year fishing-pressure fluctuation's association with its
own residual differs across the four policy-anchored time phases.

Both are read off the same fitted model rather than two separately-built analyses on
two separately-aggregated datasets, and each now has its own dedicated coefficient
rather than sharing one blended term.

## 5. Why this design was chosen — evidence from the exploratory analysis

**Spatial choice.** The pre-H3 feasibility check showed only 1 of 197 rectangles had
enough repeat annual sampling for a rectangle-level temporal test on its own (median
2 hauls/rectangle/year) — ruling out treating rectangles as independently analysable
units without some form of pooling. Two fixed-aggregation alternatives were tested as
a fix:

- Geographic block merging (blind to fishing pressure) reduced the real, meaningful
rectangle-to-rectangle variation in fishing pressure (ICC dropped from 0.59 at
individual-rectangle resolution to 0.38–0.45 at coarser block sizes).
- Fishing-pressure-based zones preserved fishing-pressure variation better (ICC
0.52–0.53) but preserved *residual* variation worse than the geography-blind blocks,
and produced highly uneven zone sizes (many zones reduced to a single rectangle).

No fixed zone scheme protected both variables' real structure at once, because
fishing pressure and prediction error aren't distributed across the map in the same
pattern — already evident from the original H2 finding that fishing pressure didn't
robustly predict residual once spatial autocorrelation was properly modelled. A
hierarchical model with partial pooling avoids this trade-off: it uses every
individual rectangle's real data without discarding resolution, and handles sparse
rectangles by borrowing strength from the overall pattern rather than forcing them
into an arbitrary group.

**Temporal choice.** The original plan (annual resolution) was ruled out by the same
feasibility check. A single before/after split at 2003 was the first proposed fix,
anchored to the 2002–2004 cod crisis reforms. The structural-break analysis tested
whether this framing — and the more detailed 4-phase visual pattern — was genuine
rather than an artefact of eyeballing the series. BIC comparison gave decisive support
for 3 breaks (large, consistent improvements at each step up to 3 breaks; only a
marginal 1.2-point improvement for a 4th; BIC increases for 5 or more), and the 3
robust breaks (1989, 2001, 2008) closely matched the original visual boundaries. The
2001 break's confidence interval (2000–2002) independently overlaps the 2002–2004
reform window, without that policy history being an input to the statistical method —
useful corroboration that a mid-study transition is real in the fishing-hours series.

The **presented primary** nevertheless uses policy-anchored `phase_v2` (1992 / 2002 /
2008) so phase labels map directly onto CFP / LTMP–MSFD milestones. The Bai–Perron
breaks remain design support and a natural sensitivity alternative; refitting the
within-between model on `phase_v2` left the headline null result unchanged (marginal
R² still ≈ 0.051; no single uniform H2/H3 disturbance signal).
