No, not quite. BPPR and the PBX rerun are the only two with work ready to start. The E8 brief still lists other open items, though most are your calls or were deliberately left for later.

**Ready to start**
- **BPPR (E36):** `issues/regression/bass/bppr-brief.md` is waiting for your approval.
- **PBX 9501 rerun:** the data is in `artifacts/datasets/PBX9501`, and you deferred it.

**Your calls, from slice 6's report**
- **The matrix check's default:** it's opt-in now. With the cache it still costs 1.3–2× the readout for scalar and functional fits, and 3–5× for basis fits.
- **A K × K matrix check for basis fits:** it would cut the basis fit's check to about the scalar's cost, which could change the default question above.
- **A public `outputLocations` on `BassFit.Gradient`:** one optional parameter; the internal path already takes it.

**Model changes, also your call:** a directional hinge.

**Deferred outside this brief**
- A robust median of subspaces as the default.
- A closed-form C for BGP.
- G4's joint mixture of products, and E25's readouts (ALE, Shapley).

**Settled since they were listed**
- Rotate-and-refit was built as slices 7a–7c. The brief's "Not in this brief" line still names it, which is stale.
- `DenseEigen`'s absolute threshold: the relative Jacobi rule (`75523c4`) replaced it, and `ScaledEigen` is gone.
- The truncated-normal constant and the matrix check's cache were built in slice 6.
- How projected fits handle column-mass support and the raw ridge: slice 7d built it, and you confirmed the defaults today.

--- 


I've recorded both changes in the brief's open-items list (commit `f978a86`):
- **PBX 9501:** deferred to your broader pass over the applications in the papers behind the regression engine.
- **The five unbuilt items:** these were the ones I'd listed as "not asked for". They're now marked as items you intend to commission, after the other loose ends and before BPPR:
  1. Co-active output weightings at output locations beyond the training ones.
  2. Poincaré constants for Student t, empirical and mixture laws.
  3. Cross-fitting with groupings or case weights.
  4. A held-out log predictive density in the projection comparison.
  5. Basis fits through a projection.

Two of these touch other decisions:
- **Item 5:** a projected refit of a basis fit is a model change, like rotate-and-refit was. Its plan would need your calls, for example whether all components share one frame.
- **Item 2:** the mixture and Student t constants interact with slice 6's open question about the matrix check's default. These constants extend which laws the check covers. They don't change its cost.

---

I'm reading the remaining code paths: the Poincaré constants, the input measures, case weights, and how basis fits reach the projection.

I have what I need for all five items. Writing slice 8's plan into the brief now.

The plan is written and checked against the code, and I'm committing it now.

Another session committed at the same moment; retrying.

Slice 8's plan is drafted in the E8 brief (commit `a6da53b`), just above "Gates common to every slice". The brief's status line and open-items list link to it. Nothing is built: it waits for you to commission it and make the decisions below. Each decision has a recommendation.

**8a: two readout items**
- **A. Co-active output weightings at new locations.** These would cover direct functional fits only. Basis fits refuse new locations, because their basis exists only at the training locations.
  - **A1, the form:** I recommend two new properties, `OutputLocations` and `OtherOutputLocations`, beside the existing weights, the same way `OutputMetric` takes its locations.
  - The alternative pairs weights and locations in one record. That would replace `Output` and `OtherOutput`, which makes it a removal and its own question.
- **B. Poincaré constants for the remaining laws.**
  - The existing bounds hold only for log-concave laws, and mixtures and Student t laws generally aren't log-concave. A two-peaked mixture can even exceed the uniform law's constant.
  - Every such law here is truncated to a bounded interval, so a constant still exists.
  - **B1:** I recommend solving for it numerically on a fine grid. The error is far below Monte Carlo error, and the readout would state the grid. The alternatives are a valid but loose bound, under which the check would rarely fire, or leaving these laws unchecked.
  - **B2:** I recommend keeping slice 6's rule for a single truncated normal, so no existing reading moves.
  - The empirical law has no such constant: its point masses let a function change between training values with zero derivative on them. It stays unchecked, and the readout says why. That part of the item closes as a finding.

**8b: two comparison items**
- **C. Held-out log predictive density.** It reuses the conformal chip's existing predictive-density code without editing it. That code already handles every noise law, case weights and random effects.
  - **C1:** I recommend scalar responses only for now. That code takes one observation of a scalar response, so a functional fit's held-out density would be refused with that reason.
- **D. Cross-fitting with groupings and case weights.**
  - **D1, how cases are split into folds:** I recommend splitting by case by default, with an opt-in that keeps each level of a named grouping in one fold. That opt-in scores prediction for unseen groups.
  - **D2, scoring under case weights:** I recommend following what the weights mean. Replication weights count each held-out case w times. Precision weights enter only the density, and each case counts once in the RMSE.

**8c: basis fits through a projection.** This changes the model, as slice 7 did.
- **E1:** I recommend one frame shared by every component, read from the basis fit's own analysis. The components then stay on one set of inputs, so readouts and Sobol indices over the reconstruction still work.
- **E2:** I recommend reusing the output basis unchanged for the full-data refit. Under cross-fitting, each fold re-estimates its basis from the training folds.
- **E3:** I recommend a new `BassProjectedBasisFit` type beside `BassProjectedFit`, rather than widening the existing type.

Each sub-slice has its own gates in the plan, including a recovery test against a function from outside the model and a seeded sweep. The plan runs 8a, then 8b, then 8c, each verified before the next. You can settle 8c's decisions at the same time as the rest, so it won't hold anything up.