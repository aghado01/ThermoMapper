# Conformal validity in the Lean harness — chip brief

CB-1 to CB-4 were proved on 2026-11-06, E35's M-1 and M-2 on 2026-11-08, and E16 slice 9's jackknife+ and asymmetric CQR on 2026-11-09. Their reports are below. E16 slices 10 and 12 were commissioned and proved on 2026-11-10, and the Lean for slices 13 to 16 and 18 was commissioned the same day.

## Slice 10's grid reading, and the Lean for E16 slices 13 to 16 and 18 (2026-11-10)

Status: commissioned 2026-11-10 by the user ("lets also send the pending lean work to conformal lean chip"), relayed by the main session.

**Why now.** These slices are Lean first (decision 7 on the [layer plan](../../../issues/regression/conformal-layer-brief.md#e16s-completion-the-plan-slices-7-onward)). The C# order of 13 to 18 waits for slice 10's BASS reading (decision 10). The Lean runs ahead of it, in the plan's order, so whichever slice the user picks first has its statements ready.

**How.**
- The method, environment, gates and reporting are as for slices 10 and 12 below: checks first, exhaustive at small n, then the Lean.
- A statement that proves false, or holds only with a further hypothesis, is a finding.
- One report per item, under "Reports". Tell the requesting session tersely when each lands.
- No `docs/` edits.

**Sources.** Before starting each item, check its sources are on disk: `gauntlet/conformal-bayes` with its `inventory.jsonl`, and `PDFs/Bayesian`.
- List any missing sources in one message to the main session, which asks the user to fetch them.
- Meanwhile, work on the items whose sources are present.
- Read a paper before stating what it says.

**Order:**
1. **Slice 10's grid reading, formalized.** Slice 10's C# is next in the C# queue, and it reports two things the realized-ordering section checks but does not formalize:
   - the margin over the grid: the largest kept reading, the largest |D(y)| over every grid value, and 2ν;
   - the PIT's linearity between grid values.

   Prove them, or record exactly what reading the maximum over the grid needs.
2. **Slice 13, the tournament correction** (Bhattacharyya, Zhang and Barber, arXiv 2605.29200, whose theorems the CB-3 prior art already read).
   - *The statements:* coverage at least 1 − 2α whatever the draws, from leave-one-out fits each symmetric in its data; and the 2ε-inflated form's 1 − α − 3√ν under their stability.
   - *What the C# may claim:* the symmetry each leave-one-out chain needs (it treats its data symmetrically but need not mix), and the comparison rule's exact indices and ties under CB-1's convention.
3. **Slice 14, shift.**
   - *Covariate shift:* weighted conformal with a stated likelihood ratio (Tibshirani, Barber, Candès and Ramdas 2019), and the loss with an estimated ratio.
   - *Beyond exchangeability:* Barber, Candès, Ramdas and Tibshirani's weights and their coverage-gap term (2202.13415).
   - *Label shift.*
   - *The check:* at equal weights on exchangeable data, CB-1 is recovered.
4. **Slice 15, approximate conditional coverage:** Zhang and Candès's posterior conformal prediction. Marginal coverage exact; the per-cluster statement with its hypotheses.
5. **Slice 16, batches** (Lee, Tchetgen Tchetgen and Dobriban): sets for functions of a batch, simultaneous sets, and selection with false claims controlled. The finite core is m test units placed uniformly among n + m (`powersetCard`).
6. **Slice 18, the guarantee read out:**
   - Snell and Griffiths's posterior over the loss a calibration set implies, for small m, and its calibration over replicated calibration sets;
   - proper scores (Gressmann et al.) where a finite propriety statement fits.

## E16 slices 10 and 12 — the full step's error, Hoff's score and the latent posterior (2026-11-10)

Status: commissioned 2026-11-10 by the user ("Yes, send it"), relayed by the main session. Slice 10 and slice 12 were proved on 2026-11-10 and reported under Reports.

These are Lean first, by the user's decision 7 on the [layer plan](../../../issues/regression/conformal-layer-brief.md#e16s-completion-the-plan-slices-7-onward), so the C# for slices 10 and 12 waits for these statements. Do slice 10 first and report it on its own, so its C# can start while slice 12's Lean runs. The method, environment, self-review and reporting are as set out below. The numeric checks come first, exhaustive at small n with exact arithmetic. Then come the Lean, `#print axioms` on the main results, the module build, `scripts/meta-ci.ps1 -NoBuild` and D:'s free space. The files are the three named under Scope, plus `conformal-searches.py` where a search needs it.

### Slice 10 — the full step's error, completed

Read first, in the [layer brief](../../../issues/regression/conformal-layer-brief.md):
- the slice 4 build report, for `ConformalDistribution` and the CDF score F̃;
- the per-unit bounds' report, for `MarginBounds.Carried`, the orderings' bounds and why each dominates.

- **CB-3 for the CDF score.** CB-3 is proved for the accept rule on a density score at one level. The full step's opt-in distribution reads the CDF score F̃ under the add-one-in weights: Q(y, τ) = (B + τ(T + 1))/m, and the equi-tailed set takes k = ⌊αm/2⌋ per tail. Pose and prove what estimated scores within εⱼ of the exact ones cost:
  - the equi-tailed set, in the margin form and the inflated form, per tail and two-sided, and whether the two tails' bands add or share units;
  - the distribution itself: a bound on P(PIT ≤ u) − u from the fragile counts, at each u or uniformly in u;
  - whether the band read from any ordering's estimates dominates, as `card_fragile_le_estimated` does for the density score.
- **The carried bound that leaves out the degenerate orderings.** This is CB-3's open follow-up, and it closes the ledger's open "margin under an ε that holds only with high probability".
  - *The problem:* the layer's carried bound takes each unit's worst error over the orderings. A unit far out degenerates the weights in its own turn as the candidate, and that widens every unit's bound. On BASS the margin's band holds 96 to 100 units of 100.
  - *What to pose:* bounds over the orderings outside an excluded set D, the turns whose weights degenerate (read by k̂ or the ESS). Each bound holds with probability 1 − ν. The coverage bound subtracts D's share and ν from the margin read under the kept orderings' bounds, of the form 1 − (k + |M_D|)/m − |D|/m − ν if that holds.
  - *Also settle:* the share's exact form, and whether the inflated rule has the counterpart.
  - *Compare with* Bhattacharyya, Zhang and Barber's 1 − α − 3√ν (arXiv 2605.29200, appendix).
  - A different or weaker form is a finding.
- *Decisions served:* whether the full step's distribution reports its Monte Carlo error, and whether the margin can read anything on fits like BASS's.

### Slice 12 — Hoff's score and the latent posterior

Read first: CB-4 and its Hoff finding in the [ledger](conformal-lemmas.md), and Hoff 2023 §4.1 (`PDFs/Bayesian/Hoff2023.pdf`).

- **Hoff's leave-one-out score,** cᵢ = p(zᵢ | every other unit, the candidate included), in CB-4's finite setting:
  - its equivariance, so that CB-1 covers with it, under CB-1's tie rule (the ledger records that Hoff's strict rule fails under ties);
  - the identity the C# computes from the n-unit draws. Each data unit takes the weight g(z_new | θ)/g(zᵢ | θ), add-one-in and leave-one-out at once, and c_{n+1} is the plain n-unit predictive. State the estimator's error in a form slice 10's bounds can carry;
  - under latents, Hoff's score with every latent integrated (CB-4's (b)) and by the latent posterior ((c)). Say whether each stays equivariant, and whether reading the units differently breaks it as (d) does;
  - Hoff's optimality: the finite statement that ranking by cᵢ gives the smallest expected set, under the model's own law, among sets valid at that level. Prove it if it fits the finite style; otherwise record what it needs. It is what would make the option worth its cost, so a counterexample is a finding.
- **CB-4's latent-posterior reading from draws.** The exact score (c) is proved (`latentScore_equivariant`, `latentScore_eq`, `sum_jointU_addOneIn`).
  - *The C# estimator:* from the n-unit draws, the data units read each draw's retained latents. The candidate's latent is drawn from p(v′ | θ, z_new), as a second candidate matrix.
  - *To prove:* the estimator's expectation over the drawn latents is (c)'s score for every unit, the candidate included.
  - *To state:* what the drawn latent adds to the Monte Carlo error that slice 10's bounds must carry.
- *Decision served:* whether the full step offers Hoff's score, and reads a conditioned observation latent.

Out of scope, as under Scope below, and slices 13 to 18. Report under "Reports", newest first, one report per slice, and tell the requesting session when each is done.

## E16 slice 9 — cross-fitted intervals and asymmetric CQR (2026-11-09)

Status: done 2026-11-09, report below.

This is the Lean prerequisite for the [conformal layer's slice 9](../../../issues/regression/conformal-layer-brief.md#e16s-completion-the-plan-slices-7-onward). Establish the finite-sample statements the implementation may claim:

- **Jackknife+.** Read the jackknife+ theorem of Barber, Candès, Ramdas and Tibshirani against the proposed rule. State the exchangeability and fitting-symmetry hypotheses, the exact order-statistic indices and endpoint conventions, and prove the marginal coverage lower bound of at least 1 − 2α. State which hypotheses and bound apply to CV+ if you include it. The engine brief's near-1 − α coverage is an empirical measurement.
- **Asymmetric conformalized quantile regression (CQR).** With a fixed fit independent of the exchangeable calibration and test units, calibrate the lower and upper scores separately at α/2. State each tail's miscoverage bound, how ties and small calibration samples affect the endpoints, and the resulting two-sided coverage bound.

Use the existing [protolemma ledger](conformal-lemmas.md) and [small-case checks](conformal-checks.py) to pin down ranks and ties before proving the statements in [ConformalValidity.lean](../../Enthymemata/ConformalValidity.lean). Run `lake build Enthymemata.ConformalValidity` and the existing meta-CI gate from the environment section. Report theorem names, assumptions, any counterexamples or weaker bounds, and the exact rank rules the C# worker can implement. The layer plan separately describes the later Lean work for slices 10 and 12.

## Context

Conformal Bayes is engine work: one layer over retained draws and pointwise densities, built into GBASS and the regression engine generally ([E16](../../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer), the user's call on 2026-11-05). The layer owns coverage and the score owns width. Four statements pin the layer's contract before it is built, and each one changes a decision, named below. A statement that turns out false, or true only with a further hypothesis, is a finding worth more than a proof.

The harness already has the sibling: `GuardedEstimator` (GE-1 to GE-3), "estimators that stay valid whatever the emulator does". Conformal is the same idea for prediction sets: valid whatever the model does. Use its finite-sum style, and its `powersetCard` and simple-random-sample machinery where they fit.

**This chip's job:** pose CB-1 to CB-4, check them numerically, then prove them as compiling enthymemata, and report what each fixes in the layer's contract.

**Read first:**
- [lean/README.md](../../README.md): the taxonomy, meta-CI's rules, import discipline.
- [guarded-estimator-lemmas.md](../BASS/guarded-estimator-lemmas.md) and `lean/Enthymemata/GuardedEstimator.lean`, for the style, the namespace and the finite sums.
- E16 in [`extensions-brief.md`](../../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer), including "What the conformal literature adds".
- Fong and Holmes, *Conformal Bayesian Computation* (`gauntlet/conformal-bayes/2106.06137v2`, TeX source), for the add-one-in step, its exchangeability statement and its remark on Monte Carlo error. Hoff 2023 (`PDFs/Bayesian/Hoff2023.pdf`), its conformal section, for the rank construction from the values cᵢ.
- `docs/Maths/VarianceMeanMixtures.md` and the pointwise section of `docs/Maths/BayesianRegression.md`, for how GBASS's pointwise density integrates each latent (OQ17's scope rule).

## Environment

- **Toolchain:** as pinned in `lean/lean-toolchain` and `lake-manifest.json`. Don't bump them.
- **Mathlib's build output may be absent** from `lean/.lake`.
  - Fetch by scoped module only, as `lake exe cache get <Module>`. A bare `lake exe cache get` restores everything; leave that out.
  - D: has been nearly full. Check `df -h /d` before each fetch and after, and **if a fetch would leave under 3 GB free, stop and report.** Delete nothing to make room: disk decisions are the user's.
  - Build this chip's module (`lake build Enthymemata.ConformalValidity`), and run `scripts/meta-ci.ps1 -NoBuild` for the taxonomy and ledger gates. Run the full `lake build` only if the disk allows, and say which ran.
- **Iterate fast:** call `lean` directly with `LEAN_PATH` from `lake env printenv LEAN_PATH`. Linters run only under `lake build <Module>`, so finish with that.
- **Quirks:**
  - `Π` is a reserved token, and identifiers with a combining macron don't parse.
  - `omit … in` and `open … in` go before the docstring.
  - `lean/` is LF, so write files from Python with `newline='\n'`.
  - Meta-CI's ledger counts the word for an apology in comments too, so keep it out of docstrings.
- **Papers:**
  - arXiv sources: `D:/aghado01/codex-scientiae/supellex-private/gauntlet/<folder>/<id>/<id>-tex/`. The conformal ones are in `conformal-bayes`, martingale posteriors and Syring–Martin in `bayesian`.
  - PDFs: `…/PDFs/Bayesian/`, with an `inventory.jsonl` by title.
  - `pdftotext` works, and `pypdf` keeps math glyphs.
- **Python:** `C:/Users/azrie/PDenv/PyenvPython/pyenv-win/versions/3.12.6/python.exe`, run with `-I -X utf8` on a script file. The bare `python` shim breaks on `|` in arguments.
- **The machine and the checkout are shared with other sessions.**
  - Run builds at below-normal priority.
  - Commit by pathspec (`git commit -- paths`), leave nothing staged, work on main, push nothing.
  - Commit messages describe the work, with no attribution or Co-Authored-By lines ([CLAUDE.md](../../../CLAUDE.md)).
  - Touch only the files named here. Leave the user's untracked notes under `issues/` and their modified `AGENTS.md` and `CLAUDE.md` alone.

## Scope

New files:
- `lean/Protolemmata/Conformal/conformal-lemmas.md`: the statements, with a status line;
- `lean/Protolemmata/Conformal/conformal-checks.py`: the numeric checks;
- `lean/Enthymemata/ConformalValidity.lean`: the proofs, added to `Enthymemata.lean`. The namespace is `Conformal`, since these serve every engine, unless the README's taxonomy says otherwise.

For each statement, the protolemma gives the statement, its hypotheses, the decision it changes and its prior art. **The numeric checks come before any Lean.** Make them exhaustive enumeration over every ordering at small n, exact rational arithmetic where it fits, plus CB-4's counterexample. Take the statements below as a first draft: fix the conventions (conformity or nonconformity, ≤ or <, the tie rule) once, in CB-1, and use them throughout.

Order, easiest and most load-bearing first, so an early stop still leaves value:

1. **CB-1, the rank guarantee.**
   - *The statement:* take n + 1 units, exchangeable, and a score that is permutation-equivariant (permuting the units permutes the scores). Then the test unit's rank among the n + 1 scores is uniform, and the accept rule at level k covers with probability at least 1 − k/(n + 1). With no ties, or with random tie-breaking, it covers with at most 1 − k/(n + 1) + 1/(n + 1).
   - *The finite form:* exchangeability as the average over all (n + 1)! orderings of a fixed multiset.
   - *The per-stratum form:* units exchangeable within strata give coverage in each stratum.
   - *Decision:* the layer's tie rule, the coverage band the C# tests assert on enumerable fixtures, and what the steps of 1/(k + 1) mean for small strata.
2. **CB-2, the split form holds for any fit.**
   - *The statement:* given the training split, any map from training data to a score function leaves the calibration and test units exchangeable, so CB-1 applies.
   - *The instance:* CQR's score max(q_lo(x) − y, y − q_hi(x)).
   - *Decision:* the split form as the guarantee that does not lean on the sampler, cited as a theorem.
3. **CB-4, latents keep the score symmetric only if every unit treats them alike.**
   - *The finite setting:* θ in a finite Θ, each unit's latent v in a finite V, prior π(θ), latent law q(v | θ), density f(y | x, θ, v), and the integrated density g(y | x, θ) = Σ_v q(v | θ) f(y | x, θ, v).
   - (a) *Add-one-in:* p(θ | Z₁…Zₙ₊₁) ∝ p(θ | Z₁…Zₙ) · g(yₙ₊₁ | xₙ₊₁, θ).
   - (b) *The integrated score* sᵢ = Σ_θ p(θ | Z₁…Zₙ₊₁) · g(yᵢ | xᵢ, θ) is permutation-equivariant.
   - (c) *The latent-posterior score*, sᵢ = E[f(yᵢ | xᵢ, θ, vᵢ)] under the joint posterior of θ and all n + 1 latents, is permutation-equivariant too. Computed from draws given n units, it needs the candidate's latent drawn from p(v | θ, yₙ₊₁) for each draw.
   - (d) *The mixed score breaks symmetry:* data units scored by (c), with latents from the n-unit posterior, and the candidate by g. Give a numeric counterexample, and a Lean one if it comes cheaply.
   - *Decision:* which density the conformal layer reads for every latent-noise model (Student t, quantile, Normal-Wald, generalized beta prime, and Pólya–Gamma later). It is the contract to fix before that interface is built.
4. **CB-3, what an imperfect score costs.**
   - *The deterministic part:* if every estimated score is within ε of a permutation-equivariant exact score, the k-th order statistic moves by at most ε. So the test's decision is unchanged unless |s_test − s₍ₖ₎| ≤ 2ε.
   - *The probabilistic part:* by exchangeability, P(|s_test − s₍ₖ₎| ≤ 2ε) = E[#{i : |sᵢ − s₍ₖ₎| ≤ 2ε}]/(n + 1), so coverage falls at most that far below CB-1's bound. Keep the ties consistent with CB-1.
   - *The engine's reading:* the count can be bounded from the estimated scores, those within 4ε of the estimated threshold. State that as a corollary if it holds.
   - *Prior art:* check whether this bound, or a sharper one, is known. Barber, Candès, Ramdas and Tibshirani, "Conformal prediction beyond exchangeability" (*Annals of Statistics*, 2023; arXiv 2202.13415), bound coverage gaps by total-variation terms, and there may be stability results for conformal methods. These are not on disk, so record what you find from abstracts and say so.
   - *Decision:* whether the full-conformal layer reports its margin, and when it falls back to the split form. ε has to come from the sampler's own error estimate, which is only as good as its mixing.

### Out of scope

- Any change to `src/` or `tests/`.
- Hoff's optimality, CBMA, posterior conformal prediction's approximate conditional coverage, weighted conformal under shift, and batch predictive inference. Where a step suggests one, record it as a finding. Batch's finite core, m test units placed uniformly among n + m, is `powersetCard` again.
- Promotion to `Lemmas/`: proof-closed enthymemata wait for semantic review.
- New design, and any measurement on a fixture.

## Self-review

After each statement, take a fresh-eyes pass before moving on, for correctness and for coverage gaps, and fix what it finds. Does the Lean say what the protolemma claims, with no hypothesis doing the work? Is the exchangeability hypothesis the one the layer can actually meet? Run `#print axioms` on the main results. Name the passes and their fixes in the report.

## Reporting

Append a report below, newest first, under "Reports", and keep the protolemma's status line current. It should say:
- what is proved, by declaration name, and what remains open;
- every statement that changed, why, and any that proved false;
- what each statement fixes in the layer's contract, as a line E16 can take. Edit E16's "What the conformal literature adds" only where a result corrects it, and say where;
- prior art found for each statement, read or from abstracts;
- which build gates ran, and the disk's free space at the end;
- follow-ups only with the decision each would change. Record a question with no decision waiting as open, not proposed;
- what stopped you, if anything did.

## E35: conformal test martingales (2026-11-08)

The user asked for this on 2026-11-08 ("Yes, start it"), relayed by the BASS extensions session. [E35](../../../issues/regression/bass/extensions-brief.md#conformal-test-martingales-exchangeability-checked-online-e35) tests exchangeability online. Its code comes after the conformal layer's slice 5, and these statements fix its contract first.
- **M-1, the smoothed p-values.** Under exchangeability, with θ₁, θ₂, … i.i.d. uniform and independent of the data, p₁, p₂, … are i.i.d. uniform. The rank part is finite. The statement makes three premises explicit: the random tie-break, ranking against every earlier score (calibration and test alike), and exchangeability of the whole sequence.
- **M-2, Ville's inequality for the betting product** Sₙ = Π f_i(p_i).
  - Each f_i is ≥ 0, integrates to 1 on [0, 1], and is chosen from the earlier p-values.
  - The link is its own lemma: the product is a martingale under i.i.d. uniform p-values.
  - P(supₙ Sₙ ≥ c) ≤ 1/c, and the capital at a stopping time is an e-value.
- **Optional:** pose the CUSUM and Shiryaev–Roberts guarantees.
- **Method as for CB-1 to CB-4:**
  - checks first: exhaustive at small n, and a seeded Monte Carlo of Ville's bound attained;
  - Lean in the finite style where it fits, with standard axioms only;
  - the module and barrel builds, meta-CI, and D:'s free space.
- **Sources:**
  - `gauntlet/conformal-bayes`: 1204.3251, 1706.03415, 2102.10439;
  - `PDFs/Bayesian`: `Vovk2003.pdf`, `Vovk_StatSci2021.pdf`;
  - text extractions in the requesting session's scratchpad, their math checked against the formulas used.

## Reports

Newest first.

### 2026-11-10 — slice 14: shift, proved

The statements are in the ledger, under "Slice 14". The new declarations are PROOF-CLOSED, with the standard three axioms only. The checks ran first.

- *The weighted quantile lemma.* The refused units carry weight at most α, for any nonnegative weights. At equal weights the rule is CB-1, which was the commission's check.
- *Covariate and label shift, per data set* (`coverage_weighted`): with the likelihood ratio at the units, the refused weight is at most α of the total.
- *The loss with an estimated ratio* (our statement; the papers' result is asymptotic): at most α + TV(p, p̂) of the data set's normalized weights.
- *Beyond exchangeability* (their Theorem 2, on a finite law): P(refused) ≤ α + Σ_K q_K · TV(P, P ∘ swap). This form uses the TV of the data laws; their score-vector form is not formalized.

### 2026-11-10 — slice 13: the tournament correction, proved

The statements are in the ledger, under "Slice 13". There are 15 new declarations, 330 in the file, PROOF-CLOSED, with the standard three axioms only. The checks ran first.

- *Theorem 1 is slice 9's cross-fitted rule* with the pair fit `tournFit` (the fit to the rest, plus the pair as an unordered multiset of values).
  - Their threshold, (1 − α)(n + 1) strict wins, is the rule at k = ⌊αm⌋ (`tournament_count_iff`).
  - Coverage is ≥ 1 − (2k − 1)/m with fits symmetric given their seeds, up to relabeling (`tournament_coverage`), or with chains that read their units in index order and need not mix (`tournament_orderReading`).
- *What the C# needs.* The pair reaches the fit as values. A fit that marks which unit is the test refuses it on every data set in the checks.
- *Theorem 2, per data set* (`coverage_tournament_stable`): for every a, the share accepted is ≥ 1 − (k + a + 2F/(a + 1))/m.
  - *Finding:* the reference R can be any scores of the data set, not only their one-out score.
  - *Finding:* optimized over a, the cost is about α + 2√(2ν), below their α + 3√ν: 0.227 against 0.300 at m = 20, ν = 0.01.

**Gates.** The module build is clean, and all 9 main results print the standard axioms.

### 2026-11-10 — item 1: slice 10's grid reading, proved

The statements are in the ledger under slice 10, "The grid reading". There are 11 new declarations, 315 in the file, PROOF-CLOSED, with the standard three axioms only. The checks ran first, and every assertion holds.

- *The PIT between grid values.* It is linear: at (k + θ)/m it equals (1 − θ) P(R ≤ k) + θ P(R ≤ k + 1), per ordering and averaged (`pitCdf_interp`, `avg_pitCdf_interp`). So sup_x |P(PIT ≤ x) − x| is the largest departure at the grid values (`interp_dev_le`).
- *The margin over the grid* (`coverage_grid_excluding`): coverage is at least 1 − (average over orderings of [largest kept reading + largest flagged count over the grid])/m − 2ν.
- *Finding: the maximum needs the grid to hold the value the test unit takes.*
  - The grid may be built from the data, for example over y at the test unit's own x.
  - On a continuous response, a fixed grid's maximum is a reading at its points, not a bound.
  - A sound reading needs a grid value in each interval where the counts are constant. Their breakpoints are where estimate comparisons or flags change, which are not slice 8's crossings of π(y).
  - In the checks, a grid that missed one value missed its reading in 222 of 1,620 arrangements.

**Gates.** The module build is clean, the barrel is green, and meta-ci −NoBuild passes. D: had 38 GB free.

### 2026-11-10 — slice 10's follow-up: D's rule, and reading from the realized ordering

The main session verified slice 10 (dfe5c99) and raised two points about the engine's reading. Both are settled in the ledger, under "Reading from the realized ordering" in the slice 10 section. There are 5 new declarations, 304 in the file, PROOF-CLOSED, with the standard three axioms only. The checks ran first, and every assertion holds.

**1. D is each unit's own turn, the candidate's included.** Your reading is right.
- D(y) is read per (data set, y) over all m units, by one per-turn rule.
  - The candidate's own turn is the held ordering at y, read by its own held weights.
  - A data unit's own turn is the ordering in which it is the candidate.
- A D from the data units' turns alone never contains the test unit, and it fails `coverage_close_excluding`'s hypothesis.
- The engine reads data units' own-turn flags only through a proxy.
  - A proxy is sound if its set contains every unit whose own turn degenerates. Extra flags only charge more.
  - A proxy that misses one is not sound.
  - At a far y, a data unit's own turn has the far candidate inside its fit. So the proxy must not inherit the candidate's held weight.

**2. At a flagged y the realized ordering reads nothing, and the reading moves to the kept orderings with D charged again.**
- At a flagged y, every unit's estimate shares the candidate's add-one-in weights, so no unit can be read there. Charging the candidate as one unit of D and reading the rest from that ordering is not sound.
- *Proved per data set* (`le_read_excluding`, `read_le_excluding`, `coverage_read_excluding`, `coverage_read_twoSided_excluding`, `pit_read_excluding`). Each ordering reads its own estimates, and a left-out ordering reads nothing.
  - *Coverage:* at least 1 − (average reading)/m − |D|/m − 2ν.
  - *The PIT:* between (average lower reading)/m − 2ν and (average upper reading)/m + |D|/m + 2ν.
  - D is charged once inside the count and once for the orderings that read nothing.
  - The second ν is the realized ordering's own bounds failing.
- *Over the grid* (as the layer already reads, not formalized):
  - the margin is 1 − max over kept y of ĉ(y)/m − max over all grid y of |D(y)|/m − 2ν;
  - the PIT's bounds are [min over kept y of (lower)/m − max_y |D(y)|/m − 2ν, max over kept y of (upper)/m + max_y |D(y)|/m + 2ν].
- *The maximum does not collapse because the grid reaches far out.* Flagged y drop out of the maximum, and a far y costs its one flagged unit in max_y |D(y)|. The report collapses only if the flagged count is large at some grid value.
- *Other readings.* A reading per y is not a coverage statement. Reading over the candidates outside D is sound only with D's share charged.

**Rules for the C# worker,** replacing slice 10's "report from the observed ordering's estimates":
- *D(y):* all m units. The candidate's own turn is read directly from its held weights. Each data unit is read through a proxy that flags at least whenever its own turn would.
- *Read at kept y only:* the band from ĉ(y) = k + |M̂_D| + |D's units the inflated rule accepts| (2k and the union band for the equi-tailed set), and the PIT's counts likewise.
- *Charge* max over the grid of |D(y)|/m, and 2ν in place of ν.

**Files.**
- `ConformalValidity.lean`, section `ReadingExcluded`.
- `conformal-checks.py`, slice 10, "reading from the realized ordering".
- `conformal-lemmas.md`: the slice 10 section and the status line.
- This brief.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- `#print axioms` on all 5: the standard three.
- D: had 38 GB free.

**Open, no decision waiting:** how the C# proxies a data unit's own-turn flag. Its soundness can be checked on fixtures against refits.

**What stopped me:** nothing.

### 2026-11-10 — slice 12: Hoff's score and the latent posterior, proved

Slice 12 is proved in `ConformalValidity.lean`: 60 new declarations, 299 in the file, PROOF-CLOSED, with the standard three axioms only. The ledger's [slice 12 section](conformal-lemmas.md#slice-12--hoffs-score-and-the-latent-posterior) has the statements. The checks ran first, in CB-4's model and in exact arithmetic, and every assertion holds.

**1. Hoff's leave-one-out score.**
- *Valid.* It is equivariant, and CB-1 covers with it, ties counted for the candidate (`hoffScore_equivariant`, `hoff_coverage`). The strict rule covers 0 when every unit ties.
- *The draw identity.*
  - The candidate's score c₍ₙ₊₁₎ = A = E_n[g(z_new | θ)] is the plain n-unit predictive (`hoffScore_none`).
  - A data unit's score is cᵢ = A/Bᵢ, with Bᵢ = E_n[g(z_new | θ)/g(zᵢ | θ)], the mean of its add-one-in and leave-one-out weight (`hoffScore_some`).
  - The candidate is accepted iff at least k data units have Bᵢ ≥ 1 (`hoff_accepts_iff`).
- *The error in slice 10's form.* Â is common to every unit of an ordering. Ranks see only differences, so CB-3 and every slice 10 form hold with bounds up to one shift per ordering (`coverage_close_shift`, `coverage_inflated_shift`). On the log scale the bounds are |log B̂ᵢ − log Bᵢ| for the data units and 0 for the candidate.
- *Latents.*
  - (b), every latent integrated, is `hoffScore`.
  - (c) read literally is the joint posterior given the others. It is (b), since unit i's latent is then at its prior (`hoffJointScore_eq`).
  - (c) read as the draws read it conditions each unit's latent on its own value (`hoffCondScore`). It is equivariant and sits above (b) unit by unit. Its draw identities: the data units read their retained latents with Hoff's weights, and the candidate draws its latent per draw, unweighted.
  - Reading the units differently breaks it, as CB-4(d) does: coverage 0 against 1/2 (`hoffMixed_counterexample`), and 0.53 against 0.75 in the checks' model.
- *Hoff's optimality: proved, finite* (`hoff_expectedSize_le`, `hoff_cb1_expectedSize_le`).
  - The setting is any exchangeable positive law and any base measure, so covariates are covered.
  - The claim: among procedures (any φ, randomized, not necessarily symmetric) that cover at least as much on every data set, a threshold in Hoff's score has the smallest expected size.
  - With distinct scores, CB-1's rule is exactly 1 − k/m, so the comparison is with every valid procedure. `hoffScore_eq_hoffP` ties `hoffP` to the model's score.
- *Finding: ties.* CB-1's rule counts ties for the candidate. So where Hoff's scores tie at the boundary, his set over-covers, and another score's set can be smaller. In the checks' model, CB-1's Hoff set is 3.719 against 3.678 for the in-sample score at m = 3, k = 1.
- *How much it buys.* At exact coverage (ties broken at random), Hoff's set is smallest, but only 0.2–0.5% below the in-sample (Fong and Holmes) score's. It is 1–4% below (c)'s.

**2. CB-4(c) from draws, the candidate's latent drawn per draw.**
- *Proved:* the estimator's expectation is (c)'s score for every unit, the candidate included (`latentScore_some_draws`, `latentScore_none_draws`).
- *What the drawn latent adds.* Over the drawn latents (`drawnLatent_mean_var`), the candidate's estimate has:
  - mean equal to the Rao–Blackwell reading Σ_s w̄_s E[f | θ_s, z_new];
  - variance exactly Σ_s w̄_s² Var(f | θ_s, z_new).

  That variance is what each unit's carried bound must take from its own turn as the candidate, and a closed-form conditional mean removes it. The common normalizer drops out as a shift.

**Rules for the C# worker.**
- *Hoff's score:*
  - each data unit's weight is g(z_new | θ)/g(zᵢ | θ);
  - the candidate gets the plain n-unit predictive, with no weight;
  - accept iff #{B̂ᵢ ≥ 1} ≥ k;
  - carry bounds on log B̂ᵢ, and 0 for the candidate.
  - The weights' k̂ is the place to read slice 10's D.
- *Any latent reading:* every unit reads alike.
  - For the conditioned reading, the candidate's latent is drawn from p(v′ | θ, z_new) for each draw, unweighted under Hoff's score and with the add-one-in weight under CB-4(c).
  - Prefer Σ_v q f²/g in closed form when the model has it.

**Decisions served.**
- *Does the full step offer Hoff's score?*
  - For: valid and optimal at equal coverage.
  - Against: the gain over the in-sample score is small here, and a per-unit weight brings PSIS instability and degenerate turns. Under CB-1's tie rule on discrete data the gain can reverse.
  - The measurements are at m ≤ 4 in one model. The decision is yours.
- *Does it read a conditioned observation latent?* For Hoff's score there is nothing to read, since the literal latent posterior is (b). The conditioned reading is valid only if every unit reads alike.

**Files.**
- `ConformalValidity.lean`, sections `ShiftedBounds`, `Hoff`, `DrawnLatents`, `HoffOptimal`, `HoffLink`.
- `conformal-checks.py`, section "Slice 12".
- `conformal-lemmas.md`: the slice 12 section, and the status line, whose count is now 299 with the Hoff finding closed.
- This brief.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- `#print axioms` on all 49 new theorems: the standard three, or fewer.
- Builds ran at below-normal priority. D: had 38 GB free.

**Open, no decision waiting:**
- Hoff's gain at larger n and in a continuous model.
- The law forms on continuous spaces, unchanged.

**What stopped me:** nothing.

### 2026-11-10 — slice 10: the full step's error, proved

Slice 10 is proved in `ConformalValidity.lean`, section `FullStepError`: 45 new declarations, 239 in the file, PROOF-CLOSED, with the standard three axioms only. The ledger's [slice 10 section](conformal-lemmas.md#slice-10--the-full-steps-error-completed) has the statements. The checks ran first, and every assertion holds.

**1. CB-3 for the CDF score F̃.**
- *The equi-tailed set.*
  - Margin form: at least 1 − |R ∪ M_lo ∪ M_hi|/m, with |R| ≤ 2k (`coverage_close_twoSided`).
  - Inflated form: at least 1 − 2k/m (`coverage_inflated_twoSided`).
  - *Do the bands add or share?* They enter as their union, so they add unless they share a unit. A unit can be fragile in both tails only when at least m − 2k + 2 units sit within its bounds, which is more than (1 − α)m at k = ⌊αm/2⌋ (`card_near_of_fragile_both`).
  - The union band read from any ordering's estimates dominates (`card_fragileTwoSided_le_estimated`).
- *The distribution.* At each grid value,

  |P(PIT ≤ k/m) − k/m| ≤ max(#{L < k} − k, k − #{U ≤ k})/m.

  - L and U are the counts shifted by the bounds (`pit_le_upper`, `pit_ge_lower`, `card_shiftedAtMost_le_le_card_shiftedBelow`).
  - With distinct exact scores, the excesses are the PIT's fragile units at level k.
  - The bound holds uniformly in u, since P(PIT ≤ u) is linear between grid values. That step is checked, not formalized.
  - The counts read from any ordering's estimates, with doubled shifts, dominate (`card_pit_counts_estimated`).

**2. The carried bound that leaves out the degenerate turns.**
- *The form holds, and its share is exact* (`coverage_close_excluding`). Coverage is at least 1 − |D ∪ R ∪ M_D|/m − ν.
  - D's share is |D \ (R ∪ M_D)|/m: a degenerate turn whose unit the exact rule already refuses costs nothing.
  - The brief's 1 − (k + |M_D|)/m − |D|/m − ν follows, and is never tighter.
- *The inflated counterpart:* at least 1 − |D ∪ R|/m − ν (`coverage_inflated_excluding`).
- *Added, as the full step needs them:*
  - the equi-tailed set with D left out, both forms;
  - the PIT with D left out (`pit_excluding`): #{u ∉ D : U ≤ k}/m − ν ≤ P(PIT ≤ k/m) ≤ #{u ∈ D or L < k}/m + ν.
- *What the engine reads from one ordering whose bounds hold* (`card_excluded_le_estimated` and its two-sided form, `card_pit_excluding_estimated`). The share is at most k (or 2k) + |M̂_D| + |{u ∈ D : the inflated rule accepts u}|. A D unit the inflated rule refuses is refused by the exact rule too.
- *Against Bhattacharyya, Zhang and Barber's 1 − α − 3√ν:*
  - Theirs takes ν per score.
  - Ours takes ν per ordering, jointly over the m units, and pays ν once. The layer's z = Φ⁻¹(1 − ν/(2m)) already buys the joint ν by a union bound.
  - For a cost of 0.01, theirs needs ν ≈ 1.1·10⁻⁵ per score (z ≈ 4.39). Ours takes ν = 0.01 (z ≈ 3.89 at m = 100, 4.42 at m = 1000).
- *An illustration at m = 100, k = 5, ν = 1/20,* with one far unit in D and in R. The margin reads 0.00 with bounds over every turn, and 0.88 with that turn left out. The bounds there are chosen; BASS's own reading is the C#'s to measure.

**Rules for the C# worker.**
- *Read D per turn, once for the data set,* from k̂ or the ESS of that turn's weights. It must not depend on which unit is the test unit.
- *Make ν joint per ordering,* through the existing per-unit z. Carry εⱼ over the kept turns only.
- *Report from the observed ordering's estimates:*
  - 1 − (k + |M̂_D| + |D ∩ inflated-accepted|)/m − ν for the margin;
  - 2k, the union band, and D's units that the inflated set accepts, for the equi-tailed set;
  - 1 − |D ∪ R|/m − ν, or at least 1 − (k + |D|)/m − ν, for the inflated rule;
  - for the PIT, the counts #{u ∈ D or L̂ < k} and #{u ∉ D, Û ≤ k}, with the shifts doubled, ± ν.
- Without D, the same rules hold with D empty and ν = 0 for the deterministic forms.

**Decisions served.**
- *The full step's distribution can report its Monte Carlo error.* Use the PIT band at each grid level and the union band for the set. Whether it reads anything on real fits is the C#'s measurement.
- *The margin can read something on BASS-like fits* when the degenerate turns are far units that the exact rule refuses.

**Files.**
- `ConformalValidity.lean`.
- `conformal-checks.py`, section "Slice 10".
- `conformal-lemmas.md`: the slice 10 section, and the status line, whose stale count of 187 is now 239.
- This brief.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- `#print axioms` on all 38 new theorems: the standard three.
- The full checks run in about 5 minutes. Builds ran at below-normal priority, and D: had 38 GB free.

**Open, no decision waiting:** the PIT's uniform-in-u step in Lean. It needs a real-valued CDF of the randomized rank.

**What stopped me:** nothing. Slice 12 is next.

### 2026-11-09 — slice 9's open items: what a fit may see, and CV+ on a fixed bag

You asked to address the two open items. The ledger has the record: "What a fit may see" under the jackknife+, and "CV+ by one fold" plus "On a fixed bag" under CV+. `ConformalValidity.lean` is 194 declarations, PROOF-CLOSED, with the standard three axioms only.

**1. Does the jackknife+ need a symmetric fit? The answer depends on what the fit sees.**
- *Its units' order: no symmetry needed* (`jackknifePlus_orderReading`, through `OrderReading`, `acceptsPairwise_cycleRange` and `coverage_orderReading`).
  - A fit that reads its training units in index order, with the left-out row's place closed up, keeps 1 − (2k − 1)/(n + 1).
  - The proof moves any unit to the front with the others kept in order, so the test unit is strange exactly as often as any other.
  - Exact MILPs had shown this worst case up to m = 6 before the proof.
- *Its units' positions: no guarantee.* The parity construction: swapping the test unit with unit i flips the ordering's parity, and a fit that knows which positions are missing can side with the even completion. Coverage is exactly 1/2 at every level, for every n.
  - It is checked exactly in `conformal-checks.py` for m = 3 to 7, with a real interval.
  - It passes 1 − 2α at k = 1 from m = 5 on.
- *The C# contract.* Each leave-out fit gets its training rows in index order, with the left-out row's place closed up. It gets no row index, no mask, and no seed keyed to a row. Then a sampler sweeping in index order with a fixed seed is covered, with no random presentation order needed.
- *CV-a, the merged bound,* holds even for fits that know positions.

**2. CV+'s tournament bound on a fixed bag of n + 1 units: partly resolved.**
- *Proved:* `cvPlus_oneFold`, P(refused) ≤ ℓ/(m + 1) on any bag. That settles K = 2, and m = 1 is the jackknife+.
- *Two-fold CV+ covers above 1 − 2α with no extra term,* since ℓ/(m + 1) < 2α. That is better than both of Barber et al.'s bounds at K = 2.
- *Settled at K = 3, m = 2, ℓ ≤ 2:* exactly by MILP at ℓ = 1 (1/7 against 1/4), and by a Markov count at ℓ = 2 (1/2).
- *Searched for K = 3, 4, 5 and m = 2, 3 at every level, with no counterexample.* The worst found stays well inside CV-b's bound.
- *Still open for K ≥ 3,* with no decision waiting: on a fixed bag the layer asserts the merged and one-fold bounds.

**Rules for the C# worker, updated.**
- The CV+ guarantee to assert is the largest of three:
  - the merged bound, 1 − 2(ℓ + K − 1)/(n + K);
  - one fold, 1 − ℓ/(m + 1);
  - the tournament, 1 − (2ℓ + m − 2)/(n + m), for a symmetric fit on i.i.d. data.
- At K = 2 that is the one-fold bound.

**Files.** A new file [`conformal-searches.py`](conformal-searches.py) holds the searches, annealing and an optional scipy MILP. This chip had not named it; it keeps the open item's evidence reproducible. `conformal-checks.py` keeps the standard-library assertions: the parity certificate, the one-fold bound over every CV+ ordering, and the K = 2 and Markov arithmetic.

**Coordinator replay.** The optional MILP in the delivered script checks the jackknife+ cases. The fixed-bag optima can now be replayed separately with `python -I -X utf8 conformal-searches.py --exact-bags`. For two folds of two, strict-order witnesses attain the one-fold upper bounds of 1/3, 2/3 and 1. For three folds of two at ℓ = 1, a seven-point counting argument gives 1/7 and a strict-order witness attains it. The [ledger](conformal-lemmas.md) gives the argument. These are finite certificates in Python.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- Builds ran at below-normal priority. D: had 38 GB free.

**Open, no decision waiting:** CV-b on a fixed bag for K ≥ 3.

**What stopped me:** nothing.

### 2026-11-09 — slice 9's follow-up: CV+ proved, both bounds

You fetched Barber, Candès, Ramdas and Tibshirani 2021 (1905.02928v3) on 2026-11-09, and I read it in TeX. CV+ is now proved in `lean/Enthymemata/ConformalValidity.lean`, which is 187 declarations and still PROOF-CLOSED. `#print axioms` on the CV+ results shows the standard three only. [conformal-lemmas.md](conformal-lemmas.md) has the record, in its CV+ section.
- *CV-a, by merging the folds' ranks* (Vovk; Theorem 4(a)).
  - `cvPlus_merged`: coverage at least 1 − 2(ℓ + K − 1)/(n + K). `cvPlus_merged_alpha` puts it at least 1 − 2α − 2(1 − α)(K − 1)/(n + K).
  - Through `merged_le`, `sum_antitone_le`, `sum_max_sub_le`, `rank_superUniform`, `crossFold_refused_le` and `sum_rankLe_lt_of_cvRefuses`.
- *CV-b, by a tournament across folds* (Theorem 4(b)).
  - `cvPlus_tournament`: coverage at least 1 − (2ℓ + m − 2)/(n + m).
  - Through `card_refused_across_le`, `ratio_prodPerm`, `coverage_across` and `not_acceptsAcross_of_cvRefuses`.

**The rules the C# worker can implement.**
- *CV+ at level α.*
  - Set ℓ = ⌊α(n + 1)⌋ and split the n training units into K folds of m = n/K.
  - Fit μ̂₋ₖ without fold k. Rᵢ is |Yᵢ − μ̂₋ₖ₍ᵢ₎(Xᵢ)|, under unit i's own fold's fit.
  - The interval is closed, from the ℓ-th smallest of μ̂₋ₖ₍ᵢ₎(x) − Rᵢ to the ℓ-th largest of μ̂₋ₖ₍ᵢ₎(x) + Rᵢ, the jackknife+'s ends with fold fits. At ℓ = 0 it is the whole line.
- *The guarantee to assert* is the larger of the two bounds:
  - 1 − 2(ℓ + K − 1)/(n + K) from the merged ranks, for any fit that never sees its own fold or the test unit;
  - 1 − (2ℓ + m − 2)/(n + m) from the tournament, for a symmetric fit on i.i.d. data.
  - At α = 0.1, n = 1000: K = 5 guarantees 0.79 (merged), K = 10 0.78 (merged), K = n 0.80 (the jackknife+).
- *Fits.* The merged bound needs no symmetric fit, so a sampler that reads its units in index order, with a fixed seed, is covered for K-fold CV+. The tournament bound and the jackknife+ need the symmetry, or the shared random order.

**Correction to the earlier slice 9 report.** It called 1 − (2k − 1)/(n + 1) "sharper than 1 − 2α", as if new. The 2k − 1 count is Barber et al.'s own: their Step 1 reaches 2α(n + 1) − 1 before Theorem 1 rounds it off. Their Theorem 2 also shows the factor of 2 is real, asymptotically. What is new here is the exact finite attainment. The ledger now says so.

**What changed from the task.**
- CV+ was a follow-up, and both of Barber et al.'s bounds are proved, in counts.
- The merged bound uses no randomization and no level adjustment. Barber et al. reach it through Vovk's modified cross-conformal set at an adjusted α′.
- CV-b's hypothesis is stronger than exchangeability of n + 1 units. It needs n + m units exchangeable under permuting folds and slots, which i.i.d. data give.

**Checks first** (`conformal-checks.py`, slice 9, CV+). All assertions hold.
- the merging inequality on 21,142 cases, and the layer cake;
- CV-a over 51,600 orderings, including two fits that read their units' order;
- CV-b over 4,392 orderings of n + m units;
- CV-b's counting core on 551,934 cross-fold outcomes. It is not attained for 2 folds of 3, where the worst is 4 against 5;
- both bounds inside Barber et al.'s α forms on 8,400 cases.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- Builds ran at below-normal priority. D: had 38 GB free.

**Open, no decision waiting:** whether CV-b holds over the orderings of a fixed n + 1 units, without the further test units. CV-a does.

**What stopped me:** nothing.

### 2026-11-09 — slice 9: the jackknife+ and asymmetric CQR proved

**Proved** in `lean/Enthymemata/ConformalValidity.lean`, now 153 declarations and still PROOF-CLOSED. `#print axioms` on the slice 9 results shows the standard three only. [conformal-lemmas.md](conformal-lemmas.md) has the record, in its slice 9 section.
- *The jackknife+:*
  - `card_refused_pairwise_le`, the counting core;
  - `coverage_pairwise_orderings`, over orderings;
  - `coverage_pairwise` and `crossFit_coverage`, under a law and for any score;
  - `not_acceptsPairwise_of_jackRefuses`, the interval contains the cross-fitted set;
  - `jackknifePlus_orderings`, `jackknifePlus_coverage` and `jackknifePlus_two_alpha`;
  - `acceptsPairwise_zero`, `pairScore_equivariant` and `one_sub_two_mul_le` as tools.
- *Asymmetric CQR:*
  - `cqrAsym_coverage`, each tail and both;
  - `cqrAsym_coverage_half`, at α/2 per tail;
  - `cqrLo_isUpperSet`, `cqrHi_isLowerSet` and `cqrAsym_ordConnected`, each tail's side;
  - `accepts_zero` and `sum_and_ge` as tools.

**The rules the C# worker can implement.** n training units, m = n + 1. Each rule cites its Lean declaration.
- *Jackknife+ at level α.*
  - Set k = ⌊α(n + 1)⌋. Refit μ̂₋ᵢ without unit i, for each i, and take Rᵢ = |Yᵢ − μ̂₋ᵢ(Xᵢ)|.
  - The interval is closed, from the k-th smallest of μ̂₋ᵢ(x) − Rᵢ to the k-th largest of μ̂₋ᵢ(x) + Rᵢ. A y at an end is covered.
  - At k = 0 the interval is the whole line; the layer says so, as at small strata in CB-1.
  - As counts (`JackRefuses`): y is refused when at least n + 1 − k upper ends are strictly below it, or at least n + 1 − k lower ends strictly above it.
- *The guarantee to assert:*
  - coverage at least 1 − (2k − 1)/(n + 1) (`jackknifePlus_coverage`), which is at least 1 − 2α (`jackknifePlus_two_alpha`);
  - on enumerable fixtures, over every ordering (`jackknifePlus_orderings`).
  - The bound is attained, so no test may assert more. At k = 1 it is 1 − 1/(n + 1), the split form's exact rate.
- *The cross-fitted set with any score* (`AcceptsPairwise`, `crossFit_coverage`), Vovk's n-fold cross-conformal set.
  - The test candidate is refused when it is stranger than at least n + 1 − k units. Unit i counts when the candidate's conformity under μ̂₋ᵢ is strictly below unit i's own conformity under μ̂₋ᵢ, and a tie does not count.
  - The same guarantee holds. For the residual score, the set lies inside the jackknife+ interval.
  - This is the form that takes conformal Bayes' predictive densities, since it is not tied to |y − μ̂|.
- *The fits* (`FitSymmetric`).
  - Every leave-out fit reads its training units as a set.
  - A sampler qualifies if it receives its units in one uniformly random order drawn once per run and shared by every leave-out fit, or if its seeds travel with the units. That is the law form's relabeling `act`.
  - A sampler that sweeps units in index order with a fixed seed is outside the statement.
  - A fit may depend on its sample size and need not be stable. The plain jackknife, centred on the fit to all n units, has no guarantee, and its coverage is 0 in the checks with such a fit.
- *Asymmetric CQR at α/2* (`cqrAsym_coverage`, `cqrAsym_coverage_half`).
  - Set k_lo = k_hi = ⌊(α/2)(n₂ + 1)⌋.
  - The lower end is q_lo(x) + c_lo(k_lo), where c_lo(k) is the k-th smallest of yᵢ − q_lo(xᵢ) over calibration. The upper end is q_hi(x) − c_hi(k_hi), from q_hi(xᵢ) − yᵢ. Both are closed, with ties counted for the test, as CB-1's rule.
  - At level 0 that end is infinite (`accepts_zero`), which happens below n₂ + 1 = 2/α: n₂ < 19 at α = 0.1.
  - Assert each tail missing at most k/m, exactly k/m on distinct scores, and two-sided coverage at least 1 − (k_lo + k_hi)/m. With equal tail levels, the excess of this lower bound over 1 − α lies in [0, 2/m). Its excess over the symmetric finite-rank bound 1 − ⌊αm⌋/m is either zero or 1/m.

**What changed from the task.**
- *The jackknife+'s bound is 1 − (2k − 1)/(n + 1),* sharper than 1 − 2α by at least 1/(n + 1), and attained.
  - A symmetric fit built from a regular tournament covers exactly that over every ordering, for every (m, k) checked up to m = 7.
  - At α = 0.1 and n = 19 the worst case is 0.85.
  - Nothing proved false.
- *The law form is stated up to a relabeling of the fits' randomness,* so the "fitting-symmetry hypothesis" is one a sampler can meet in law. The task asked only for symmetry.
- *CV+ is not included,* and neither are cross-conformal predictive distributions. The reasons are in the ledger; CV+ is a follow-up below.
- *The engine brief's "near 1 − α" is a measurement under stability, not a guarantee,* as the assignment said.

**Checks first** (`conformal-checks.py`, section "Slice 9"). All assertions hold.
- *The jackknife+:*
  - the counting core on 92,577 comparison outcomes, with 2k − 1 attained;
  - the order statistics, Barber et al.'s quantiles and the counts agreeing on 304,500 cases;
  - 4,920 orderings over 4 data sets and 5 symmetric fits;
  - exact laws, i.i.d. and Pólya, at m = 4;
  - the attaining fit at 10 values of (m, k).
- *Each premise broken:* the plain jackknife, and ties counted as strange.
- *Order-reading fits:* exhaustive at m = 4, with none beyond the bound.
- *Asymmetric CQR:*
  - each tail and both under an i.i.d. law at every pair of levels;
  - each tail exact over orderings;
  - 16,384 accept sets matching the closed ends and Romano et al.'s quantile;
  - symmetric CQR putting all its miscoverage above the set (0.363), against 0.127 per tail for the asymmetric form.

**Prior art, read.**
- Barber, Candès, Ramdas and Tibshirani 2023 (2202.13415v5), its jackknife+ appendix: the definition and quantile conventions, the background theorem, the proof via Landau's theorem, and the cross-conformal remark.
- Romano, Patterson and Candès 2019 (1905.03222v1), Theorem 2 and their quantile convention.
- Barber et al. 2021 is not on disk. It is read through the 2023 restatement, and nothing was fetched.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean, no warnings.
- `lake build Enthymemata`: green; the only warnings come from apologies already present in other files.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- Builds ran at below-normal priority, and nothing was fetched. D: had 39 GB free at the end.

**Follow-ups, each with the decision it would change.**
- *CV+ with K folds.* It decides whether the layer offers K fits instead of n.
  - It needs Barber et al. 2021 for the added term; fetching the paper is your call.
  - It also needs a finite statement of its own, because the test unit belongs to no fold.
- *Open, no decision waiting:*
  - whether the jackknife+ alone needs a symmetric fit, since the search found no counterexample;
  - the law forms on a continuous data space, as before.

**What stopped me:** nothing.

### 2026-11-08 — E35: M-1 and M-2 proved

**Proved** in `lean/Enthymemata/ConformalValidity.lean`, which is now 119 declarations and still PROOF-CLOSED. `#print axioms` on the E35 results shows the standard three only. [conformal-lemmas.md](conformal-lemmas.md) has the record, in its E35 section.
- *M-1:*
  - `seqRank_uniform`, the finite core: the randomized ranks of the smoothed p-values are independent and uniform over every ordering of a fixed data set, ties included, for any scores symmetric at each time;
  - `sum_tieRankProb`: exactly one tie class covers each rank;
  - `seqRank_law`: the same under any exchangeable law.
- *M-2:*
  - `lintegral_capital_snoc` and `lintegral_capital`, the link: given any history the next capital integrates to the current one, so E Sₙ = 1;
  - `lintegral_capital_le`, its supermartingale form;
  - `lintegral_stoppedCapital`: the capital at a bounded stopping time is an e-value;
  - `ville_horizon` and `ville`: P(∃ n, Sₙ ≥ c) ≤ 1/c.

**Not formalized, checked exactly instead.**
- M-1's lift from ranks to continuous p-values: θT = (J − 1) + V. 11,400 exact joint CDF values equal x₁⋯x_N.
- The e-value at an unbounded stopping time, by Fatou.
- CUSUM and Shiryaev–Roberts: posed, with their route. E Rₙ = n follows from `lintegral_capital`; the rest needs optional stopping at unbounded times.

**Checks first** (`conformal-checks.py`, section E35). All assertions hold.
- *M-1:*
  - the ranks are exactly uniform in 68 cases, inductive and full scores;
  - each premise broken in turn fails: ties counted for the unit; stale full-form scores; tests sharing one calibration set, where P(both lowest) = 1/10, not 1/16.
- *M-2, exact on a 4-point grid,* for Simple Jumper and two power strategies:
  - the martingale property holds on every history;
  - E S₄ = 1, and c·P(max S ≥ c) ≤ 1;
  - a bet integrating to 11/10 breaks it, and so does one reading the current p.
- *Seeded Monte Carlo:*
  - going all in on p₁ ≤ 1/c attains 1/c;
  - small bets approach it (0.0995 against 0.1);
  - Simple Jumper and the power martingale stay under it.

**What changed from the task.**
- *Ville is not proved by Mathlib's maximal inequality.* The bet is switched off at the alarm: the stopped capital is again a betting product with expectation at most 1, and Markov's inequality gives the bound at each horizon. The limit over horizons is a monotone union. The route needs Tonelli and Markov only, and states the code's premises as they are. Mathlib's route needs the product as a conditional-expectation martingale, which needs a freezing lemma, and modules that were not cached. Nothing was fetched.
- *Ville's premise is weaker than the task's.* Each bet need only integrate to at most 1, and the product is then a supermartingale. The martingale link and the e-value's equality keep the exact premise.

**Gates for E35's code**, each a premise of a statement above:
1. The whole monitored sequence is exchangeable. In the inductive form, it is exchangeable given the training set proper, which fits the score and is not ranked.
2. Each new score is ranked among every score so far, itself included, calibration and test alike, never against a fixed calibration set.
3. Ties are broken by a fresh θₙ, uniform and independent of the data and of earlier draws: pₙ = (L + θₙT)/n. Neither CB-1's tie rule nor the arrival order.
4. The scores at time n are a symmetric function of the n units so far: a fixed A, or all recomputed from the bag. Never scores kept from earlier times.
5. Each betting function is ≥ 0 and measurable, reads only p₁, …, p_{n−1}, and integrates to at most 1 on [0, 1]. It integrates to exactly 1 where the e-value's equality is wanted; a quadrature-normalized plug-in rounds down. Simple Jumper and the power mixture qualify, because their step ratios are predictable mixtures of betting functions.
6. The capital starts at 1. An alarm or a stop rule reads only the capital and p-values so far.

**Corrections to E35's entry**, written into its Lean bullet:
- Ville needs ≤ 1, not = 1.
- "The full form, each score recomputed" means every score, at every time.
- The finite argument covers the ranks, and the lift is measure-theoretic.
- Ville is proved without Mathlib's martingale library.
- CUSUM and Shiryaev–Roberts are posed, not proved.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- D: had 40 GB free, and nothing was fetched.

**Open, no decision waiting:**
- formalizing the lift and the unbounded stopping time;
- the CUSUM and Shiryaev–Roberts proofs.

The CUSUM and Shiryaev–Roberts guarantees are cited from Vovk 2021, where they are proved.

### 2026-11-06 — CB-3 with an error bound per unit

The follow-up below on per-unit bounds is done, at your call. CB-3 now takes a bound εⱼ for each unit. With every εⱼ = ε it reads as before (2ε, 4ε), so no earlier claim changes.

**Proved.** `ConformalValidity.lean` has 89 declarations and remains PROOF-CLOSED. `#print axioms` on the CB-3 results shows the standard three only.
- `accepts_or_fragile`, `acceptsInflated_of_accepts`, `accepts_shift_of_accepts_close` and `fragileEstimated_of_fragile` are restated with per-unit bounds.
- `coverage_close` now takes bounds carried by the units. `coverage_inflated` takes bounds that may differ by ordering and by position.
- `card_fragile_le_estimated` reads the band with each ordering's bounds.
- New: `fragileEstimated_mono`, so any dominating bounds serve; and `ratio_permsOn`, the counting step that `coverage_eq` and `coverage_close` share.

**Checks first.** `conformal-checks.py` runs the per-unit forms on 150,435 deterministic cases, exhaustive with per-unit bounds to m = 3, plus the adversary's attainment, over orderings with bounds by unit and by position, and in the conformal-Bayes case. Every assertion holds.

**Finding: where per-unit bounds pay.** Monte Carlo error attaches to a position, the candidate's or a data unit's, more than to a unit.
- The inflated rule takes bounds by position, so it gains directly.
- The margin needs bounds carried by the units, and each unit's worst error over the orderings includes its turn as the candidate. In the conformal-Bayes check that equals the uniform bound.
- So the per-unit margin pays only where units differ in themselves.
- This sharpens the decision. The layer inflates with per-position bounds, the candidate's error and each data unit's apart. It reports the margin with unit-carried bounds, and the uniform bound is the safe choice when no unit stands apart.

**Gates.**
- `lake build Enthymemata.ConformalValidity`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- D: had 41 GB free.

### 2026-11-06 — CB-3's prior art, read in full

The four arXiv papers are now in `gauntlet/conformal-bayes`, and Ndiaye's is in `PDFs/Bayesian/ndiaye22a.pdf` (ICML 2022, PMLR 162). Each was read at its theorems. [conformal-lemmas.md](conformal-lemmas.md) has the detail under CB-3's prior art.
- *No statement changes.* The checks, the Lean and the earlier report's findings stand.
- *The inflated rule is known.* It is Ndiaye's stable set (Prop. 3.2 to 3.4) with a uniform ε in place of per-unit sure bounds τᵢ. The RKHS paper gives the same construction with τᵢ(y). `coverage_inflated` is the finite, sure-bound case of both, and `accepts_shift_of_accepts_close` is Ndiaye's deflated side.
- *The margin was not found.* The fragile count bounding the uninflated rule (`coverage_close`), and its attainment, are in none of the five. Bhattacharyya, Zhang and Barber call uncorrected approximations unguaranteed in general. Abad Martinez et al. leave their scores' finite-sample validity open.
- *The probabilistic-ε follow-up is half answered.* For the inflated rule, Bhattacharyya et al.'s appendix gives 1 − α − 3√ν when each unit's error is within ε with probability 1 − ν, test unit included, and their Example 3 is Fong and Holmes's add-one-in score. The margin's form under such an ε stays open.
- *A new option for the fallback decision: the tournament correction.* It uses full data and covers at least 1 − 2α with no assumption on the approximation, randomized scores included. For conformal Bayes it needs draws given each leave-one-out set (n chains), and those chains need to treat their data symmetrically but not to mix. Their rejection shortcut from the n-unit draws needs exact draws and a likelihood bounded below. Barber, Candès, Ramdas and Tibshirani's swap randomization would need a refit per candidate.

**Follow-ups, each with the decision it would change.**
- *Per-unit bounds.* Ndiaye's τᵢ and the RKHS paper's τᵢ(y) suggest stating CB-3 with a per-unit εⱼ, the layer's per-unit MCSE in place of the largest one. That decides how wide the inflated set is and how large the reported margin is.
- *E16's sampler bullet* ("What the conformal literature adds") weighs the full step on good draws against the split form. The tournament form is a third choice: full data at 1 − 2α, without mixing. That decides the fallback when the chain fails its diagnostics. At your call it was added to the bullet on 2026-11-06.

### 2026-11-06 — CB-1 to CB-4 posed, checked and proved

**Proved** in `lean/Enthymemata/ConformalValidity.lean`, namespace `Conformal`, added to `Enthymemata.lean`. There are 88 declarations, no apologies, and the file is PROOF-CLOSED. `#print axioms` on the 28 main results shows `propext`, `Classical.choice` and `Quot.sound` only. The statements, the checks' figures and the decisions are in [conformal-lemmas.md](conformal-lemmas.md), and the exact-arithmetic checks in [conformal-checks.py](conformal-checks.py), run before any Lean.
- CB-1: `card_rankLe_le` (the counting core), `coverage_orderings`, `coverage_orderings_eq`, `coverage_strict_le`, `coverage_orderings_le`, `card_rank_eq`, `rank_uniform`, `coverage_tieBreak`, `coverage_given`, `coverage_law`, `coverage_stratum`.
- CB-2: `split_coverage`, `splitScore_equivariant`, `exchangeable_given_training`, `cqr_ordConnected`.
- CB-4: `sum_jointU`, `postU_addOneIn`, `post_addOneIn`, `integratedScore_equivariant`, `latentScore_equivariant`, `latentScore_eq`, `sum_jointU_addOneIn`, `integratedScore_le_latentScore`, `accepts_latent_of_mixed`, `mixed_counterexample`.
- CB-3: `accepts_or_fragile`, `accepts_shift_of_accepts_close`, `acceptsInflated_of_accepts`, `fragileEstimated_of_fragile`, `coverage_close`, `coverage_inflated`, `card_fragile_le_estimated`.

**Open.**
- The law forms on a continuous data space. They follow by conditioning on the multiset, but are not formalized.
- CB-3 when ε holds only with high probability.
- The k-th order statistic moving at most ε. It is checked, not formalized; the count form replaces it.

**What changed, and what proved false.** Nothing proved false.
- *CB-1.* The brief's "at most 1 − k/(n + 1) + 1/(n + 1)" is Lei et al.'s bound for a level off the grid. At a grid level, with distinct or randomly tie-broken scores, coverage is exactly 1 − k/m. With ties, the liberal rule's excess is the share of units whose ties straddle the threshold, and it can reach 1. Hoff's strict rule, which he states without ties, covers at most 1 − k/m and refuses every candidate when all scores tie. The per-stratum form is stated for permutations inside the test unit's stratum (Fong and Holmes's form, where the score may read other strata in any way), with `coverage_stratum` as the stratum-map corollary.
- *CB-3.* The band is one-sided and tight. A decision is lost only for accepted units within 2ε above the highest refused score, M = {u : s₍ₖ₊₁₎ ≤ s_u < s₍ₖ₎ + 2ε} with distinct scores, and an adversary within ε attains the bound for every unit. The estimated band is ŝ₍ₖ₊₁₎ − 2ε ≤ ŝ_u < ŝ₍ₖ₎ + 4ε, the brief's "within 4ε" on the high side. A third policy came out of the same lemma: inflating the threshold by 2ε keeps 1 − k/m exactly.
- *CB-4.* Posed as in the brief, with the direction of (d) proved. The mixed score accepts no candidate the latent score refuses, so it undercovers and never over-covers. It covers 0 against 1/2 in the Lean counterexample, and 0.51 to 0.53 against 0.75 in the checks' model.

**Lines E16 can take.**
- *CB-1:* the accept rule counts ties for the candidate, π(y) > k/m. It covers at least 1 − k/m whatever the ties, and exactly 1 − k/m with distinct or tie-broken scores; the strict form is never used. A stratum of m units attains only 1 − k/m, so below m = ⌈1/α⌉ the set is the whole space and the layer says so. Ranks are taken within the stratum.
- *CB-2:* the split form holds for any fit, the sampler's included, and given the training split (`split_coverage`, `coverage_given`). It is the guarantee that does not lean on the sampler. CQR's set is an interval.
- *CB-4:* every unit, the candidate included, is read by one density. The default is latents integrated (OQ17's scope rule, GBASS's pointwise content now); the alternative is the latent-posterior density, with the candidate's latent drawn per draw. Data units at their retained latents with the candidate integrated undercover, so the layer refuses a latent-conditioned matrix.
- *CB-3:* an estimate within ε costs at most |M|/m, and |M| is read from the estimates (`card_fragile_le_estimated`). The layer reports the margin, inflates by 2ε, or falls back to the split form. ε must be a sure bound, or its failure probability adds to the miscoverage.

E16's "What the conformal literature adds" is unchanged: no result corrects it. One result corrects a line outside it, and at your call that line was corrected on 2026-11-06. The bullet "Why the posterior predictive is the score" said Hoff's score "is the add-one-in step above". Hoff's cᵢ is each unit's predictive given every other unit, a leave-one-out. Fong and Holmes's step scores in-sample, p(zᵢ | all units). The two rank differently: on 12 of 216 data sets in the checks' model, and in Hoff's own setting with no covariate, y = (2, 2, 2, 0) puts the lone 0 lowest in-sample and highest by Hoff's. Both cover; Hoff's optimality attaches to cᵢ only.

**Prior art.**
- Read: Fong and Holmes's §2–3 (the in-sample score's invariance, the group proposition, the remark on MC samples and exchangeability); Hoff 2023, §4.1.
- Cited through Fong and Holmes, not read: Vovk, Gammerman and Shafer (2005); Lei et al. (2018).
- CQR from its abstract.
- For CB-3, from abstracts found by web search (none on disk): Barber, Candès, Ramdas and Tibshirani (arXiv 2202.13415), with total-variation gap terms and randomization for asymmetric algorithms; Ndiaye (ICML 2022, arXiv 2112.10224), a stability bound added to keep coverage, which is the inflated rule; Razafindrakoto, Celisse and Lacaille (arXiv 2601.13102), per-point slack whose region contains the full one; Bhattacharyya, Zhang and Barber (arXiv 2605.29200), the tournament correction; Abad Martinez et al. (arXiv 2202.01315), influence functions.
- None of those abstracts states the one-sided fragile count or its attainment.

**Self-review passes and fixes.**
- *CB-1.*
  - The upper bound was sharpened to an equality.
  - The strict rule's failure under ties and the straddle band were added.
  - The per-stratum hypothesis was weakened to the test unit's stratum.
  - Checks: the first position-reading score favoured the test and showed no violation, so it was replaced by one that marks the test down. The first across-strata contrast did not undercover, so it was replaced by the partially pooled case (0.472 against 0.6). The within-stratum equivariance check covered 20 sorted multisets and now covers all 243 tuples.
- *CB-2.* The guarantee was marginal only; `coverage_given` now states it given the training split.
- *CB-4.* "Breaks symmetry" was sharpened to the proved direction, undercovering. Reviewing E16's claim that the score is Hoff's turned up the Hoff finding.
- *CB-3.*
  - `rankShift_le_rankLe` first lacked u ∈ C, which left a unit outside C with no closeness bound. The hypothesis was added.
  - The draft band's threshold was s₍ₖ₊₁₎ + 2ε. The tightness check fixed it at s₍ₖ₎ + 2ε, the candidate's self-count included.
- *Lint:* long lines, unused decidability instances and `show` were fixed under `lake build`.

**Build gates.**
- `lean` on the file, iterating.
- `lake build Enthymemata.ConformalValidity`: clean, no warnings.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- The full `lake build`: 3,180 jobs, green. A `--no-build` dry run first showed only the barrel stale. A final `lake build Enthymemata` ran after the docstring edit.
- Everything ran at below-normal priority. Every import was already cached, so nothing was fetched.
- D: had 41 GB free at the start and the end.

**Follow-ups, each with the decision it would change.**
- *Hoff's score in the layer.* This decides whether the layer offers Hoff's leave-one-out cᵢ, which needs a weight of g(z_new | θ)/g(zᵢ | θ) for each data unit and brings PSIS's instability, E16's Hoff bullet was corrected on 2026-11-06 and no longer claims Bayes-optimal width for the add-one-in step.
- *ε as a probabilistic bound.* This decides how the layer sets ε (a concentration bound on the self-normalized add-one-in estimate, or a multiple of its MCSE) and what miscoverage it budgets for ε failing. It changes the margin's reported form.
- *Open, no decision waiting:* whether the one-sided fragile count is in the literature; the continuous-space law forms; promotion to `Lemmas/`, which waits for semantic review.

**What stopped me:** nothing.
