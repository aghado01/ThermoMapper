# Conformal validity in the Lean harness — chip brief

CB-1 to CB-4 were proved on 2026-11-06 and E35's M-1 and M-2 on 2026-11-08. Their reports are below. The next assignment is E16 slice 9.

## E16 slice 9 — cross-fitted intervals and asymmetric CQR (2026-11-09)

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
