# Conformal validity in the Lean harness — chip brief

Date: 2026-11-05. Status: open. Written at the user's request, from E16's reading of the conformal literature.

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

## Reports

Newest first.
