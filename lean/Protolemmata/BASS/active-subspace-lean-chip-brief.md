# Active-subspace geometry in the Lean harness — chip brief

Commissioned by the user on 2026-11-13 through the coordinating session. The four targets were proposed on 2026-11-10 ([handoff](../../../issues/regression/coordination-handoff-2026-11-11.md#waiting-on-the-user), item 11) and exist only as names. Posing them as statements is the first job.

## Context

E8 reads a posterior active subspace from any fit ([active-subspace-brief.md](../../../issues/regression/bass/active-subspace-brief.md)). Its slice 5 readouts are built on formulas the brief states, two of them cited from memory. Its slice 7 refits in the subspace (rotate-and-refit). E36, BPPR, learns the directions while it fits ([bppr-brief.md](../../../issues/regression/bass/bppr-brief.md)); its brief is a draft waiting for the user, and the coordinating session recommended that its moves' balance go Lean first.

A statement that turns out false, or true only with a further hypothesis, is a finding worth more than a proof.

**Read first:**
- [lean/README.md](../../README.md): the taxonomy, meta-CI's rules, import discipline.
- [guarded-estimator-lemmas.md](guarded-estimator-lemmas.md) and [structure-chain-lemmas.md](structure-chain-lemmas.md) with their modules in `lean/Enthymemata/`, for the style, the finite sums and the reversible-jump balance already proved.
- E8's brief: slice 5 (the settled Grassmann items), slice 7's plan and the reports of 7a to 7d.
- The BPPR brief: decision 5 and slice 0's row.
- The [Grassmann discussion](../../../issues/regression/discussion/grassmann-stiefel-subspaces-brief.md), for the derivations slice 5 cites.

## Targets

For each: the statement, its hypotheses, the decision it changes and its prior art, in the protolemma. **Numeric checks come before any Lean:** exact rational arithmetic where it fits, small cases enumerated, seeded random cases beside chosen ones. The order puts first what already-built code rests on.

1. **P̄ over mixed ranks** (E8 slice 5). With P_s orthogonal projections of ranks k_s and weights w_s ≥ 0 summing to 1:
   - P̄ = Σ w_s P_s is symmetric with 0 ≼ P̄ ≼ I and tr P̄ = Σ w_s k_s.
   - For a projection P of any rank, Σ_s w_s ‖P − P_s‖²_F = tr P + Σ w_s k_s − 2 tr(P P̄). Over projections of every rank it is least at the P that keeps P̄'s eigenvectors with eigenvalue above ½, with an eigenvalue of exactly ½ indifferent.
   - Across ranks, ½‖P_s − P‖²_F = |k_s − k|/2 + Σ sin²θ_i over the min(k_s, k) principal angles.
   - *Prior art:* the brief cites Ye and Lim 2016 for the distance and Barbieri and Berger 2004 for the median probability reading, both from memory. Say what each actually states, read or from abstracts.
   - *Decision:* whether the readout's rank summary and its cross-rank distance stand as documented.
2. **The matrix Poincaré check** (E8 slice 5; its decision 6 names the marginals).
   - *The reduction:* Cov_ρ(f) ≼ E[J D Jᵀ] holds exactly when, for every v, Var(vᵀf) ≤ Σ_i D_i E[(∂_i vᵀf)²].
   - *The product law's constant:* Efron–Stein gives Var(g) ≤ Σ_i E[Var_i g]. Each coordinate's one-variable inequality then puts its own constant on the diagonal D, coordinate by coordinate, where a scalar bound would take the largest.
   - *The uniform's constant,* (b − a)²/π², and so Cov_ρ(f) ≼ G/π² under the box metric, and a largest generalized eigenvalue of at most 1.
   - *Form:* Efron–Stein over finite product spaces is the finite core. The one-variable inequality is a hypothesis unless Mathlib carries it; say which.
   - *Decision:* whether a reading above 1 beyond Monte Carlo error marks a defect for every marginal decision 6 covers, or only for some.
3. **E36's direction moves** (the BPPR brief, decision 5).
   - *The change move.* The power spherical density is a function of μᵀx alone, so on one sphere q(θ′ | θ) = q(θ | θ′) and the acceptance ratio carries no proposal term. State it on the active coordinates' sphere, with the active set and J kept.
   - *Birth and death.* A direction drawn from its prior on the active coordinates' sphere, with the number of active coordinates varying between ridge functions: the reversible-jump ratio and its dimension matching. Reuse the structure chain's balance where it fits.
   - *Decision:* the acceptance ratios BPPR's first build slice would code. Give the BPPR brief a line it can take.
4. **Rotate-and-refit's non-invariance** (E8 slice 7). No record states it beyond the name. Pose both readings below, say which the record supports, and prove what holds.
   - *The frame.* For d ≥ 2, the tensor-product hinge class in u = W₁ᵀx is not closed under an orthogonal change of u, so the refit depends on the frame chosen within the active span. Where the active block's eigenvalues are near-equal that choice is arbitrary (slice 7's decision 2). A counterexample is the natural form.
   - *The scale.* C_z = LᵀC_xL, so the subspace is equivariant under orthogonal changes of the inputs and moves under any other. For the projected fit: the g-form's marginal is invariant to a column's scale; the column-mass support and the raw ridge are not; and under the training-range rule the evaluated hinge is invariant to the declared box while κ is not (slice 7a's decision 4 gate and slice 7d).
   - *Decision:* for the frame, whether the eigengap report within the active block is the right warning. For the scale, whether the training range as BASS's default (the handoff's item 2) rests on a proved invariance.

## Added targets: a penalized system's admission (2026-11-14)

Commissioned by the user on 2026-11-14 through the coordinating session. They sit outside active subspaces and come to this chip for its harness. The environment, working rules and reporting below hold unchanged.

**Context.** BAPS factors A(c) = Σₖ cₖSₖ at every sampler step: ZᵀZ + λP, ZᵀZ + Σⱼ λⱼPⱼ, and ZᵀZ + Σᵢ λᵢdᵢdᵢᵀ for the locally adaptive penalty. Every Sₖ is positive semidefinite and fixed by the design, and every cₖ > 0. A gate admits a factor that is positive definite with an estimated κ₁ of the unit-diagonal scaling E(c) at most 10¹⁴ (`BapsAdmission`).
- [Slice 6 of the defects brief](../../../issues/regression/shared-code-defects-brief.md#slice-6-the-samplers-condition-check-as-a-setting) builds a screen on a bound relative to one reference point c_ref. The defects chip (`local_0e31b37f`) is building it now and derives the screen's safety factor S.
- The [BAPS follow-up's admission input](../../../issues/regression/bass/extensions-brief.md#input-for-the-follow-up-admission-and-singular-systems-2026-11-14) states the two failures behind the gate.
- Every statement below was derived by the coordinating session on 2026-11-14 and checked against no run. Target 7's proof sketch is its own, not read from van der Sluis's paper.

**Hypotheses an instance must meet.** Each result carries a concrete instance in the checks and in the report: a family with at least one positive definite combination, strictly positive coefficients, and a nonempty index set for the max and the min. S₁ = diag(1, 0), S₂ = diag(0, 1) is the smallest.

5. **Identification does not read the coefficients.** For c > 0 entrywise:
   - xᵀA(c)x = 0 exactly when Sₖx = 0 for every k, so ker A(c) = ⋂ₖ ker Sₖ at every such c;
   - A(c) is positive definite at one c > 0 exactly when it is at every c > 0, and exactly when the Sₖ share no null vector but 0;
   - BAPS's instance: ZᵀZ + Σⱼ λⱼPⱼ is singular exactly when some x ≠ 0 has Zx = 0 and Pⱼx = 0 for every j.
   - *Decision:* whether one check of the design at construction stands on a proved footing (the follow-up's direction 1). It also says what a per-step positive-definiteness check can still catch: a loss in the computed factor only.
6. **The bound.** With ρₖ = cₖ / c′ₖ, m = minₖ ρₖ and M = maxₖ ρₖ:
   - m · xᵀA(c′)x ≤ xᵀA(c)x ≤ M · xᵀA(c′)x for every x;
   - λ_min(A(c)) ≥ m · λ_min(A(c′)) and λ_max(A(c)) ≤ M · λ_max(A(c′));
   - where A(c′) is positive definite, κ₂(A(c)) ≤ (M / m) · κ₂(A(c′));
   - M / m is attained (the diagonal instance above, with c′ = (1, 1)), so no smaller constant serves every family;
   - a coefficient held at 1 in both, as ZᵀZ's is, enters with ratio 1.
   - *Form:* the quadratic-form line is the core. Say whether the eigenvalue lines are proved through Mathlib's ordered eigenvalues or through Rayleigh bounds.
   - *Decision:* slice 6's screen rests on it. Sharpness says the window of c the solve carries widens only with another form of the solve (the follow-up's direction 2), never with a better inequality.
7. **From the bound to what the gate reads.** Restated on 2026-11-14 from the defects chip's derivation (`76851ee`, its In-progress entry). That chain is tighter than the first statement here and needs no scaling theorem.
   - *The diagonal.* Every (Sₖ)_ii ≥ 0, so A(c)_ii lies between m and M times A(c′)_ii.
   - *The scaled matrix's smallest eigenvalue.* λ_min(E(c)) ≥ (m / M) · λ_min(E(c′)).
   - *Norms.* ‖E(c)‖₁ ≤ q, from a unit diagonal, entries of at most 1 in absolute value and at most q nonzeros in a column. For any n × n matrix ‖A‖₁ ≤ √n · ‖A‖₂, and for a symmetric one ‖A‖₂ ≤ ‖A‖₁.
   - *Composed.* κ₁(E(c)) ≤ q · √n · (M / m) · ‖E(c′)⁻¹‖₁.
   - *Outside these hypotheses:* one hierarchical coefficient, ξ², can reach 0. The screen hands that step to the estimate.
   - *Outside Lean:* rounding. The screen multiplies the bound by 2 for it, and the defects chip's sweep tests that. The per-step estimate is a lower bound on κ₁ by its construction, which is the direction the screen needs. The reference's ‖E(c′)⁻¹‖₁ is read exactly from n solves, so no estimator's shortfall enters.
   - *Off the path:* the Jacobi-scaling bound κ₂(B) ≤ q · κ₂(FBF) for a unit-diagonal B (van der Sluis 1969, Numer. Math. 14:14–23) was the first statement's route. No decision rests on it now, so it is recorded here and not queued.
   - *Decision:* whether the screen admits only what the estimate admits by proof in exact arithmetic, with rounding alone left to the sweep. Give the defects brief a line it can take, and tell the defects chip by message when the chain is settled either way.

**Order.** 5, then 6, then 7. Numeric checks before any Lean, as for the first four: exact rational arithmetic for the quadratic-form lines, seeded random families with rank-deficient Sₖ with and without a common null vector, and small banded cases built from a B-spline Gram and difference penalties.

**Files.** `lean/Protolemmata/BAPS/penalized-admission-lemmas.md` and `penalized-admission-checks.py`; the module under `lean/Enthymemata/`, named and placed by the README's taxonomy. Outside it: `src/`, `tests/`, `docs/`, the defects brief and the extensions brief, each of which gets a line it can take in the report.

## Scope

New files:
- `lean/Protolemmata/BASS/active-subspace-lemmas.md`: the statements, with a status line per target;
- `lean/Protolemmata/BASS/active-subspace-checks.py`: the numeric checks;
- the Lean modules under `lean/Enthymemata/`, named and placed as the README's taxonomy says, with their lines in `Enthymemata.lean`. Say in the report where each target landed.

Outside it: `src/`, `tests/` and `docs/`; promotion to `Lemmas/`, which waits for semantic review; and edits to E8's or BPPR's brief, each of which gets a line it can take in the report.

## Environment

- **Toolchain:** as pinned in `lean/lean-toolchain` and `lake-manifest.json`, unchanged.
- **Mathlib's build output may be absent** from `lean/.lake`.
  - Fetch by scoped module, `lake exe cache get <Module>`.
  - Check `df -h /d` before and after each fetch. If a fetch would leave under 3 GB free, stop and report: disk decisions are the user's.
  - Build each new module with `lake build <Module>` and run `scripts/meta-ci.ps1 -NoBuild`. Run the full `lake build` only if the disk allows, and say which ran.
- **Iterate fast:** `lean` directly, with `LEAN_PATH` from `lake env printenv LEAN_PATH`. Linters run only under `lake build <Module>`, so finish with that.
- **Quirks:** `Π` is a reserved token; identifiers with a combining macron fail to parse; `omit … in` and `open … in` go before the docstring; `lean/` is LF, so write files from Python with `newline='\n'`; meta-CI's ledger counts the word for an apology in comments too.
- **Python:** `C:/Users/azrie/PDenv/PyenvPython/pyenv-win/versions/3.12.6/python.exe`, run with `-I -X utf8` on a script file.
- **Papers:** `D:/aghado01/codex-scientiae/supellex-private/gauntlet/<folder>/<id>/`, each folder with an `inventory.jsonl`; the power spherical law is `bayesian/2006.04437v2`. PDFs are under `…/supellex-private/PDFs/`. `pdftotext` works, and `pypdf` keeps math glyphs.
- **Active modules import declaration-owning Mathlib modules,** never `import Mathlib` ([AGENTS.md](../../../AGENTS.md#environment)).

## Working rules

- **Git.** On `main`, committed by pathspec, nothing left staged, nothing pushed, no attribution lines.
- **The machine and the checkout are shared.** Builds run at below-normal priority. Three chips are at work in `src/` and `tests/`. The Conformal Lean chip (`local_21649695`) is at work in `lean/Enthymemata/ConformalValidity.lean` and `lean/Protolemmata/Conformal/`: this chip owns its own new files, tells that chip before editing `Enthymemata.lean`, and agrees with it by message when both need `lake build` in the one workspace.
- **Self-review.** After each target, a fresh-eyes pass before the next: does the Lean say what the protolemma claims, with no hypothesis doing the work? Run `#print axioms` on the main results. Name the passes and their fixes in the report.
- **Verification.** The coordinating session checks each landed target from a clean archive ([scripts/verify-lean.ps1](../../../scripts/verify-lean.ps1)), so the report names the modules and the declarations to print axioms for.
- **Sub-agents** need the user's permission.
- **Events.** Commit the protolemma and its checks for all four targets before the first proof, and end that turn with one line and its link. After that, end a turn with one line and the report's link when a target lands or needs a decision.

## Reporting

Append a report under "Reports", newest first, and keep the protolemma's status lines current. It says:
- what is proved, by declaration name, and what remains open;
- every statement that changed, why, and any that proved false;
- the line each of E8's and BPPR's briefs can take;
- prior art for each statement, read or from abstracts;
- which build gates ran, and the disk's free space at the end;
- follow-ups with the decision each would change; a question with no decision waiting is recorded as open.

Hand-over state between turns is an "In progress" entry here.

## Reports

### 2026-11-14 — the disproof verified from a clean archive, and rebuilt apart (coordinating session)

At the user's request, ahead of the chip's write-up.

- **Commit:** `fc793df`, through `scripts/verify-lean.ps1` with `-Module Enthymemata.PenalizedAdmission` and the BAPS checks.
- **Result:** 26 of 26 checks; 2,634 jobs with no warnings and no errors; meta-ci PASS.
- **Axioms:** every theorem of the module, 89 in all, prints `propext`, `Classical.choice` and `Quot.sound` only, `inverse_part_unbounded` among them.
- **The statement, read.** For every C there is a two-term positive semidefinite family with the bracket at m = 1 and M = 2 and A(c′) positive definite, whose ‖E(c)⁻¹‖₁ exceeds C·(M/m)·‖E(c′)⁻¹‖₁. These are `cond1_unitScale_le_sharp`'s own hypotheses, and the spike family meets them.
- **Rebuilt apart from the chip's proof and script,** in exact rationals: v = (3/5, 4/(5r), …) on r² + 1 coordinates, A(c) = 2I − vvᵀ with inverse (I + vvᵀ)/2 by Sherman–Morrison, checked entry by entry.
  - ‖E(c)⁻¹‖₁ reads 1.47, 3.27, 5.45 and 18.50 at r = 1, 5, 10 and 40.
  - Each lies above the proved floor 6r/25 and below the proved ceiling (M/m)(1 + √n)/2: 2.41, 6.10, 11.05 and 41.01.
  - κ₁(E(c)) / ((M/m)·κ₁(E(c′))) grows as well: 0.97, 3.8, 10 and 107. So no constant serves the whole κ₁ bound over every family either.
- **What it does not say.** The spike family is dense: every column of vvᵀ is full. A banded family with q fixed is the open case, and the chip's search is on it.
- **Receipts:** `artifacts/verify-lean-fc793df/`, with the axioms file, the summary and the rebuild's script and output tracked in [receipts](../../../issues/regression/receipts/verify-lean-fc793df/).

### 2026-11-14 — the two follow-ups verified from a clean archive (coordinating session)

- **Commit:** `dd57650` (the Lean is `7f535b3`), through `scripts/verify-lean.ps1` with `-Module Enthymemata.PenalizedAdmission` and the BAPS checks.
- **Result:** 23 of 23 checks; 2,634 jobs with no warnings and no errors; meta-ci PASS.
- **Axioms:** every theorem of the module, 78 in all, prints `propext`, `Classical.choice` and `Quot.sound` only. The nine new names and the reproved `cond1_unitScale_le` are among them.
- **Statements, read for an instance that meets them:**
  - *The bridge,* `lamMin_eq_iInf_eigenvalues` and `lamMax_eq_iSup_eigenvalues`, asks only that A be symmetric.
  - *The column step,* `norm1_le_half_one_add_sqrt`: a symmetric positive semidefinite P with form at most μ‖z‖² has ‖P‖₁ ≤ μ(1 + √n)/2. Checked by hand on rank-one projectors, whose largest column sum is exactly (1 + √n)/2.
  - *Its witness,* `norm1_columnWitness`: vvᵀ/12 with v = (3, 1, 1, 1) meets every hypothesis and attains the constant at n = 4.
  - *The composed chain,* `cond1_unitScale_le_sharp`, keeps the earlier hypotheses. `cond1_unitScale_le` keeps its statement and is its corollary.
- **For slice 6:** the screen may take (1 + √ν)/2 in place of √ν. The defects chip holds that change until the user answers the default question.
- **Receipts:** `artifacts/verify-lean-dd57650/`, with the axioms file and the summary tracked in [receipts](../../../issues/regression/receipts/verify-lean-dd57650/).

### 2026-11-14 — the penalized admission's two follow-ups: the eigenvalue bridge and the sharp column step

At the user's request ("do them anyway"), after the targets were verified (`af5cd65`).

**Landed:** `7f535b3`, the module `lean/Enthymemata/PenalizedAdmission.lean`. The checks gained six items and pass 23 of 23. The new [penalized-admission-searches.py](../BAPS/penalized-admission-searches.py) gives the readings.

**What is proved.**
- *The bridge.* `lamMin_eq_iInf_eigenvalues` and `lamMax_eq_iSup_eigenvalues`: for a symmetric A, the Rayleigh extremes are Mathlib's least and greatest `Matrix.IsHermitian.eigenvalues`. So `cond2` is κ₂, and the first departure from the wording is closed.
  - The proof goes through `Matrix.IsHermitian.spectral_theorem`: xᵀAx = Σᵢ λᵢ(Uᵀx)ᵢ² with (Uᵀx)ᵀ(Uᵀx) = xᵀx (`form_eq_sum_eigenvalues`, `eigenvector_facts`).
  - The module now imports `Mathlib.Analysis.Matrix.Spectrum` (2,634 jobs). The earlier results still use no spectral theorem.
- *The sharp column step.* `norm1_le_half_one_add_sqrt`: a symmetric positive semidefinite P with form at most μ‖·‖² has ‖P‖₁ ≤ μ(1 + √n)/2. `norm1_columnWitness` attains it at n = 4 with P = vvᵀ/12, v = (3, 1, 1, 1).
- *The composed bound with it.* `cond1_unitScale_le_sharp`: κ₁(E(c)) ≤ q·((1 + √n)/2)·(M/m)·‖E(c′)⁻¹‖₁. `cond1_unitScale_le` keeps its statement and is now this theorem's corollary.

**What the searches read** (the protolemma's findings 4 and 5):
- (1 + √n)/2 is sharp for the column step, even for a unit-diagonal E (n ≥ 4).
- q is reached by c₁J + c₂I.
- The inverse part's reading stays near 1 up to n = 16, at most 1.23, while its proved bound grows to 2.50. The search runs over every family at once: every family gives m·A(c′) ≼ A(c) ≼ M·A(c′), and every such pair is a two-term family.
- Conjecture, open: κ₁(E(c)) ≤ C·q·(M/m)·‖E(c′)⁻¹‖₁ with C a constant. Each step alone is sharp, so a proof needs the coupling between E(c) and E(c′).

**Self-review.**
- *Pass 1, statements against the claims.* The bridge's two equalities are stated with Mathlib's own `eigenvalues`, over `⨅`/`⨆` on a nonempty finite index.
- *Pass 2, can the hypotheses be met.* The column witness meets every hypothesis of `norm1_le_half_one_add_sqrt` and attains its bound, so the theorem is not vacuous and its constant is exact.
- *Pass 3, proof repairs.* The first draft's key inequality failed under `nlinarith`. It now goes through an explicit identity, (wμ/2)² − (x + sr)² = (sx − r)², by `linear_combination`. A check's tolerance was set below rounding at condition 10⁹ and was loosened to 10⁻⁶. The readings themselves were unchanged.

**Declarations to print axioms for** (10 new): `form_inv_le`, `norm1_le_half_one_add_sqrt`, `columnWitness` (a definition), `form_columnWitness`, `norm1_columnWitness`, `cond1_unitScale_le_sharp`, `form_eq_sum_eigenvalues`, `eigenvector_facts`, `lamMin_eq_iInf_eigenvalues`, `lamMax_eq_iSup_eigenvalues`. Also `cond1_unitScale_le`, which was proved again. Each prints only `propext`, `Classical.choice` and `Quot.sound`.

**The line the defects brief can take** (slice 6), replacing the earlier one's constant:
> The screen's bound is proved with a sharp column step (`cond1_unitScale_le_sharp`): κ₁(E(c)) ≤ q·((1 + √n)/2)·(M/m)·‖E(c_ref)⁻¹‖₁ in exact arithmetic, so the screen admits when q·(1 + √n)·(M/m)·‖E_ref⁻¹‖₁ ≤ 10¹⁴, with the 2 for rounding. The column constant is attained, so a smaller factor would have to use how E(c) and E(c_ref) are coupled. Searches over every family read that part near 1, which is open as a conjecture.

**Gates:**
- The checks pass 23 of 23.
- `lake build Enthymemata.PenalizedAdmission`: 2,634 jobs, with no warnings.
- `scripts/meta-ci.ps1 -NoBuild`: PASS, with `PenalizedAdmission.lean` `PROOF-CLOSED` at 89 declarations.
- No cache fetch. D: had about 33 GB free.

**Follow-ups:**
- *The conjecture:* a constant C with κ₁(E(c)) ≤ C·q·(M/m)·‖E(c′)⁻¹‖₁ for every family. It would shrink S's proved part from q·(1 + √n) to about 2C·q. Open.
- *A banded pattern's q.* A column whose entries all approach 1 needs its support to be a clique of the pattern. A band of half-width w has cliques of only w + 1, so ‖E‖₁ ≤ q may be loose on BAPS's patterns. Not posed.
- *The promotion review* of `PenalizedAdmission`.

### 2026-11-14 — targets 5 to 7 verified from a clean archive (coordinating session)

- **Commit:** `21096ad` (the Lean is `212752f`), through `scripts/verify-lean.ps1` with `-Module Enthymemata.PenalizedAdmission` and `Protolemmata\BAPS\penalized-admission-checks.py`.
- **Result:** 17 of 17 checks; 2,055 jobs with no warnings and no errors; meta-ci PASS.
- **Axioms:** every theorem of the module, 65 in all, prints `propext`, `Classical.choice` and `Quot.sound` only. The report's 23 are among them.
- **Statements, read for an instance that meets them:**
  - *Target 5* asks for positive semidefinite terms and positive coefficients, nothing else. `penalized_det_eq_zero_iff` holds for any Z and any λ > 0.
  - *Target 6* asks for the bracket m·c′ₖ ≤ cₖ ≤ M·c′ₖ and positive definiteness at c′. `comb_diagPair_posDef` shows the smallest instance meets them, and `cond2_diagPair` that it attains M / m.
  - *Target 7,* `cond1_unitScale_le`, concludes κ₁(E(c)) ≤ q · √n · (M / m) · ‖E(c′)⁻¹‖₁ for any q bounding the nonzeros in a column. That is the quantity slice 6's screen computes (`SamplerAdmission`, `5f7fac5`).
- **The two departures stand.** The composed result is stated in norms alone, so it does not rest on reading the Rayleigh extremes as eigenvalues. The bracket is how the screen takes its ratios.
- **Receipts:** `artifacts/verify-lean-21096ad/`, with the axioms file and the summary tracked in [receipts](../../../issues/regression/receipts/verify-lean-21096ad/).
- **Verified again at `94f6b37`,** after `0be098f` added the smallest instance as Lean `example`s and four helper declarations: 17 of 17 checks, 2,055 jobs, meta-ci PASS, and all 69 theorems on the three standard axioms. The seven examples discharge each main theorem's hypotheses at S₁ = diag(1, 0), S₂ = diag(0, 1), c′ = (1, 1), c = (2, 3), and BAPS's form at Z = (1 0), P = diag(0, 1). Receipts in [verify-lean-94f6b37](../../../issues/regression/receipts/verify-lean-94f6b37/).

### 2026-11-14 — targets 5 to 7 proved: a penalized system's admission

**Landed:** `212752f`.
- The module is `lean/Enthymemata/PenalizedAdmission.lean`, namespace `PenalizedAdmission`, with its line in `Enthymemata.lean`.
- The statements, checks and prior art are in [penalized-admission-lemmas.md](../BAPS/penalized-admission-lemmas.md), with [its checks](../BAPS/penalized-admission-checks.py) (17 of 17). They were posed at `a78d3d0`, and target 7 was restated at `6ec7a8c`.

**What is proved.**
- *Target 5:*
  - `form_comb_eq_zero_iff` and `comb_mulVec_eq_zero_iff`: for c > 0, ker A(c) = ⋂ₖ ker Sₖ;
  - `comb_posDef_iff` and `comb_posDef_iff_of_pos`: positive definite at one c > 0 exactly when at every one, exactly when the Sₖ share no null vector but 0;
  - `penalized_det_eq_zero_iff`: ZᵀZ + Σⱼ λⱼPⱼ is singular exactly when some x ≠ 0 has Zx = 0 and every Pⱼx = 0.
- *Target 6:*
  - `form_comb_bound`, the quadratic-form line;
  - `lamMin_comb_ge` and `lamMax_comb_le`, the eigenvalue lines;
  - `cond2_comb_le`, κ₂(A(c)) ≤ (M/m)·κ₂(A(c′));
  - `cond2_diagPair`: on S₁ = diag(1, 0), S₂ = diag(0, 1) with c′ = (1, 1), κ₂(A(c′)) = 1 and κ₂(A(c)) = max c / min c, so M/m is attained;
  - `bracket_one`, the coefficient held at 1.
- *Target 7, as restated:*
  - `diag_comb_bound`, the diagonal;
  - `lamMin_unitScale_ge`, λ_min(E(c)) ≥ (m/M)·λ_min(E(c′));
  - the norm lines `norm1_le_card` (‖E(c)‖₁ ≤ q), `norm1_le_sqrt_card` (‖A‖₁ ≤ √n·‖A‖₂) and `mulVec_sq_le_norm1` (‖A‖₂ ≤ ‖A‖₁ for a symmetric A);
  - `cond1_unitScale_le`, κ₁(E(c)) ≤ q·√n·(M/m)·‖E(c′)⁻¹‖₁.
- *Nothing remains open* among the three targets.

**What changed or proved false.** Nothing proved false.
- *Target 7's first statement held, and the same pieces gave a constant √q smaller.*
  - Posing found that ‖E(c)‖₁ ≤ q directly, and that the diagonal's ratio bound replaces the Jacobi step. It is the protolemma's finding 1.
  - The defects chip had the same chain, and the brief restated target 7 to it (`90fa6ce`).
  - The first statement's q·√(qn)·κ₁(E(c′)) follows from the restated one, through ‖E(c′)‖₁ ≥ 1 and q ≥ 1. It is not stated in Lean: the two steps do not come free.
- *The eigenvalue lines are proved through Rayleigh bounds,* not Mathlib's ordered eigenvalues.
  - `lamMin` and `lamMax` are the infimum and supremum of xᵀAx / xᵀx, and `cond2` is their ratio.
  - For a symmetric A these are the extreme eigenvalues and κ₂. The Lean does not prove that identification.
  - ‖·‖₂ enters only as a bound ‖Ax‖ ≤ β‖x‖, in squares, and no operator norm is defined.
- *The ratios enter as a bracket,* m·c′ₖ ≤ cₖ ≤ M·c′ₖ for every k. With c′ > 0 that is m ≤ ρₖ ≤ M, and min ρ and max ρ are the tightest choice.
- *q counts the nonzeros in every column of A(c).* E(c) has no more.

**Self-review.**
- *Pass 1, the draft against the protolemma.* The composed proof's first draft left a placeholder for ‖E(c′)⁻¹‖₁ > 0. `inv_norm1_le_lamMin` already proved that positivity on the way, so it now returns it, and no placeholder reached the repo.
- *Pass 2, can the hypotheses be met.* `comb_diagPair_posDef` shows the positive-definite hypothesis met by a family whose every term is singular. The bracket is met by m = min c and M = max c at c′ = 1, and the column count by q = 1 on that diagonal family. The smallest instance is kept as Lean `example`s at the module's end (`0be098f`). Each instantiates a main theorem on S₁ = diag(1, 0), S₂ = diag(0, 1), c′ = (1, 1), c = (2, 3), m = 2, M = 3 and q = 1, and BAPS's form on Z = (1 0), P = diag(0, 1) (`rowZ`, `penaltyP`).
- *Pass 3, the build's linters.* Six statements carried a `DecidableEq` or `Fintype` hypothesis their types do not use. Each now takes it inside the proof (`classical`, `Fintype.ofFinite`). Also one line over 100 characters.
- No hypothesis does the work. The two differences from the protolemma's wording, the Rayleigh extremes and the bracket, are stated above and in its Lean lines.

**Declarations to print axioms for** (23). Each prints only `propext`, `Classical.choice` and `Quot.sound`:
- `form_sq_le`, `mulVec_eq_zero_of_form`, `lamMin_pos`, `posDef_of_lamMin_pos`;
- `form_comb_eq_zero_iff`, `comb_mulVec_eq_zero_iff`, `comb_posDef_iff`, `comb_posDef_iff_of_pos`, `penalized_det_eq_zero_iff`;
- `form_comb_bound`, `lamMin_comb_ge`, `lamMax_comb_le`, `cond2_comb_le`, `bracket_one`, `cond2_diagPair`, `comb_diagPair_posDef`;
- `diag_comb_bound`, `lamMin_unitScale_ge`, `norm1_le_card`, `norm1_le_sqrt_card`, `mulVec_sq_le_norm1`, `inv_norm1_le_lamMin`, `cond1_unitScale_le`.

**The line the defects brief can take** (slice 6):
> The screen's bound is proved (`cond1_unitScale_le`, `lean/Enthymemata/PenalizedAdmission.lean`). With A(c_ref) positive definite, m·c_ref,ₖ ≤ cₖ ≤ M·c_ref,ₖ for every k, and at most q nonzeros in any column of A(c): κ₁(E(c)) ≤ q·√n·(M/m)·‖E(c_ref)⁻¹‖₁ in exact arithmetic. The screen leaves rounding alone to the sweep, the factor 2. A coefficient at 0 has m = 0, outside the hypotheses, and the screen hands that step to the estimate.

**The lines the extensions brief can take** (the admission input):
> - *Direction 1.* Identification does not read the coefficients (`comb_posDef_iff_of_pos`, `penalized_det_eq_zero_iff`). A(c) is positive definite at one c > 0 exactly when it is at every one. ZᵀZ + Σⱼ λⱼPⱼ is singular exactly when some x ≠ 0 has Zx = 0 and Pⱼx = 0 for every j. So one check of the design at construction stands on a proof, and a per-step positive-definiteness check can catch only a loss in the computed factor.
> - *Direction 2.* The bound's M/m is attained (`cond2_diagPair`). So the window of c a solve carries widens only with another form of the solve, never with a better inequality.

**Prior art,** read or from secondary sources, in the protolemma:
- Targets 5 and 6 are textbook linear algebra. Target 7's chain is the defects chip's (`76851ee`).
- *van der Sluis (1969), Numer. Math. 14:14–23: read.* The PDF is on the shelf as `PDFs/other/VanderSluis1963.pdf`. Theorems 4.1 (factor n) and 4.3 (factor q, for at most q nonzeros per row) are the Jacobi-scaling bound that is now off the path.
- *Bradley and Murray, arXiv 1110.2805: read* in its LaTeX source. It cites both theorems with those numbers.
- *Higham (2002), Corollary 7.6: read at second hand,* through NAG's F07FFF (DPOEQU) document (`PDFs/other/LapackF07.pdf`). It gives the factor n for the unit-diagonal scaling of a symmetric positive definite matrix. The book is not on the shelf.

**Gates:**
- The checks pass 17 of 17.
- `lean` on the draft: clean.
- `lake build Enthymemata.PenalizedAdmission`: 2,055 jobs, with no warnings.
- `scripts/meta-ci.ps1 -NoBuild`: PASS, with `PenalizedAdmission.lean` `PROOF-CLOSED` at 73 declarations.
- The barrel's full `lake build` did not run, since other chips' modules share the workspace. meta-CI's reachability check covers the new import.
- No cache fetch. D: had 33.3 GB free at the end.
- The Conformal Lean chip was told before the barrel edit.

**Follow-ups:**
- *A bridge from `lamMin` and `lamMax` to Mathlib's `Matrix.IsHermitian.eigenvalues`.* No decision waits on it, since the screen reads κ₁. Open.
- *Whether a factor below q·√n holds for every family.* The banded families used at most 1.007 of (M/m)·κ₁(E(c′)) (the protolemma's finding 2). A smaller proved factor would change S's proved part. The sweep measures the slack. Open.
- *The promotion review* of `PenalizedAdmission`.

**Doc changes owed:** none from this chip. The defects and extensions briefs take the lines above.

### 2026-11-14 — Payne–Weinberger verified from a clean archive (coordinating session)

- **Commit:** `e683f48`, through `scripts/verify-lean.ps1` with `-Module Enthymemata.ActiveSubspace` and this folder's checks.
- **Result:** 44 of 44 checks; 3,120 jobs with no warnings and no errors; meta-ci PASS, the seven modules proof-closed.
- **Axioms:** every theorem of the seven modules, 142 in all, prints `propext`, `Classical.choice` and `Quot.sound` only. The five new and five changed declarations are among them.
- **Statements:** the new and changed ones were read. No declaration name of `dd06a48` is gone, and five are added. The four changed statements are each the stronger or more general form.
- **Receipts:** `artifacts/verify-lean-e683f48/`, with the axioms file and the summary tracked in [receipts](../../../issues/regression/receipts/verify-lean-e683f48/).

### 2026-11-14 — Payne–Weinberger: the truncated normal's min, proved whole

**Landed:** `1ea5d6b`.
- `ContinuousPoincare`: new `wirtinger_dirichlet`, Wirtinger with zero ends, through the odd extension.
- `WeightedPoincare`:
  - new: `exists_dual`, `poincare_of_dual`, `integral_bochner` and `poincare_logConcave`;
  - `poincare_weighted` keeps its statement. Its proof now goes through the first three.
  - Stronger statements: `poincare_truncNormal`, `hasPoincareOn_truncNormal` and `matrix_poincare_box_mixed` carry min(σ², ((b − a)/π)²), E8's own D. `hasPoincareOn_of_density` takes any weighted inequality's constant as a hypothesis, in place of Bakry–Émery's.

**What changed or proved false.** Nothing proved false.
- *Payne–Weinberger needs no eigenfunction in one dimension.*
  - W = ∫ₐˣ g e^{−V} is zero at both ends. A zero-end inequality for W under the weight e^{V} gives the one with no end condition for h, with the same constant.
  - For V convex, z = e^{V/2}W has ∫ z′² ≤ ∫ g²e^{−V}. Wirtinger with zero ends then bounds ∫ W²e^{V} = ∫ z² by ((b − a)/π)² ∫ z′².
- *Bakry–Émery shares the duality.* κ ∫ W²e^{V} ≤ ∫ V″W²e^{V} ≤ ∫ g²e^{−V} replaces its pointwise step.
- *Four statements changed,* each to a stronger or more general one:
  - The earlier `hasPoincareOn_of_density` is the new one applied to `poincare_weighted`.
  - The earlier `poincare_truncNormal` and `hasPoincareOn_truncNormal` follow since min ≤ σ².
  - The earlier `matrix_poincare_box_mixed`, with Dᵢ = σᵢ², follows because E[J diag(D) Jᵀ] is monotone in D. That monotonicity is not stated in the Lean, and the earlier forms are not kept as declarations.
- The lemmas file poses this as AS-2(g) and updates AS-2(e) and (f).

**Declarations to print axioms for:**
- the new: `wirtinger_dirichlet`, `exists_dual`, `poincare_of_dual`, `integral_bochner` and `poincare_logConcave`;
- the changed: `poincare_weighted`, `poincare_truncNormal`, `hasPoincareOn_of_density`, `hasPoincareOn_truncNormal` and `matrix_poincare_box_mixed`.

Each prints only `propext`, `Classical.choice` and `Quot.sound`.

**The line E8's brief can take,** replacing the product step's line:
> The matrix check is proved on E8's laws with E8's own D (`matrix_poincare_box_mixed`, `lean/Enthymemata/ActiveSubspace/WeightedPoincare.lean`). Each input is uniform on its range or a normal N(μᵢ, σᵢ²) truncated to it. Dᵢ = ((bᵢ − aᵢ)/π)² for a uniform input and Dᵢ = min(σᵢ², ((bᵢ − aᵢ)/π)²) for a truncated normal. Every hypothesis is asked only on the box. A reading above 1 beyond rounding is a defect by theorem. The (bᵢ − aᵢ)²/π² half is Payne–Weinberger (`poincare_logConcave`).

**Prior art.** Payne and Weinberger (1960), for convex domains under the flat weight, from memory. The one-dimensional route was derived here, with no source read.

**Gates:**
- `lean` on the drafts: clean.
- `lake build Enthymemata.ActiveSubspace`: 3,120 jobs, with no warnings.
- `scripts/meta-ci.ps1 -NoBuild`: PASS, with all seven modules `PROOF-CLOSED`.
- No cache fetch. D: had 33.4 GB free at the end.

**Follow-ups:** *the promotion review* of the seven modules. Nothing else in this chip's targets stays open.

**Doc changes owed:** `docs/Maths/Regression/ActiveSubspaces.md` should say:
- the truncated normal's min(σ², (b − a)²/π²) is proved whole;
- the matrix check holds on E8's mixed laws with E8's own D.

### 2026-11-13 — the product step verified from a clean archive (coordinating session)

- **Commit:** `dd06a48`, through `scripts/verify-lean.ps1` with `-Module Enthymemata.ActiveSubspace` and this folder's checks.
- **Result:** 44 of 44 checks; 3,120 jobs with no warnings and no errors; meta-ci PASS, the seven modules proof-closed.
- **Axioms:** every theorem of the seven modules, 137 in all, prints `propext`, `Classical.choice` and `Quot.sound` only. `matrix_poincare_prod_box`, `hasPoincareOn_truncNormal` and `matrix_poincare_box_mixed` are among them.
- **Statements:** `HasPoincareOn`, the product theorems, `truncNormal` (Mathlib's `gaussianReal` conditioned on the interval) and `matrix_poincare_box_mixed` were read; the uniform theorems keep their statements.
- **Receipts:** `artifacts/verify-lean-dd06a48/`, with the axioms file and the summary tracked in [receipts](../../../issues/regression/receipts/verify-lean-dd06a48/).

### 2026-11-13 — the product step: any product of one-variable laws, and E8's mixed laws

**Landed:** `e79624a`.
- `ContinuousPoincare`:
  - new: `HasPoincareOn`, `hasPoincareOn_uniformIcc`, `ae_mem_uniformIcc`, `poincare_prod`, `matrix_poincare_prod`, `ae_mem_box_pi` and `matrix_poincare_prod_box`;
  - `poincare_pi`, `matrix_poincare_pi`, `ae_mem_box` and `matrix_poincare_box` keep their statements. Each is now a one-call proof of the uniform case.
- `WeightedPoincare` now imports `ContinuousPoincare`, plus Mathlib's Gaussian and conditional probability. It adds `hasPoincareOn_of_density`, `truncNormal`, `truncNormal_isProbabilityMeasure`, `ae_mem_truncNormal`, `integral_truncNormal`, `hasPoincareOn_truncNormal` and `matrix_poincare_box_mixed`.

**What changed or proved false.** Nothing proved false.
- *The one-variable inequality is now a predicate on a coordinate's law,* `HasPoincareOn ν a b c`. The product theorems take one per coordinate, each with its own constant. The proof is AS-2(d)(iii)'s, with the predicate where Wirtinger was.
- *The predicate asks h′ bounded on all of ℝ.* For a law on [a, b] that loses nothing. The box theorems reach it through `clampBox`.
- *The truncated normal is Mathlib's `gaussianReal μ v` conditioned on [a, b],* so its variance is Mathlib's v : ℝ≥0, which is σ².
- The lemmas file poses this as AS-2(f).

**Declarations to print axioms for:** the 16 theorems named above, twelve new and four proved again. `HasPoincareOn` and `truncNormal` are definitions. Each theorem prints only `propext`, `Classical.choice` and `Quot.sound`.

**The line E8's brief can take,** replacing what the earlier lines say of truncated normals:
> The matrix check is proved on E8's laws (`matrix_poincare_box_mixed`, `lean/Enthymemata/ActiveSubspace/WeightedPoincare.lean`). Each input is uniform on its range or a normal N(μᵢ, σᵢ²) truncated to it, with Dᵢ = ((bᵢ − aᵢ)/π)² for a uniform input and σᵢ² for a truncated normal. Every hypothesis is asked only on the box. A reading above 1 beyond rounding against that D is a defect by theorem. Against Dᵢ = min(σᵢ², (bᵢ − aᵢ)²/π²) it rests on Payne–Weinberger, which the checks support and the Lean does not prove.

**Gates:**
- `lean` on `ContinuousPoincare`: clean in 46 s.
- `lake build Enthymemata.ActiveSubspace`: 3,120 jobs, with no warnings. Mathlib's Gaussian modules were already in the cache.
- `scripts/meta-ci.ps1 -NoBuild`: PASS, with all seven modules `PROOF-CLOSED`.
- No cache fetch. D: had 33.5 GB free at the end.

**Follow-ups:**
- *Payne–Weinberger.* It decides whether E8 keeps min(σ², (b − a)²/π²) by theorem, or σ² alone. Open.
- *The promotion review* of the seven modules.

**Doc changes owed:** `docs/Maths/Regression/ActiveSubspaces.md` should say the matrix check holds on mixed uniform and truncated-normal inputs.

### 2026-11-13 — AS-2(e) verified from a clean archive (coordinating session)

- **Commit:** `6d6b446`, through `scripts/verify-lean.ps1` with `-Module Enthymemata.ActiveSubspace` and this folder's checks.
- **Result:** 44 of 44 checks; 2,968 jobs with no warnings and no errors; meta-ci PASS, the seven modules proof-closed.
- **Axioms:** every theorem of the seven modules, 125 in all, prints `propext`, `Classical.choice` and `Quot.sound` only. `poincare_weighted` and `poincare_truncNormal` are among them.
- **Statements:** both were read. Their hypotheses on h are asked on the interval, with kinks on a countable set, and the constant is κ⁻¹, σ² for the truncated normal.
- **Receipts:** `artifacts/verify-lean-6d6b446/`, with the axioms file and the summary tracked in [receipts](../../../issues/regression/receipts/verify-lean-6d6b446/).

### 2026-11-13 — AS-2(e), Bakry–Émery on an interval

**Landed:** `1aaa85a`, `ActiveSubspace/WeightedPoincare`, with `poincare_weighted` and `poincare_truncNormal`. It joins the barrel.

**The statement.** V has V′ and V″ everywhere, with V″ continuous and V″ ≥ κ > 0 on [a, b]. h is continuous on [a, b], with a measurable derivative off a countable set, bounded on [a, b]. If ∫ₐᵇ (h − m)e^{−V} = 0, then ∫ₐᵇ (h − m)²e^{−V} ≤ κ⁻¹ ∫ₐᵇ h′²e^{−V}. The truncated normal N(μ, σ²) on [a, b] gets σ², on every interval.

**What changed or proved false.** Nothing proved false.
- The proof is one-dimensional, with no ODE or semigroup: by parts with kinks, one identity for Ψ = V′W²e^{V}, then AM–GM pointwise. It passes through Brascamp–Lieb's form, Var h ≤ E[h′²/V″], before bounding 1/V″ by κ⁻¹. That form is not stated separately.
- The bound on h′ is asked on [a, b] only, as the review asked of the uniform's form.
- The hand-over's three errors were a sign in a `linear_combination`, and a parse. `c * ∫ x in a..b, f x + d` reads `+ d` into the integrand, so `hstep3` had stated something else. The first integral is now parenthesized.

**Declarations to print axioms for:** `poincare_weighted` and `poincare_truncNormal`. Each prints only `propext`, `Classical.choice` and `Quot.sound`.

**The line E8's brief can take,** adding to the earlier ones:
> For one normal component N(μ, σ²) truncated to [a, b], σ² is proved as that coordinate's one-variable constant (`poincare_truncNormal`, Bakry–Émery, `lean/Enthymemata/ActiveSubspace/WeightedPoincare.lean`). The (b − a)²/π² half of min(σ², (b − a)²/π²) still rests on the checks (Payne–Weinberger, open). The matrix check that mixes uniform and truncated-normal coordinates still rests on AS-2(b)'s finite form and the checks, until the product step with each coordinate's inequality as a hypothesis lands.

**Gates:**
- `lean` on the module.
- `lake build Enthymemata.ActiveSubspace`, 2,968 jobs, with no warnings.
- `scripts/meta-ci.ps1 -NoBuild`: PASS. All seven modules are `PROOF-CLOSED`.
- No cache fetch. D: had 33.5 GB free at the end.

**Follow-ups:**
- The product step with each coordinate's Poincaré inequality as a hypothesis, so uniform and truncated-normal coordinates mix in one matrix check. It is the next step here.
- Payne–Weinberger. It decides whether E8 may use min(σ², (b − a)²/π²) by theorem, or σ² alone. Open.
- The promotion review of the seven modules.

**Doc changes owed:** `docs/Maths/Regression/ActiveSubspaces.md`: the truncated normal's σ² is proved.

### 2026-11-13 — the extension verified from a clean archive (coordinating session)

- **Commit:** `36bd90e`, through `scripts/verify-lean.ps1` with `-Module Enthymemata.ActiveSubspace` and this folder's checks.
- **Checks:** 44 of 44.
- **Build:** 2,967 jobs, with no warnings and no errors.
- **Axioms:** every theorem of the six modules, 123 in all, prints `propext`, `Classical.choice` and `Quot.sound` only. The axioms file was generated from the modules' own declarations.
- **meta-ci:** PASS; the six modules read as proof-closed.
- **Statements:** the extension's headline statements were read at `779c1e0`, and the two points of the entry below were the outcome. `matrix_poincare_box`, the per-section kink set, `birth_death_flux` without `hmD` and `IsHingeSum.comp_monomial` were read again at `36bd90e`.
- **Earlier run:** `779c1e0`, before those points, passed the same gates with 114 theorems.
- **Not covered:** `WeightedPoincare.lean` (AS-2(e)), uncommitted and outside the barrel at this commit.
- **Receipts:** `artifacts/verify-lean-36bd90e/`, with the axioms file and the summary tracked in [receipts](../../../issues/regression/receipts/verify-lean-36bd90e/).

### 2026-11-13 — the review's two points on the matrix check, and two follow-ups closed

**The review (coordinating session, at `779c1e0`).** Both points are taken in the Lean.
1. *The bounds were asked on all of ℝⁿ,* where a hinge is unbounded. `matrix_poincare_box` now asks every hypothesis only on the box `Πᵢ[aᵢ, bᵢ]`. It reads f and J through `clampBox`, and `ae_mem_box` shows the box has full measure under the product of uniforms, so the covariance and energy do not move.
2. *The kink set did not depend on x,* so a rotated hinge was not covered. `poincare_pi`, `matrix_poincare_pi` and `matrix_poincare_box` now take `s : Fin n → (Fin n → ℝ) → Set ℝ`, countable for each i and x.

**Follow-ups closed:**
- `birth_death_flux` no longer assumes the death flux's measurability. `measurable_deathFlux` derives it from measurable pr, q, L_M and L_{M+1}.
- AS-4(iii)'s class-level converse is proved: `IsHingeFactor.monomial` and `IsHingeSum.comp_monomial`. Every permutation with signs and positive scales keeps the hinge class. With the informal invariance of BASS's prior, this is what makes the frame's sign convention (routed to E8's owner) free in law.

**Declarations to print axioms for:** `measurable_deathFlux`, `birth_death_flux`, `IsHingeFactor.monomial`, `IsHingeSum.comp_monomial`, `poincare_pi`, `matrix_poincare_pi`, `ae_mem_box`, `matrix_poincare_box`. Each prints only `propext`, `Classical.choice` and `Quot.sound`.

**Gates:** `lake build Enthymemata.ActiveSubspace` clean, and `scripts/meta-ci.ps1 -NoBuild` PASS.

**E8's line, corrected:**
> The matrix Poincaré check is proved under the product of uniform inputs (`matrix_poincare_box`, `lean/Enthymemata/ActiveSubspace/ContinuousPoincare.lean`). It covers every fit whose sections along each input are continuous on the box and differentiable off a countable set (which may move with the other inputs), with bounded measurable partials on the box. That includes hinges of any degree, axis-aligned or rotated, as a projected refit reads them in x.

### 2026-11-13 — the extension: four open items closed

On 2026-11-13 the user asked for four of the report's open items. They are posed in the lemmas' "Extension" section (`1ef4846`) and taken smallest first. The frame's sign convention went to E8's owner through the coordinating session; it is a `src/` item.

**Landed:**
- `42a449c`, AS-1(d): `IsOrthProj.eigenvalues_mem`, `IsOrthProj.rank_eq_trace`.
- `8d9837b`, AS-4(viii) in general: `leading_moves_of_not_conformal`, with `outer_mulVec`, `bessel_pair` and `eigvec_dot`.
- `12a63da`, AS-3(c): `mh_flux_symm`, `birthFlux_eq_deathFlux`, `birth_death_flux`.
- `0af5b2d`, AS-2(d)(i): `fourierCoeffOn_of_hasDerivAt_off_countable`, `wirtinger_periodic`, `integral_reflect`, `wirtinger_interval`.
- `0f5cdd5` and `6f08bfc`, AS-2(d)(ii) and (iii):
  - `ProductVariance`: `sq_integral_le`, `var_prod_le`, `varM_nonneg`, `varM_le`, `abs_varM_le`, `varM_comp_equiv`, `measurable_varM_section`, `coordVarPi`, `measurable_coordVar`, `integral_pi_succ`, `efron_stein_pi`.
  - `ContinuousPoincare`: `abs_integral_le_of_abs_le`, `uniformIcc`, `uniformIcc_isProbabilityMeasure`, `integral_uniformIcc`, `varM_uniformIcc`, `integral_integral_update`, `poincare_pi`, `varM_eq_sub`, `covM`, `energyM`, `integral_sq_lin`, `covM_form`, `energyM_form`, `matrix_poincare_pi`.
  - `ProductVariance` joins the barrel.

**Declarations to print axioms for:** `IsOrthProj.rank_eq_trace`, `leading_moves_of_not_conformal`, `mh_flux_symm`, `birth_death_flux`, `wirtinger_interval`, `efron_stein_pi`, `poincare_pi`, `matrix_poincare_pi`. Every declaration printed only `propext`, `Classical.choice` and `Quot.sound`.

**What changed or proved false.** Nothing proved false. Two things were added that the plan did not foresee:
- *A Fourier integration by parts with a derivative off a countable set.* A hinge's kinks need it, and Mathlib's version asks a one-sided derivative at every point.
- *Efron–Stein proved directly by induction over coordinates.* The plan had it through `CommutingMeans`, but its full mean, P_∅, needs a marginal for signed integrals, which Mathlib carries only for nonnegative ones (`lmarginal`).

**Hypotheses that remain, by statement:**
- `birth_death_flux` takes the death flux's measurability (`hmD`), which follows from measurable pr, q and L_M and is not derived.
- `wirtinger_interval` and `matrix_poincare_pi` take a bounded measurable derivative, which covers hinges of any degree.
- The truncated normal's constant stays with the checks.

**The line E8's brief can take,** adding to the earlier one:
> The matrix Poincaré check is proved on the continuous product of uniforms (`matrix_poincare_pi`, `lean/Enthymemata/ActiveSubspace/ContinuousPoincare.lean`). No hypothesis is left but the fit's regularity, which hinges of any degree meet. A reading above 1 beyond rounding is a defect under uniform inputs by theorem, and under one truncated normal component by the checked bound. Also proved: a projection's rank is its trace, and every non-conformal change of units moves some C's leading direction.

**The line BPPR's brief can take,** adding to the earlier one:
> The moves balance as measures (`mh_flux_symm`, `birth_death_flux`). Dimension matching is proved as measure preservation of inserting a ridge (`measurePreserving_piFinSuccAbove`). What the code must keep is that prior and proposal densities are taken against one base measure, which slice 2's Geweke gate exercises.

**Gates:**
- `lean` on each module.
- `lake build Enthymemata.ActiveSubspace`, 2,967 jobs, clean with no warnings after the linter passes.
- `scripts/meta-ci.ps1 -NoBuild`: PASS. All six modules are `PROOF-CLOSED`.
- One scoped cache fetch, `Mathlib.Analysis.Fourier.AddCircle`: 14 files, with D: at 34 GB free before and after.
- The machine ran at 100% CPU with 3 GB free throughout, so compiles took 0.5 to 22 minutes.
- My build script had shown only the last 60 lines of `lake` output, which hid early errors for two rounds. It now shows them all.

**Follow-ups,** none with a decision waiting:
- Derive `hmD` from measurable pr, q and L_M.
- Bakry–Émery for the truncated normal.
- AS-4(iii)'s class-level converse.
- The promotion review of the six modules.

**Doc changes owed,** beyond the earlier entry's:
- `docs/Maths/Regression/ActiveSubspaces.md`: the matrix check proved on the continuous uniform law.
- The page that will carry BPPR's sampler: the measure form of its moves.

### 2026-11-13 — AS-2's hypothesis corrected after review

The coordinating session reviewed the statements at `6a3268d` and returned two points, both taken.
1. **`matrix_poincare`'s `hcoord` was too strong.**
   - It bounded Varᵢ(vᵀf) at x by cᵢ(Σₐ vₐJₐᵢ(x))², pointwise in xᵢ. `condVar` does not read xᵢ, so the hypothesis demanded the bound at the least squared partial over xᵢ.
   - The one-variable inequality Varᵢ g ≤ cᵢEᵢ[(∂ᵢg)²] cannot supply that: g = xᵢ² on a symmetric box has zero partial at 0 and nonzero variance. The theorem was true, but its hypothesis could not be discharged by the inequality it names.
   - The right side is now cᵢ · `coordMean μ i (fun z => (Σₐ vₐ J z a i)²) x`. The proof passes it to `var_le_of_coord_poincare` and returns E[Eᵢ ·] to the energy form by `meanOf_coordMean`.
   - The file header, both docstrings and the lemmas' AS-2 wording now say the integrated form.
   - *Which form the checks exercised:* the integrated one throughout. AS-2c tests Var_U(h) ≤ c·E_U[h′²] over the coordinate's law, AS-2d tests the matrix form against E[J D Jᵀ] by quadrature over the whole law, and AS-2b tests Efron–Stein exactly. None tested the pointwise form, which is why the checks did not catch it.
2. **`powerSpherical_antipode` now takes 0 < κ.** At κ < 0 the zero came from rpow's convention, not from the law.

Declaration names are unchanged: only `matrix_poincare`'s hypothesis and `powerSpherical_antipode`'s hypothesis changed. Both print only `propext`, `Classical.choice` and `Quot.sound`.

Gates:
- `lean` on the module, then `lake build Enthymemata.ActiveSubspace`, clean.
- The full `lake build` (Lemmas, Theorems, Enthymemata, 3,185 jobs) passed. It ran after the Conformal chip confirmed it was idle, and was taken on the tree just before this fix. The only warnings are the existing apologies in `BARS`, `Stability` and `CumulativeField.Local`.

### 2026-11-13 — all four targets proved

**Landed.**
- Posing: `21cfcd2`, [active-subspace-lemmas.md](active-subspace-lemmas.md) and [active-subspace-checks.py](active-subspace-checks.py), 44 of 44 checks passing.
- The Lean: one commit per target, then the barrel and `Enthymemata.lean`'s line, then this report and the lemmas' status lines.
- The coordinating session read the statements at `21cfcd2` and asked for the four in order without stopping.
- AS-4 stands as posed, the frame reading first.

**Modules and the declarations to print axioms for** (namespace `ActiveSubspace`; every one prints only `propext`, `Classical.choice`, `Quot.sound`):
- `Enthymemata.ActiveSubspace.ProjectionMean` (AS-1): `projMean_posSemidef`, `one_sub_projMean_posSemidef`, `trace_projMean`, `projLoss_eq`, `trace_sub_two_trace_mul_ge`, `summaryProj_trace_sub`, `summaryLoss_le`, `half_frobSq_cross_rank`, `compress_eigenvalues_mem`.
- `Enthymemata.ActiveSubspace.MatrixPoincare` (AS-2): `loewner_iff`, `covForm_eq_var`, `energyForm_eq`, `meanOf_mul_sub_coordMean`, `var_le_of_coord_poincare`, `coordwise_le_scalar`, `matrix_poincare`, `energy_box`.
- `Enthymemata.ActiveSubspace.DirectionMoves` (AS-3): `mhRatio_of_symm`, `powerSpherical_symm`, `powerSpherical_antipode`, `change_ratio`, `change_detailed_balance`, `change_ratio_powerSpherical`, `birth_green_ratio`, `poisson_count_ratio`, `double_count_ratio`, `prior_draw_cancel`, `birth_death_balance`.
- `Enthymemata.ActiveSubspace.Refit` (AS-4):
  - the frame: `ridgeGrid_det`, `ridgeHinge_not_separable`, `hinge_isHingeSum`, `hinge_comp_not_hingeSum`, `hingeSum_not_closed`, `hingeFactor_neg`, `hingeFactor_scale`;
  - the scale: `outer_pullback`, `activeMatrix_pullback`, `pullback_eigen_of_orthogonal`, `conformal_pullback`, `pullback_ker`, `pullback_range`, `leading_C_coord`, `leading_pullback_coord`, `leading_moves`, `hat_invariant`, `rawRidge_moves`, `trainingRange_factor`, `trainingRange_hinge`, `kappa_moves`, `kappa_ratio`.
- The barrel is `Enthymemata.ActiveSubspace`, and `Enthymemata.lean` imports it.

**What remains open, as posed.**
- AS-2's one-variable Poincaré inequality is a hypothesis (`hcoord`). Mathlib carries no Poincaré or Wirtinger inequality, and the constants (b − a)²/π² and min(σ², (b − a)²/π²) rest on checks AS-2c and AS-2e.
- AS-3's measure-theoretic form on S^{a−1} × ℝ is not stated; check AS-3d reads it.
- AS-4(viii) is proved on its example; the general non-conformal L stays with the check.
- AS-4(iii)'s converse (monomial maps keep the class) is proved at the factor level only.
- That a projection's trace is its rank is used, not formalized.

**Statements that changed or proved false.** None is false. The lemmas' head lists seven findings.
- *From posing:*
  - the cross-rank distance is ½‖P_s − P‖²_F, not Ye and Lim's chordal metric;
  - Barbieri and Berger prove optimality for a different loss;
  - the change move never flips a sign at a = 1;
  - BPPR's eq. (mhb) slips at a\* = 1;
  - the scale's non-invariance is a truncation's only;
  - the hinge class is kept by monomial maps alone.
- *From the proofs:* finding 7, that a birth ratio with a further 1/(M + 1) still balances, but against the count weight λ^M/(M!)², not Poisson's (`double_count_ratio`). The protolemma's "breaks balance" wording was corrected.
- *Lean-side widenings:*
  - `half_frobSq_cross_rank` holds at every k;
  - `hat_invariant` needs no invertible BᵀB.

**The line E8's brief can take** (slice 5's items, slice 7, and the handoff's item 2):
> Lean (`lean/Enthymemata/ActiveSubspace/`, 2026-11-13):
> - **P̄ over mixed ranks:** 0 ≼ P̄ ≼ I and tr P̄ = E[k_s] are proved. The rank summary (P̄'s eigenvalues above ½, a tie at ½ either way) is least over projections of every rank (`summaryLoss_le`).
> - **The cross-rank distance** is ½‖P_s − P‖²_F = |k_s − k|/2 + Σ sin²θᵢ (`half_frobSq_cross_rank`). Ye and Lim's chordal metric weights the rank term 1. Barbieri and Berger's optimality is for predictive loss; here the ½ rule is optimal for Hamming loss outright.
> - **The matrix check:** Cov ≼ E[J D Jᵀ] follows from Efron–Stein and each input's one-variable inequality (`matrix_poincare`, the inequality a hypothesis). A reading above 1 is a defect under the uniform and one truncated normal, the laws the check covers.
> - **The frame:** the hinge class in u is kept only by signed permutations (`hinge_comp_not_hingeSum`), so the refit depends on the frame at every frame. The eigengap report warns where that choice is arbitrary; it does not say the frame suits the class.
> - **The scale:**
>   - The exact active subspace is carried by every invertible change of units (`pullback_range`). Only the leading k of a C of rank above k moves, as in `leading_moves`.
>   - Under the training-range rule the evaluated hinges do not read the declared box (`trainingRange_hinge`); κ alone does (`kappa_moves`). The g-form reads no column scale under either rule (`hat_invariant`). So the training range as default rests on a proved invariance.

**The line BPPR's brief can take** (decision 5, and slice 2's ratios):
> Lean (`lean/Enthymemata/ActiveSubspace/DirectionMoves.lean`, 2026-11-13):
> - **The change move** accepts on the collapsed likelihood ratio alone (`change_ratio_powerSpherical`), with t₀ drawn from its prior given the proposed direction. At a = 1 it keeps the direction's sign (`powerSpherical_antipode`), so the sampler treats a = 1 apart.
> - **A birth** accepts on L′/L · λ/(M + 1) · P_D(M + 1)/P_B(M) · π(a)π(𝒥|a)/(q(a)q(𝒥|a)) (`birth_green_ratio`, `poisson_count_ratio`, `prior_draw_cancel`); **a death** on its reciprocal.
> - **No further 1/(M + 1):** with it the chain targets λ^M/(M!)² (`double_count_ratio`).
> - **The paper's slip:** its eq. (mhb) divides by W_{a\*} at a\* = 1, where its Algorithm 2 draws 𝒥 uniformly.

**Prior art,** in each target's section of the lemmas:
- *Ye and Lim 2016,* read in its LaTeX source (v3; Table 2, Theorem 12, Table 3, §7).
- *Barbieri and Berger 2004,* from the abstract.
- *De Cao and Aziz 2020, and Collins, Francom and Rumsey §2.4,* read in their TeX.
- *From memory, unread:* Efron–Stein, Wirtinger, Bakry–Émery, Payne–Weinberger, Green 1995, Constantine 2015, and Draper et al. 2014. Efron–Stein is also proved here as SG-1(d).
- *The rank argument and the matrix Poincaré form* had no source found.

**Self-review, by target.**
- *AS-1:*
  - (c)'s Lean holds at every k, wider than posed; the docstring says so.
  - `compress_eigenvalues_mem` took `DecidableEq` out of its type (linter).
  - The "½ indifferent" clause is covered by `summaryProj_trace_sub`'s range of S.
- *AS-2:* this pass judged `hcoord` to be the posed hypothesis. It was not: its pointwise form was stronger than the one-variable inequality. The review entry above records the correction.
- *AS-3:*
  - The pass found finding 7 and added `poisson_count_ratio` and `double_count_ratio`.
  - `change_ratio` assumed only a symmetric kernel, so `change_ratio_powerSpherical` ties it to the power spherical law on the sphere.
  - The two-point balance device became a constant kernel.
- *AS-4:*
  - At the coordinating session's request, `leading_moves` states leading eigenvectors (`IsLeading`) and not a product. `leading_C_coord` and `leading_pullback_coord` say which.
  - `rawRidge_moves` dropped an unused `Fintype`.
  - The flexible `simp`s became `simpa`s.

**Build gates.**
- `lean` on each module.
- `lake build Enthymemata.ActiveSubspace` (the four modules and the barrel), clean with no warnings after the linter pass, at below-normal priority.
- `scripts/meta-ci.ps1 -NoBuild`: PASS. The four modules are `PROOF-CLOSED`.
- The full `lake build` was not run at first: `ConformalValidity.lean` held the Conformal chip's uncommitted work. It ran once that chip confirmed it was idle, and passed (see the entry above).
- No Mathlib fetch was needed. D: had 34 GB free at the end.

**Follow-ups,** none with a decision waiting, so recorded as open:
- A one-variable Wirtinger inequality in Lean would turn AS-2's hypothesis into a proof. Nothing in the readout changes.
- BPPR's supplement, which holds the paper's change ratio, is not on disk. Reading it would confirm the paper's own form; the derived ratio stands either way.
- The four modules are candidates for the promotion review.

**Doc changes owed.**
- `docs/Maths/Regression/ActiveSubspaces.md`, slice 5's section: the cross-rank distance as ½‖·‖²_F with Ye and Lim's weighting beside it, the median model's reason, and the matrix check's covered laws.
- The same page, slice 7's section: the frame's non-invariance (signed permutations only), the scale's (truncation only), and the training range's invariance.
