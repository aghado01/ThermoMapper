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

### 2026-11-13 — verified from a clean archive (coordinating session)

- **Commit:** `42c980d`, through `scripts/verify-lean.ps1` with `-Module Enthymemata.ActiveSubspace` and this folder's checks.
- **Checks:** 44 of 44.
- **Build:** 2,642 jobs, with no warnings and no errors.
- **Axioms:** every theorem of the four modules, 73 in all, prints `propext`, `Classical.choice` and `Quot.sound` only. The axioms file was generated from the modules' own declarations, not from the report's list.
- **meta-ci:** PASS; the four modules read as proof-closed.
- **Statements:** read against the lemmas at `6a3268d`. The two points of the entry below were the outcome, and the archived commit holds both corrected hypotheses.
- **Receipts:** `artifacts/verify-lean-42c980d/`, with the axioms file and the summary tracked in [receipts](../../../issues/regression/receipts/verify-lean-42c980d/). The run on `6a3268d`, before the correction, passed the same gates (`artifacts/verify-lean-6a3268d/`).

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
