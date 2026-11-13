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
