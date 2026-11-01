# The lineage gaps' Lean leads — chip brief

Date: 2026-11-01. Status: open, for a separate thread. Written at the user's request from the Francom-lineage review.

## Context

The [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md) records where the Francom lineage and the literature around it stop short. This folder's protolemmata pose the statements behind it. They have passed numeric checks ([`lineage-gap-checks.py`](lineage-gap-checks.py)) and nothing more: no statement has been compiled, and every paper-dependent claim comes from an abstract or a search snippet.

**This chip's job:** turn the self-contained statements into compiling, proved enthymemata, then run the proof-walks now that the papers are on disk, and correct the record with what it finds. A statement that turns out false, or true only with a further hypothesis, is a finding worth more than a proof. The harness exists for those.

**Read first:**
- [lean/README.md](../../README.md): the taxonomy, meta-CI's rules, import discipline.
- This folder's four protolemmata, in full: [sensitivity-generator-lemmas.md](sensitivity-generator-lemmas.md) (SG), [structure-chain-lemmas.md](structure-chain-lemmas.md) (SC), [guarded-estimator-lemmas.md](guarded-estimator-lemmas.md) (GE), [contraction-walk.md](contraction-walk.md).
- [docs/Maths/PosteriorContraction.md](../../../docs/Maths/PosteriorContraction.md), and `lean/Enthymemata/BgpCovering.lean` for the namespace, the import style and how cited theorems enter as labelled axioms.

## Environment

- **Toolchain** Lean 4.32.2 and Mathlib v4.32.2, as pinned. Don't bump them.
- **Mathlib's build output is absent** from `lean/.lake` (`.lake/packages/mathlib/.lake/build` is empty; the compressed cache in `~/.cache/mathlib` is on C:), and **D: is 98% full, with about 6.7 GB free.**
  - Fetch by scoped module only: `lake exe cache get Mathlib.RingTheory.Polynomial.ShiftedLegendre`. A bare `lake exe cache get` restores everything; leave that out.
  - Check `df -h /d` before each fetch and after. **If a fetch would leave under 3 GB free, stop and report.** Keep to what the modules need, and delete nothing to make room: disk decisions are the user's.
  - Whether the rest of `lean/` still has its build output is unknown. Build only this chip's modules (`lake build Enthymemata.BassStructureChain`), and run `scripts/meta-ci.ps1 -NoBuild` for the taxonomy and ledger gates. Run the full `lake build` only if the disk allows, and say in the report which of these ran.
- **Iterate fast:** call `lean` directly with `LEAN_PATH` from `lake env printenv LEAN_PATH`. Linters run only under `lake build <Module>`, so finish with that.
- **Quirks:** `Π` is a reserved token; identifiers with a combining macron don't parse; `omit … in` and `open … in` go before the docstring; `lean/` is LF, so write files from Python with `newline='\n'`; meta-CI's ledger counts the word for an apology in comments too, so keep it out of docstrings. `π` as a variable name is fine unless `Real` is open.
- **Papers:** arXiv sources in `D:/aghado01/codex-scientiae/supellex-private/gauntlet/bayesian/<id>/<id>-tex/`; PDFs in `…/PDFs/Bayesian/`; Ghosal & van der Vaart 2017 by chapter in `…/PDFs/GhosalVanDerVaart2017/`. `pdftotext` works, and `pypdf` keeps math glyphs.
- **Python:** `python` is a pyenv `.bat` shim that breaks on `|` in arguments, so write scripts to files. It has no numpy; `C:/Users/azrie/PDenv/PyenvPython/pyenv-win/versions/3.12.6/python.exe` does.
- **The checkout is shared with other sessions.** Commit by pathspec (`git add -- paths`, then `git commit -- paths`), leave nothing staged, work on main, push nothing. End each commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Never touch `issues/regression/opus-bayesian-evals.md`, `perplexity-*.md` or anything outside the files named here.

## Scope

New modules, flat in `lean/Enthymemata/`, namespace `BASS`, each added to `Enthymemata.lean`: `BassStructureChain`, `GuardedEstimator`, `BassSensitivity`, `BassContractionCover`. Take the statements in the protolemmata as a first draft. Their Lean has never been compiled, so expect to adjust instances and hypotheses.

### Part 1: self-contained statements, in this order

Easiest and most load-bearing first, so an early stop still leaves value.

1. **SC-1:** `mhRatio_mul_swap`, `mh_detailed_balance`.
2. **SC-2:** `gap_le_cut`, through the two computations for the indicator of A.
3. **SC-3:** `khaos_precision_entry`, then the matrix form with `Matrix.hadamard` and positive definiteness of the sum. `Matrix.PosSemidef.hadamard` exists, though the sum form should not need it.
4. **GE-1:** `crossfit_gap`, then the conditional expectation it comes from, over simple random samples of the pool outside A, on a finite pool.
5. **GE-3:** telescoping unbiasedness, λ* as the minimizer of the stated quadratic, and Horvitz–Thompson unbiasedness under independent inclusion, all as finite statements. GE-2 follows the same pattern with an abstract expectation; take it if it comes cheaply.
6. **SG-1:**
   - `coordMean` idempotent, self-adjoint and commuting;
   - `generator_eigen`;
   - the components P_u as a complete orthogonal family;
   - Σ_v total(v) = Σ_u |u|·V_u;
   - `efron_stein`;
   - (e), the factorization of U_ρ on product-form functions.

   Prefer proving (b)–(d) over an abstract commuting family of self-adjoint idempotents, then instantiating, so the continuous case can follow.
7. **SG-2(a):** the 1-Lipschitz bound on √V_u.
8. **SG-3:** `mixture_product_integral`; then the continuous form from `MeasureTheory.integral_fintype_prod_eq_prod`, if cheap.
9. **SG-4:** `shiftedLegendre_jacobi`, from `coeff_shiftedLegendre`. Then orthogonality by integration by parts and the weighted identity, as far as they go without Legendre normalization.
10. **SG-5:** the 2 × 2 counterexample, with r = 1 by Rayleigh quotient.
11. **Contraction cover:** a union of finitely many subspaces' balls, in an M-dimensional Euclidean norm, covered at radius ξr/4 by a count free of r and of any norm bound. This is the walk's central claim about local entropy, and it is checkable with no paper. Build on `BgpCovering`'s lemmas where they fit.

### Part 2: proof-walks, with the papers

Method: step through the proof against BASS's setting. The first step that doesn't transport is the lemma to pose. Ask whether the proof goes through, not only whether the statement's hypotheses hold, and add new items to the protolemma with the walk's table.

1. **Ronen et al. (`2210.09352v1`) onto BASS's structure chain,** at fixed τ (SC-0, SC-4). What is the bottleneck set, which step uses trees, and what would the barrier's growth in N have to be? Note what Kim & Roćková (`2306.00126v1`) add on the upper-bound side.
2. **Woodard, Schmidler & Huber (`Woodard2009.pdf`).** State persistence exactly, and what measuring it on BASS would take from a flat-histogram run's ln g. Resolve the question SC-4 raises: the Potts toy fits tempering's replicas at N³ for q = 3 while the theorem, as reported in abstracts, says torpid for q ≥ 3. Read the toy's tests and the docs to tell whether the sizes are too short or the settings differ. Report; change no toy.
3. **The contraction walk against Ghosal & van der Vaart 2017, Chapter 8.** Confirm or correct:
   - that Theorem 8.19(b) as printed asks only for local entropy;
   - that the Bgp fixed-constant axiom's global cover is a strengthening, so a local restatement would stay a consequence of the printed theorem (propose it; don't edit `BgpUnknownVarianceTheorem`);
   - each row of the walk's table.

   *MARS via LASSO* is not on disk yet; if it arrives during the chip, read it for an M-term approximation rate. Otherwise leave prior mass conditional, as written.
4. **Prior-art checks for posed statements.** Read Roustant et al. (`1906.09883v1`, `1612.03689v1`) for whether they state SG-1's frame; Yoo & Ghosal (`PDFs/Bayesian/Ghosal2015-1411.6716v3.pdf`) and Odin et al. (`2403.03540v2`) for SG-2's limits; Bürkner et al. (`2204.06043v2`) for SG-4. Correct each "not found" and each abstract-level summary that reading changes, in the protolemma and the register, with what the paper actually says.

### Out of scope

- Any change to `src/` or `tests/`, including the exact-gap measurement in the register's G1.
- G6 (orthogonal discrepancy), which is a derivation for later, and the KHAOS package's source.
- Promotion to `Lemmas/`: proof-closed enthymemata wait for semantic review.
- New design. Where a walk suggests one, record it as a finding.

## Self-review

After each module, and after each walk, take a fresh-eyes pass before moving on, for correctness and for coverage gaps, and fix what it finds. For a module: does the statement say what the protolemma claims, with no hypothesis doing the work? Run `#print axioms` on the main results. For a walk: does each verdict rest on a line you read? Name the passes and their fixes in the report.

## Reporting

Append a report below, newest first, under "Reports". Update each protolemma's status line as its items land. It should say:
- what is proved, by declaration name, and what remains apologizing;
- every statement that changed, why, and any that proved false;
- each walk's table, and the lemma each walk posed;
- which build gates ran, and the disk's free space at the end;
- corrections made to the register, and any that change what a gap claims;
- what stopped you, if anything did.

## Reports
