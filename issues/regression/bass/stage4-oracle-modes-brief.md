# BASS stage 4: the oracle modes (chip brief)

Date: 2026-10-28. A chip of the [fitting segment's stage-4 plan](port-fitting-brief.md#stage-4-plan-decided-2026-10-28), which the user decided on 2026-10-28. The plan names new R oracle modes the C# steps check against:
- an exact posterior over an enumerable structure space with data, for S6's chain;
- that space's evidence and tempered normalizers, for S8's ladder;
- R's own sufficient statistics, for S4's and S5's parity.

This chip builds them in R only. The main session runs S4 in C# meanwhile, and another chip builds the shared C# additions.

Read first: the fitting segment's D01–D06 (each decided; D03's defaults and its proper-defaults subsection matter most), [the hub's parity rule](../bass-full-port-brief.md#validation-and-performance-qualification), `r/README.md`, and the oracles you extend: `r/oracles/bass_enumeration_oracle.R` (its `laws` mode already enumerates D01's term measure under D02's count prior with D06's centred angle gate), `bass_marginal_oracle.R`, `bass_intercept_oracle.R`, `bass_tau_prior_oracle.R`, and `r/scripts/bass_stage1_checks.R`.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Switch no branch. Run oracles with `r/` as the working directory; Rscript is `C:\Users\azrie\PDenv\rlang\R-4.6.1\bin\x64\Rscript.exe`. R BASS's source is sourced at the pinned commit through `bass_source` (`r/oracles/_common.R`); evaluate R's expressions as written (`source_calls`) rather than transcribing them.
- **Touch only:**
  - `r/oracles/bass_enumeration_oracle.R`, adding modes and leaving every existing mode's output unchanged;
  - new files under `r/oracles/` and `r/scripts/`;
  - your own rows in `r/README.md`'s oracle map;
  - this brief.

  The main session edits `src/**`, `tests/**`, the fitting segment and its checklist; the other chip edits `src/maths/{distributions,samplers}/**`, `src/maths/regression/gprior/**` and its own brief. You build nothing in C#.
- **`r/scripts/bass_stage1_checks.R` must still pass** with no FAIL.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with the `Co-Authored-By` line. Do not push.
- **Report** under "Reports" at the end of this brief, newest first: what landed, the commits, every emitted field's meaning and conventions (the C# tests will read them), the checks and their numbers, and anything the main session must know.

## Conventions the C# side will match

- **Indices are 0-based in JSON.** A term is emitted as `{ cont: [[var, sign, knot]], cat: [[var, [levels]]], func: [[var, sign, knot]] }`, each part sorted by variable, `knot` an index into that variable's ascending distinct values, the scaled values C#'s knots are (`BassInputs`). That is C#'s `BassTerm`.
- **A state** is a set of terms, emitted as sorted indices into a `terms` list.
- **Inputs are R's scaled inputs** in [0, 1], as the existing modes take them. Hinges are degree 1. The marginal is invariant to scaling a column (a g-prior depends on the column space only), so R's normalizing constant does not matter here; say so where it applies.

## Task A: the posterior mode (for S6)

- **What.** `mode = posterior`: the exact posterior over F, given data y and continuous inputs X.
  - F is the `laws` mode's, as decided: the support count per term, D06's order-free centred angle gate at s_min, and the cap M_max.
  - A state's prior is D01's measure under global truncation, times D02's count prior with λ integrated out.
  - Its likelihood is D03's collapsed marginal under the decided defaults, with a and σ² integrated out and τ integrated by quadrature on log τ against its prior:
    - level centring (`GPriorCentre.Mean`): the centred basis [1, Z_c], ỹ = y − ȳ, k = M + 1;
    - σ² ~ IG(c, c·s²) with c = 10⁻³ and s² the mean square about ȳ;
    - τ ~ Gamma(½, rate N/2), Zellner–Siow.
- **Options.**
  - the fixed-g mode (τ = 1/N);
  - the Liang flat intercept (k = M, N − 1 degrees of freedom);
  - R's zero centre;
  - the Jeffreys σ², with D03's exact-fit exclusion (‖r‖ ≤ c·ε·‖y‖, the oracle's c = 10³).

  Derive each closed form from D03, and check each against quadrature over σ² (and over τ where fixed).
- **Emits.**
  - Per state: its terms, log prior (normalized over F), log marginal, and posterior probability.
  - The posterior laws of M and of term orders.
  - Each term's inclusion probability.
  - τ's posterior mean and quantiles (quadrature mixed over states).
  - log Z, the evidence, which the proper defaults make absolute; report none under the improper options.
- **Size.** Keep the fixture enumerable in seconds: p = 2, three or four distinct values per input, K_max = 2, M_max = 2 or 3, N around 20 to 40, and y a known surface plus noise at a fixed seed. Emit the state count.

## Task B: R's sufficient statistics (for S4 and S5)

- **What.** A new oracle, `bass_statistics_oracle.R`, giving what R BASS computes from a design, by R's own code (parity where R is right).
  - **Scalar:** the design B (intercept first, R's uncentred convention) and y; R's `getQf` (`mcmc_funcs.R`) gives XᵀX, Xᵀy, qf = (Xᵀy)ᵀ(XᵀX)⁻¹Xᵀy and its factor's diagonal.
  - **Direct functional:** the case-side factors A and the output-side factors H. R's own expressions in `birth_des_func` (`rjmcmc_des_func.R`) give the Hadamard XᵀX = (AᵀA)∘(HᵀH) and the cross-products aᵀYh. Evaluate them as written through `source_calls`.
  - Emit also yᵀy and log|XᵀX|.
- **Not** R's acceptance terms: R's structure moves condition on σ², and D04 integrates it, so only the statistics are parity material.
- **Checks.** Against a direct solve and an explicit nm × k design; R's `getQf` refusal (its pivot-ratio rule, a defect fixture under D06) is emitted as a flag, not used.

## Task C: the evidence and tempered laws on the posterior space (for S8)

- **What.** In the posterior mode, for a list of rungs b, the exact tempered law and log normalizer log Z_b = log Σ_S prior(S)·p_b(y | S). This covers:
  - the integrated form p_b (D05: σ² integrated against the tempered likelihood; §3's closed form with τ by quadrature);
  - the linear form, b·log p₁.

  b = 0 is the prior restricted to F, whose law must equal the `laws` mode's target. Emit log(Z_b) for each b and each form, so the main session can check a ladder's stepping-stone and path-sampling estimates and each rung's visit frequencies.

## Task D: the checks script

`r/scripts/bass_stage4_checks.R`, in the style of the stage-1 script: run each mode on fixed fixtures and print PASS or FAIL for the self-checks:
- the posterior normalizes;
- under a constant likelihood it equals the `laws` mode's target law;
- one structure's marginal equals `bass_marginal_oracle.R`'s and `bass_intercept_oracle.R`'s quadratures where they overlap;
- b = 1 gives the posterior and b = 0 the prior on F;
- log Z equals a direct sum;
- the statistics oracle's qf equals a direct solve.

## Task E, if time allows: a categorical predictor in the posterior mode

One factor of three or four levels beside the continuous inputs, with D01's categorical part:
- the size-uniform subset measure, and the paper's uniform as the option;
- variables uniform at L = 1 and Wallenius at L ≥ 2;
- the case-side support counting the continuous part times the categorical part, as R's `birth_des_cat` does.

This lets S6 check a mixed-predictor chain.

## Reports
