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

### Report: the oracle modes, Tasks A–E (2026-10-28)

**Landed.** All five tasks, in R only, one commit each: `6acc2ef` (A, the posterior mode), `e76da3c` (B, R's statistics), `bf608b5` (C, the rungs), `6582edb` (D, the checks script), `4a563f7` (E, a categorical predictor), and this report. `r/scripts/bass_stage4_checks.R` prints PASS throughout (39 checks, about five minutes). `bass_stage1_checks.R` still passes (97 pass, no FAIL), and the existing `swap`, `count` and `laws` modes write byte-identical JSON to HEAD's on stage 1's fixtures. `r/README.md` gains three map rows and the script's line.

**Running them** (cwd `r/`, as the bridge does):
- `Rscript oracles/bass_enumeration_oracle.R <data.csv> <out.json> posterior [key=value ...]`. The data are y, then the continuous inputs in [0, 1] (R's scaled ones), then `categorical` columns of factor codes. Keys, with their defaults:
  - the laws mode's: `M_max` 2, `K_max` 2, `npart` 3, `h1` = `h2` = 10, `s_min` 1e-10;
  - `coef` (`sampled` | `fixed`, with `tau` = 1/N), `centre` (`mean` | `flat` | `zero`), `sigma2_prior` (`scaled` | `jeffreys`), `sigma2_shape` 1e-3, `tau_shape` 0.5;
  - `likelihood` (`collapsed` | `constant`), `exact_fit_c` 1e3, `rungs` (0,0.05,0.2,0.5,1);
  - `categorical` 0, `L_max` 3, `subsets` (`size` | `uniform`).

  An unknown key is refused. A zero response (zero centre) or a constant one is refused unless `likelihood=constant`.
- `Rscript oracles/bass_statistics_oracle.R <design.csv> <out.json> [scalar] [tau=t]`, with the CSV as y then B, the intercept's ones first; or `... <A.csv> <out.json> functional <H.csv> <Y.csv> [tau=t]`.

**The posterior mode's fields** (0-based throughout):
- *Shared:* `mode`, `n` (N), `p` (the continuous predictors), `options` (every key as resolved), `n_terms`, `n_supported`, `n_states`.
- *`terms`*: every enumerated term, supported or not, as `{ cont: [[var, sign, knot]], cat: [[var, [levels]]], func: [] }`.
  - `var` is the position among the continuous or among the categorical predictors; `sign` is ±1.
  - `knot` indexes that variable's ascending distinct values, which are `BassInputs.ContinuousKnots`. This mode sorts the knots; the older modes keep R's first-appearance order, unchanged.
  - `levels` indexes the factor's ascending distinct codes.
  - Each part is sorted by variable. The term's prior and the posterior are invariant to R's column normalization, so C# columns compare by span.
- *`states`*: sorted term indices per state, in enumeration order (by M, then lexicographic). Every per-state array below is aligned with it.
- *`posterior`*:
  - `supported` (per term: nonzero count of the whole product column ≥ `npart`);
  - `data` `{ y, x, z }` (z only with factors);
  - `response` `{ centre, level, s2, sigma2_shape, sigma2_rate, n_eff, absolute }`: `level` is ȳ (0 under zero), and s² = Σ(y − level)²/N, as `BandedGPrior.Level` and `Sigma2Prior.ScaledToResponse` read it;
  - `categorical` `{ count, levels, L_max, subsets }`, null without factors;
  - `tau_prior` `{ shape, rate = shape·N }` or `{ fixed }`;
  - `n_gate_refused`, `n_exact_fit_refused`.
- *Per state:*
  - `M`; `k` (M + 1, M under flat);
  - `rss` and `Q`: ‖r‖² and the fit energy of the convention's response on its design, the residual taken directly. Under mean, y − ȳ on [1, Z_c]; under zero, y on [1, Z]; under flat, y − ȳ on Z_c;
  - `log_prior`: the laws target, Γ(h₁ + M)/(h₂ + 1)^M Π π(t), normalized over F;
  - `log_marginal`: log p(y | S) with a, σ² and τ integrated. It is absolute under the proper priors. Under flat the intercept's density is taken as 1 and under Jeffreys p(σ²) = 1/σ², so it is relative there;
  - `posterior`;
  - `tau_mean`: E[τ | S, y], or the fixed τ.
- *Laws:*
  - `count` `{ M, prior, posterior }`;
  - `order` `{ types, levels, prior, posterior }`: `types` holds each term type's [J, L], sorted by J + L and then J. `levels` holds the distinct per-type count vectors, sorted by total and then as strings (the laws mode's rule). Without factors, `types` = [[1, 0], …, [K_max, 0]] and the levels are (n₁, …, n_Kmax);
  - `inclusion`: per term, P(t ∈ S | y).
- *τ:* `tau` `{ mean, log_mean, log_sd, levels, quantiles }`, the posterior mixed over states. The levels are 0.025, 0.05, 0.25, 0.5, 0.75, 0.95 and 0.975, and the quantiles are of τ itself, solved to 10⁻¹³ in log τ. Null under fixed g.
- *Evidence:* `log_Z`, null under flat or Jeffreys.
- *`quadrature`* `{ states, sigma2_gap, tau_gap }`: the states checked (the most probable of each size) and the largest gaps (below). Null under a constant likelihood.
- *`tempered`* `{ rungs, integrated, linear }`. Each rung carries `{ b, defined, log_Z, law, count, order, tau_mean, mean_slope, var_slope }`:
  - `law` is per state; `order` is aligned with `order.levels`;
  - the slope moments are E_b and Var_b over (S, τ) of ∂_b ℓ_b: the closed form on the integrated form, log p₁ on the linear one;
  - an undefined rung carries only `{ b, defined: false, log_Z: null }`.

**The forms, as computed.** Here n_b = Nb, or Nb − 1 under flat, and d = g₁ + n_b/2:
- log p_b(y | S, τ) = C − (n_b/2) log 2π − [flat]·½ log(bN) + lgamma(d) + (k/2) log(τ/(τ+b)) − d·log(g₂ + b·S_b/2);
- S_b = ‖r‖² + Qτ/(τ+b), and C = g₁ log g₂ − lgamma(g₁), 0 under Jeffreys;
- τ ~ Gamma(a, rate aN) by 16-node Gauss–Legendre on unit panels of log τ, from log(1/N) − 37/a to log(1/N) + log(80/a);
- each rung's law is prior(S)·∫p(τ) e^{ℓ_b} dτ, with ℓ_b = log p_b (integrated) or b·log p₁ (linear), τ's prior untempered.

**The statistics oracle's fields:**
- `mode`, `n`, `m`, `N` (= nm, R's `data$n`), `k`;
- `XtX` (k × k, symmetrized from R's upper triangle), `Xty`, `yty` (R's `data$ssy`);
- `qf`, `bhat` (in the order of the design's columns), `R_diag` (diag of R's `chol(XtX)`), `log_det` = 2 Σ log R_diag;
- `positive_definite`, `getqf_refuses`, `getqf_ratio`, `getqf_refuses_at` (the first birth, counted 1-based over the columns after the intercept, at which getQf refused; 0 if none);
- `tau`, `s2_rate` (R's `curr$s2.rate`, (yᵀy − qf/(1+τ))/2, which cancels), `s2_rate_direct` (½(‖r‖² + qf·τ/(1+τ)) from the residual);
- `check` `{ xtx_rel_error, xty_rel_error, qf_rel_error, rss_direct }`.

The functional design's rows are (case, location), case outer, and y = vec(Yᵀ) (the admission oracle's convention). Each R statement is found by its left-hand side through `source_calls` and evaluated as written: `bass`'s `data$n`, `data$ssy` and the intercept's start, `birth_des` or `birth_des_func`'s `ata`, `Xta`, `aty` and their three assignments, `addBasis*`'s bookkeeping, and `getQf`'s R, `bhat` and `qf` with its rule skipped. So the statistics exist when R refuses.

**The checks** (`bass_stage4_checks.R`, all PASS). The fixtures:
- *A:* x₁ and x₂ on {0, ½, 1}, the 9-point grid three times (N = 27), y = 1 + 2x₁ + 3(x₁ − ½)₊(½ − x₂)₊ + N(0, 0.3²) at `set.seed(11)`;
- *B:* x₁ and x₂ on {0, ⅓, ⅔, 1}, the 16-point grid twice (N = 32), y = 1 + 2(x₁ − ⅓)₊ + 3(x₁ − ⅓)₊(⅔ − x₂)₊ + N(0, 0.25²) at seed 7;
- *C:* A's grid crossed with a three-level factor (N = 27), y = 1 + 1.5x₁ + 1.5(x₂ − ½)₊·1(z = 0) + N(0, 0.4²) at seed 21;
- *E:* A's inputs with y = 1 + 2(x₁ − ½)₊ exactly.

The C# bridge writes its own fixture CSV and runs the oracle live, so these recipes need no cross-language RNG.
- *A at M_max = 3* (48 terms, 24 supported, 2,263 states, 62 refused by the gate):
  - the count law M = 0..3 is 0.0000, 0.7719, 0.1757, 0.0525 against a prior of 0.6037, 0.3049, 0.0781, 0.0132;
  - log Z = −18.81207237;
  - τ's mean is 0.009950, its median 0.007527 and its 95% interval [0.000714, 0.033204];
  - normalization and views 6·10⁻¹⁶; log Z against the direct sum 0; the prior against D01 and D02 recomputed from the emitted terms 0;
  - σ² quadrature 9·10⁻¹³, the τ rule 2·10⁻¹⁵.
- *Every option on A at M_max = 2* (mean, flat, zero, Jeffreys, flat + Jeffreys, fixed g): each normalizes (≤ 2·10⁻¹⁵), its closed form holds against σ² quadrature (≤ 6·10⁻¹³) and its τ rule (≤ 4·10⁻¹⁵), and log Z is null exactly for the improper three.
- *Against the stage-1 oracles:* marginals against `bass_marginal_oracle.R`, under mean, flat, zero and Jeffreys, at the most probable states of sizes 1 and 2, within 6·10⁻¹⁴. Bayes factors against `bass_intercept_oracle.R`'s grid, under Liang and level centring with the Jeffreys σ² on six states, within 10⁻¹³.
- *Invariances:*
  - y + 10 moves neither level centring (3·10⁻¹³, log Z included) nor Liang (3·10⁻¹⁵), while R's zero centre moves by TV 0.076;
  - y × 1000 under the scaled σ² leaves the posterior (3·10⁻¹⁵) and moves log Z by exactly −N log 1000.
- *The exact-fit exclusion (E):* the default keeps 299 states, whose 27 exact fits hold all the posterior with finite marginals. Jeffreys refuses exactly those 27.
- *A constant likelihood on B's inputs:* the count and order laws equal the laws mode's target (10⁻¹⁶, 989 states), and τ's law is its prior (mean 1/N, qgamma's quantiles, 5·10⁻¹⁵ relative).
- *A categorical predictor (C, K_max = 1, M_max = 2):* 90 terms, 62 supported, 1,949 states (5 refused by the gate). The count law M = 0..2 is 0.0000, 0.7340, 0.2660. States holding a categorical main effect carry 0.1466 and those holding a mixed term 0.0794; log Z = −29.202469. It normalizes, its closed forms hold to 5·10⁻¹⁴, and the prior equals D01 recomputed from the emitted terms (2·10⁻¹⁵). Every term in F has at least 3 nonzero values of the whole product, and the 28 refused have fewer. Complementary main effects have equal marginals (0) and never share a state. The most probable state with a categorical part, rebuilt from its emitted terms, matches `bass_marginal_oracle.R` to 4·10⁻¹⁴.
- *The subset measures* (constant likelihood, M_max = 1): both measures equal D01's recomputed from the terms, at three and four levels (9·10⁻¹⁶); they agree at three levels (0) and differ at four (largest log-prior gap 0.249).
- *The rungs (B):*
  - on the integrated form, log Z_b at b = 0, 0.1, 0.5, 1 is 0, −11.332868, −22.118758, −27.113381, and on the linear form 0, −4.532283, −16.281052, −27.113381;
  - b = 1 gives the posterior and b = 0 the prior on F on both forms (0), and b = 0's count law is the laws target (10⁻¹⁶);
  - d log Z_b/db by central differences equals E_b[slope] at b = 0.1 and 0.5, to 5·10⁻⁸ (integrated) and 2·10⁻⁸ (linear);
  - the two forms differ between the ends (TV 0.234 at b = ½).
- *R's statistics:*
  - scalar, at A's most probable state of size 2: XᵀX 0, Xᵀy 10⁻¹⁶, qf 10⁻¹⁵, and qf and yᵀy − qf equal the zero-centre posterior's Q and rss (10⁻¹⁵);
  - log|XᵀX| matches `determinant`, and R's s2.rate matches the residual form (4·10⁻¹⁵);
  - functional (n = 15, m = 11, k = 5, with design-only and functional-only terms): XᵀX 2·10⁻¹⁶, aᵀYh 2·10⁻¹⁶, qf 3·10⁻¹⁵;
  - R's getQf refuses a 20-point local term of height 10⁻³ at n = 2,000 (ratio 2,903, at birth 2), which is flagged while the statistics stand (2·10⁻¹⁶).

**For the main session:**
1. **S6.**
   - *The fixture.* A is the chain fixture: its mass is spread over M = 1–3, with 5% at the cap. C is the mixed one.
   - *Ties.* Both posteriors have exact ties: x and (1 − x) span the same space with the intercept, as do a factor main effect on D and one on its complement. A's top two states are {(x₁, −, knot 2)} and {(x₁, +, knot 0)}, at 0.3859 each. Compare a chain's visits to tied states summed, or test their equality.
   - *τ.* `tau_mean` per state is E[τ | S, y] for checking the τ slice on a fixed structure, and the mixed quantiles check the chain's τ draws.
2. **S8: the integrated form near b = 0.**
   - *Its slope is singular.* At b = 0 it carries −(k/2)·E[1/τ], which is infinite under Zellner–Siow (shape ½ ≤ 1; its variance diverges up to shape 2), so the oracle emits null there. Path sampling on that form cannot use the b = 0 endpoint. This supports D05's stepping-stone default.
   - *log Z_b falls steeply there.* Even at fixed g the slope at b = 0 is about (N/2)ψ(c) ≈ −N/(2c): −13,471 on A at c = 10⁻³, against about −48 on the linear form, and on B log Z_0.05 is −8.7 of −27.1. A ladder on the integrated form needs rungs placed geometrically down to about c/N.
   - *Where the integrated form is undefined.* Under flat it needs b > (1 − 2c)/N, and under Jeffreys b > 0. These are the `defined: false` rungs, so a ladder's b_min must respect them.
3. **S4 and S5.**
   - *Parity scope.* The statistics oracle is R's uncentred convention, with the intercept a column. Use it for parity at fixed τ under the zero centre (D06: R's uncentred Gram serves comparisons with R).
   - *R's s2.rate cancels.* `s2_rate` is R's expression as written and cancels when qf ≈ yᵀy, so compare C#'s residual-based S against `s2_rate_direct`, and against `s2_rate` only where they agree.
4. **A decision for the user before S3's categorical measure: the mixed order prior.**
   - *What D01 leaves open.* D01 gives each part its own cap, but not the prior over (J, L).
   - *R's has no single target.* It draws J from I(·) over {0..maxInt} beside a categorical part, but scores J ≥ 1 with a prior of 1/maxInt and a zero-order part with `lbmcmp = 0`. Its J = 0 prior is then the proposal's count-dependent weight, and the order masses sum past 1. This is a defect, under the parity rule's incoherence test.
   - *This oracle's reading.* J uniform on {0..K_max} and L uniform on {0..L_max}, with J = L = 0 outside F. Without factors this renormalizes to D01's 1/K_max. The header flags it, and C's checks recompute the prior that way. If the user decides otherwise, `cont_lp` and `cat_lp` in the oracle are the two lines to change.
   - *Out of the posterior.* Wallenius enters only the proposal at L ≥ 2, never the posterior, so a single factor does not exercise it.
5. **Cost.** The posterior mode costs O(states × about 1,300 nodes) per rung and form. Fixture A at M_max = 3 takes about 13 s wall with R's start-up, and B about 10 s. Keep bridge fixtures at this size.
