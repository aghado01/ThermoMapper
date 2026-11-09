# Proposal kernels: E4's history-based and residual-informed births

Date: 2026-11-08. Status: slice 1, the history birth, landed (`b1df14e`) at the user's word, opt-in, and its measurement makes no default candidate; waiting for the coordinating session's verification, then slice 2, the residual birth.

## Why

E4 ([the extensions brief](extensions-brief.md#additional-proposal-kernels)) lists three kernels. KHAOS's slice 1 built the first, the Bernoulli-inclusion birth with both mutations, through F-X4 ([the KHAOS brief](khaos-brief.md)). The user queued the other two on 2026-11-08, through the coordinating session:
- **A history-based variable proposal** (Nott and Kohn 2005, *Biometrika* 92). Variables are weighted by inclusion frequencies accumulated over the run, aimed at many inert predictors.
- **Residual-informed proposal mixtures,** which keep exploratory support as BASS's own births keep it. This is also the third remedy KHAOS's slice 5 recorded for SMC's trap: a birth proposed from the data's residual, against births that accept about 2% on the trap fixture.

Both are engine kernels, for hinge and Legendre terms alike. Each preserves the target when its forward and reverse probabilities are exact. No kernel becomes a default without its measurement and the user's decision.

## What is already built

- **F-X4's contract,** `IBassProposalKernel` (`LegendreTermProposal.cs`): a birth's draw and its exact ln q_b(term | counts), the birth's and death's within-move ratios, and a change. Two members: `BassTermProposal` (D01's hinge proposals, after Nott, Kuk and Duc 2005) and `LegendreTermProposal` (KHAOS's). Each reads only the structure's counts, `BassTermCounts`, and serves every chain of a fit.
- **Per-chain move state.** `BassStructureMoves` lives per chain, and per rung under a ladder. Its `BassMovesSnapshot` carries it through checkpoints. `BassChainRun` has a burn-in-end hook, where the ladder freezes (`BassLadder.EndBurnIn`). SMC gives each particle its own moves.
- **The current residual.** `BassDesign.Residual` holds r = ỹ − Z_cβ̂ at the current structure, recomputed on every acceptance.
- **The gates' harnesses.**
  - Exact detailed balance against the enumerated law (`KhaosChainTests.Kernel_IsInDetailedBalance_WithTheEnumeratedLaw`), with a sampled chain's law and Geweke beside it.
  - The tempered laws, SMC's evidence and the flat histogram against the enumeration (`KhaosReachTests`, `RandomInterceptTests`).
  - The trap fixture and its harness (`KhaosSmcTrapQualificationTests`): Ishigami at n 120, p 4, d_max 8, q_max 3, under the raw ridge.

## The design

### The birth mixture (shared)

A fit opts in by declaring a mixture on its births:

q(t | S) = (1 − ω_H − ω_R)·q_base(t | S) + ω_H·q_H(t | S) + ω_R·q_R(t | S).

- q_base is the fit's current kernel, hinge or Legendre, unchanged. Its share stays positive, so every structure the base birth reaches keeps a positive proposal probability, and the mixture is bounded below by it.
- *Exactness.* The mixture's ln q is the log-sum-exp of its components' exact ln q. A death's reverse is the same sum at the reduced structure.
- *Scope.* Deaths stay uniform, and changes and mutations stay the base kernel's.
- *Off by default.* With no mixture declared, the base kernel is called as today, so every pinned hash and baseline holds.
- *Config.* A record on `BassSampling` beside `LocalSteps`: `BirthMixture(HistoryBirth? History, ResidualBirth? Residual)`, each component with its share. The shares are validated: each in (0, 1), and together below 1.

### Kernel 1: the history-based birth (slice 1)

**The weights.**
- ĥ_u is variable u's inclusion frequency: the share of the chain's past states in which some term uses u. This is Nott and Kohn's running mean ĉ, read over BASS's terms.
- A variable's weight is max(ĥ_u, δ), with δ = 0.01. This is their truncation to [δ, 1 − δ]; the upper bound is moot for a weight.
- **Hinge births.** The weights replace z_u ∝ w₂ + (terms using u) for every order, K = 1 included. K = 1 draws its variable uniformly in the base kernel.
- **Legendre births.** The weights take the place of each variable's powered count plus ε, (aⱼ^α + ε), in KHAOS's `make_weights`, whose passes bring η to sum q₀ as before. The rest of the birth (q₀, the redraws, the degree and the composition) is the base kernel's law, and so is the forward probability's form, with the Poisson-binomial over the new η.
- The order law is the base kernel's in both cases. The history touches variables only, as E4 states it.

**D10's policy.**
- The weights are recomputed at a fixed cadence of sweeps during burn-in, from the past only, and frozen at burn-in's end. Production is then a fixed Metropolis–Hastings kernel.
- So every step leaves the posterior invariant (Nott and Kohn's condition (3)), and the base share bounds the kernel below (condition (4)). Their condition (5), diminishing change, is not needed once the kernel is frozen.
- Before the first refresh, the component is the base kernel's own weights.
- *Per chain, and per rung under a ladder:* each rung learns from the states it holds. Chains stay independent. *(Built otherwise: one history per chain, the b = 1 rung's, which every rung reads; slice 1's report says why.)*
- *Checkpoints* carry the running sums and the frozen weights, so resume replays bit for bit.
- *SMC.* At each stage, the weights come from the stage's weighted population and stay fixed through that stage's mutations. This is the population-tuned kernel of adaptive SMC; the stage's kernel leaves π_b invariant given the population.
- *The flat histogram.* The weights learn in the learning phase and freeze with the multicanonical weights.
- *A readout:* the frozen ĥ, per chain, on the fit's diagnostics. It is Nott and Kohn's estimate of the marginal inclusion probabilities.

**Not built, recorded with the decision each would change.**
- *Nott and Kohn's best linear predictor,* ĉ_u + the history covariance's regression on the current γ_{−u}. It is the step after the marginal if the measurement shows chains holding different representations that the marginal weights do not move between, as Friedman's do.
- *Diminishing adaptation past burn-in.* Its stated condition is their (5): running means move by O(k/n) over k steps. It is built only if a frozen kernel is measured to lose what continued adaptation would keep.
- *History weights in the swaps and variable changes.* They follow if the history births measure well.

### Kernel 2: the residual-informed birth (slice 2)

**The component:** a parent, then one new factor chosen by its fit to the residual.
- *The parent* is drawn uniformly over the eligible members of the structure plus the intercept, P(S). A member is eligible if it can take one more continuous factor within q_max (and, on a Legendre basis, within d_max). The intercept as parent gives a main effect.
- *The new factor* f is drawn ∝ exp(½·b·Δℓ(f)), where:
  - Δℓ(f) = −(n/2)·ln(1 − ΔQ(f)/‖r‖²) is the single column's profile-likelihood gain. Cauchy–Schwarz bounds ΔQ by ‖r‖², and the gain is capped where they meet;
  - ΔQ(f) = (cᵀr)²/‖c‖², with c = B_p ∘ φ_f, the candidate's column, centred as the design reads it;
  - r = `BassDesign.Residual`;
  - b is the rung's inverse temperature: 1 for a chain and the flat histogram, and the stage's b for SMC.
  - The exponent ½ is the locally balanced √ weighting (Zanella 2020). At b near 0 the component is near uniform over candidates, as the flat tempered law is.
- *The candidates.*
  - Hinge: every variable the parent does not use, both signs, and every knot of the prior's knot set. They are scored per variable by cumulative sums along its sort order, O(n) per variable and sign.
  - Legendre: every (variable, degree) the parent does not use, within d_max − d_p, at O(n·d_max) per variable.
  - Categorical and functional factors are left to the base component.
- *The forward probability, exactly.*
  - q_R(t | S, r) = Σ over p ∈ sub(t) ∩ P(S) of (1/|P(S)|)·w_p(t ∖ p)/W_p.
  - sub(t) are t's terms with one continuous factor removed, and W_p is parent p's normalizer.
  - At most q_max parents, each one scoring pass.
- *A death's reverse* is the same sum at the reduced structure S′, with r(S′) from the candidate's factor, O(nM).
- *What it reaches.* An interaction whose parent is in the structure: Ishigami's x₁·x₃⁴ through a term in x₁. A pure interaction with no main-effect parent (y = x₁x₂ on centred inputs) stays with the base component, whose support keeps it reachable, as the modernization review's caution asks.

**F-X4's change.**
- The birth's members take a context: the counts, the residual's provider and b, in place of the counts alone. A death's context is built at S′ from the candidate, lazily, so a kernel that reads only the counts costs nothing more.
- The two existing members read the counts as before, so their draws and hashes hold.

**The reach.** It is stated per path in the slice, as slice 3 stated E21's: gated where it runs, refused with its reason where it does not.
- *Planned to run:*
  - chains, the ladder, SMC and the flat histogram;
  - any coefficient form;
  - GBASS's noise and the random-intercept block. There the score reads the working residual against raw candidate columns: an approximate score, with exact probabilities.
- *Planned to be refused:*
  - a direct functional response: its residual is a matrix over locations, and a functional factor needs a score of its own;
  - E34's K outputs, until that chip's core lands and its residual's form is read.

**Not built, recorded.**
- Residual-informed knots in the change moves.
- A sharpness other than ½.
- Candidate pools, or multiple-try Metropolis over them.
- A reference nearer the posterior for SMC, the remedy table's third row.

## The slices

One commit per concern; tests with the code they test. Each slice ends with its gates, its measurement, a fresh-eyes review and a report, and waits for the coordinating session's verification note.

| Slice | Delivers | First needs | Its own gates |
| --- | --- | --- | --- |
| 1. The mixture and the history birth | `BirthMixture`, `HistoryBirth`; the mixture as a member of F-X4; the history's state in `BassStructureMoves` and its snapshot; the refresh and the freeze; SMC's population weights; the flat histogram's learning; the readout; the docs | This plan, verified | Common gates; the measurement below on Friedman, Ishigami and the many-inert fixture |
| 2. The residual birth | `ResidualBirth`; F-X4's context; the candidates' scores; the death's reverse at S′; the reach table; the docs | Slice 1's mixture | Common gates; the score against dense inner products; the measurement below, with the trap's reach and MCMC's climb from empty |

## Gates common to both slices

- **The kernel's own law.** On term spaces small enough to list, hinge and Legendre:
  - Σ_t q(t | S) plus the empty draw is 1, at every enumerated S;
  - the drawn terms' frequencies against q, by χ².
- **Detailed balance, exactly.** The sweep's flows against the enumerated law, as `Kernel_IsInDetailedBalance_WithTheEnumeratedLaw` reads them, with the test's own enumeration of each forward probability.
  - The history kernel with a fixed, deliberately skewed ĥ, and with one learned in a burn-in.
  - The residual kernel at each enumerated S, its r recomputed densely.
- **The chain's law** by enumeration in production, for each kernel on a hinge fit and a Legendre fit; for the history kernel, after a learned freeze.
- **Geweke,** with proper priors, on a hinge design and a Legendre design.
  - For the history kernel, a fixed skewed ĥ.
  - For the residual kernel, the redrawn y moves r every step. That makes Geweke its sharpest check of the reverse probabilities.
- **Resume bit for bit:** mid-burn-in, across the freeze, and in production.
- **Off is untouched.** The pinned tier at each commit, the eight hashes and the SMC baseline (`Smc_KeepsTheBaseline`), and the whole assembly at the slice's end ([targeted runs](../../../CONTRIBUTING.md#testing-standards)).
- **Primitives** (slice 2). Each candidate's score from the cumulative-sum pass against the direct inner product, on chosen cases and a seeded sweep: tied values at knots, far ranges, a column that is zero on the cases, and r = 0.
- **The cost row.** A birth's and a death's cost with each component against the base kernel, at the cost report's base and at p 50. All on the cycle clock, marked indicative.
- **`docs/Maths`.**
  - The kernels go on BayesianRegression's BASS page, and the adaptation as built on Samplers' D10 section, each with its alternatives and implications.
  - New symbols (ω_H, ω_R, ĥ_u, Δℓ) and their crossings go in `issues/notation-ledger.md`: ω beside the local steps' window fraction, δ beside the factor border's.

## The measurements

They are pre-registered. Each states the decision it changes and what each outcome does. The settings are fixed before any run: ω = ¼ for the component measured, δ 0.01, a refresh every 100 sweeps, and the exponent ½.

**Mixing per second, against the current kernels.**
- *Fixtures.*
  - Friedman under hinge BASS, the g-form, at S10's size (p 10, five inert).
  - Ishigami under KHAOS, the raw ridge and the modified g-prior.
  - For the history kernel, also Friedman with 45 inert predictors (p 50): its own regime, which p 10 barely holds. This fixture is added by this plan, for the coordinating session to keep or strike.
- *The design.* Four seeds of four chains at equal sweeps.
- *The readings.*
  - ESS per second for D10's three groups (the log posterior, M and σ², and the monitored fitted values);
  - the verdict's R̂ by group;
  - births' acceptance;
  - seconds per sweep.
- *The decision it changes:* whether a kernel goes to the user as a default candidate for a problem class.
  - A gain beyond the seeds' spread on a fixture is reported as a candidate for that fixture's class, with the regime named.
  - No gain leaves the kernel opt-in, and the report says so.
  - No default is proposed from these fixtures alone, since efficiency winners depend on the regime.

**The trap's reach (slice 2).** Production SMC on the trap fixture, under the raw ridge: N 500 at 5 sweeps a stage, the engine's schedule and the 1.25 cap, 16 replicates each.
- *The readings,* against slice 5's table:
  - replicates reaching σ² within twice the chains';
  - M at the end;
  - ln Ẑ's spread and the flags;
  - births' acceptance.
- *MCMC's climb from empty* on the same fixture: sweeps until M and σ² reach the long chains' band.
- *The decision it changes:* the trap's remedy, beside the flags already built.
  - If the residual birth reaches the posterior in most replicates under the engine's schedule, the user decides whether SMC on Legendre fits uses it.
  - If it does so only with the cap, the pair goes to the user.
  - If neither, the remedy table's second row reads measured and insufficient, and its third row, a reference nearer the posterior, is next.

## Not in this brief

- The recorded alternatives above, each with its decision.
- **Locally balanced proposals over a whole neighbourhood** (Zanella 2020) and a change move informed by the residual. Each is a kernel of its own.
- **The direct functional response's residual birth,** until a functional fit asks for it.

## Working rules

The KHAOS brief's:
- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)).
- Other sessions are writing in the checkout: E34 in the collapsed core, the GBASS chip on E21's remaining work, and the conformal chip on E35. Before editing a file, check `git diff -- <file>`. Build from a snapshot (`git archive HEAD` plus this chip's files), and commit only compiling states.
- Builds set `ArtifactsRoot` inline into `artifacts\chip-kernels\`, which is deleted at the end.
- Heavy runs go at below-normal priority with capped threads and MSBuild node reuse off. Tests run through `scripts/test-select.ps1`: the selected classes per edit, the pinned tier per commit, and the whole assembly once at each slice's end.
- New code goes beside old. A removal is its own question to the user, with what is lost.
- A follow-up goes in a report with the decision it would change.

## Report

Under "Reports" below, newest first, per slice:
- what landed, by commit;
- the gates' numbers and the measurements' readings;
- refusals and acceptance;
- the cost row;
- decisions taken beyond this brief;
- what was left unbuilt, and why.

## Reports

### Slice 1: the birth mixture and the history birth (chip, 2026-11-08)

At the user's word ("Start slice 1"), on the plan as written.

**Landed** on `main`, nothing pushed:
- `b1df14e`: the mixture, the history birth and their gates.
  - `BirthMixture` and `HistoryBirth` on `BassSampling.Births`, validated, left out of the serialized form when null, so every existing fingerprint stands.
  - `BassBirthMixture`, `BassInclusionHistory`, `BassVariableWeights` and `IBassWeightedBirth` (`BassBirthMixture.cs`).
  - The weighted births in `BassTermProposal` and `LegendreTermProposal`. Their own paths are unchanged, and `make_weights` takes raw weights through the same arithmetic.
  - A kernel per chain's moves (`BassStructureMoves`), the history in its snapshot, the refresh and freeze in the run's burn-in, a ladder's shared history, SMC's `Prepare`, the flat histogram's learning, and the readout (`BassChainAccount.Inclusion`).
  - The gates: `HistoryBirthTests`, and the measurement's harness, `ProposalKernelMeasurementTests`.
- This commit: `BayesianRegression.md`'s section, its API, config and source-map rows; `PolynomialChaos.md`'s row; the notation ledger's entry; the measurement's table, `proposal-kernels-history.csv`; the harness's warm-up and alternating order; the slow tier; the status line and this report.

**What is built,** as an engine kernel for hinge and Legendre terms alike ([the section](../../../docs/Maths/BayesianRegression.md#bass-the-birth-mixture-and-the-history-birth-e4)).
- *The mixture:* q = (1 − ω_H)·q_base + ω_H·q_H. One uniform picks the component, ln q is the log-sum-exp of exact components, and the base share keeps support. Deaths and changes are the base kernel's.
- *The history birth:* each variable's weight is max(ĥ_u, δ), ĥ_u the share of the chain's recorded states that use it.
  - Hinge births take the weights through Wallenius' law at every order, one included, in place of z_u ∝ w₂ + a_u.
  - Legendre births take them in place of `make_weights`' (aⱼ^α + ε).
- *D10's policy:* recorded every burn-in sweep, refreshed every `Refresh` sweeps, frozen at burn-in's end. SMC reads each stage's weighted population in `Prepare`. The flat histogram learns with its weights and freezes at production.

**The gates.**

| Gate | Reading |
| --- | --- |
| The history | Records, refreshes at the cadence, the freeze's final refresh and a mid-learning copy, against counts kept in the test, exactly. A population's weighted share, exactly. A freeze before any record leaves the base law |
| The kernel's own law | Σ_t q(t \| S) plus the empty draw is 1 within 1.3·10⁻¹⁴, for the mixture and for the history component alone, at four structures each of four fixtures (three continuous inputs at K ≤ 2, categorical, functional, Legendre), with a skewed ĥ whose main first variable sits at its floor. Each ln q against the statement written out in the test, to 8.9·10⁻¹⁶. 200,000 draws against q by χ²: z −0.03, −0.57, 0.29 and −1.39 |
| Detailed balance, exactly | The sweep's flows against the enumerated law, births from the statement and ratios from the production mixture: worst gaps 1.4·10⁻¹⁴ to 5.7·10⁻¹⁴ in logs over seven cases (hinge continuous, skewed and learned; categorical; functional; three inputs; Legendre, skewed and learned). A mixture whose ratio reads another share fails it |
| The chain's law after a learned freeze | 5,000 learning sweeps, then 400,000 against the enumeration: worst \|z\| 2.07 (hinge, 2,263 structures) and 2.21 (Legendre, 466). The frozen history records nothing more |
| Geweke | Proper priors, a skewed ĥ, 200,000 marginal draws against 300,000 successive steps: worst \|z\| 2.32 on three hinge inputs at order 2, 1.45 on Legendre |
| The ladder | Its rungs share the b = 1 history through every resize of a self-sizing ladder (Legendre: 8 rungs to 5). Each rung's count law and E[ln τ] against the enumeration's tempered law after a learned freeze, within 5 errors at every rung: hinge (b 1, 0.5, 0.2; the hinge ladder resized 8 rungs to 4) and Legendre (b 1, 0.4, 0.1 and the reference rung) |
| SMC | `Prepare`'s weights are the population's weighted share to 10⁻¹⁴. ln p(y) −72.2912 ± 0.0137 against the enumeration's −72.2910 (raw ridge), and −69.5992 ± 0.0166 against −69.5869 (modified g-prior), flags clear |
| The flat histogram | Learned ln w within 0.017 of SAMC's limit, and every inclusion and structure within its allowance of the enumeration |
| Resume | Bit for bit for a chain, a ladder, SMC across stages, and the flat histogram across learning and production; each chain's frozen ĥ equal |
| Hinge BASS untouched | The Commit scope on `b1df14e`'s tree before its commit, rebuilt after HEAD moved: 1,874 of 1,874, the pinned tier among them. The whole regression assembly once at the slice's end, on `b1df14e` with the harness's warm-up: 2,676 of 2,676, the eight pinned hashes and the SMC baseline among them |

**The measurement.** Pre-registered as the brief fixed it (the table: [`proposal-kernels-history.csv`](proposal-kernels-history.csv); the page's [section](../../../docs/Maths/BayesianRegression.md#bass-the-birth-mixture-and-the-history-birth-e4) holds the full readout).
- *The design:* four master seeds of four chains, 5,000 + 5,000 sweeps, each fit on one thread, ω_H = ¼, δ 0.01, a refresh every 100 sweeps. The reading is bulk ESS per CPU second, history over base.
- *Friedman, p 10:* M 0.73 (z −4.2) and one fitted value 0.74 (z −2.2), both lower; the rest within the seeds' spread.
- *Friedman, p 50,* the history birth's own regime: one fitted value higher, 1.50 (z 2.1); M 1.19, the log posterior 0.59 and the rest within the spread.
- *Ishigami, the raw ridge:* every ratio 1.27–1.51, none beyond the spread (|z| ≤ 0.95). *The modified g-prior:* 0.70–2.20, none beyond it.
- *Friedman's structural group fails D10's verdict under both kernels* (R̂ of M 1.3–3.0). The history birth does not move chains between its representations, as the local-move mixing brief's reframe expects.
- *τ's ESS* under the g-form falls 13× and 25× on Friedman, within a base spread wide enough that neither is beyond it (z −1.7, −1.0). Recorded, not explained.
- *Births' acceptance* moves little: 1.86% → 2.03% (p 10), 2.08% → 1.49% (p 50), 1.25% → 1.49% (raw ridge), 13.8% → 13.9% (modified g).
- *A discarded run.* The first run left tiered compilation on, and its first timed fit carried the process's compilation: 6.8 CPU s, against 1.4 for the same fit with tiering off. A short discarded warm-up fit did not absorb it. The measurement was rerun whole with `DOTNET_TieredCompilation=0`, and the harness now warms each fixture and alternates the arms' order by seed. The discarded table stays in the run's artifacts.

**The decision it changes: none to a default.** No fixture gains across D10's groups, and the one gain beyond the spread is one comparison of 27, beside a loss on p 10. The history birth stays opt-in. Nott and Kohn's best linear predictor, the step after the marginal frequencies, is not pursued on this reading: it targets correlated inclusions, and Friedman's failure is the representations, which neither birth moves.

**The cost.** On the cycle clock, indicative (`HistoryBirth_CostRow`, tiered compilation off). A birth's kernel work, its draw and its ratio, rises from 1.10 to 1.81 µs at the cost report's base (the rich surface, M 20), from 1.32 to 2.31 µs at p 50, and from 6.32 to 10.71 µs on Ishigami (M 10). A death's ratio rises 2.0–3.2×, from 0.45 to 1.11 µs, 0.52 to 1.65 µs and 4.51 to 8.82 µs. The design's work, the factor and the marginal, is the same under both. A whole fit took 4–9% more CPU in the measurement.

**Decisions beyond the brief.**
- *One history per chain under a ladder,* not one per rung. Every rung reads the b = 1 rung's, which records the b = 1 state. A self-sizing ladder rebuilds rungs as it resizes, and a shared history survives that; a frozen kernel leaves each rung's law invariant either way.
- *The marginal frequencies only,* as E4's entry states the kernel. Nott and Kohn's best linear predictor is recorded below with its decision.
- *The SMC readout.* An SMC fit reports no per-chain ĥ; its kernel's last stage's weights are on `BassSmcKernel.History`, not on the fit.
- *The flat histogram's pilot ladder* does not learn: it runs the base birth, which is valid, and its pilot only places the nodes.
- *Samplers' page* has no D10 section to extend, so the adaptation as built sits in the BASS section.
- *The slow tier* was rewritten from the slice's profile at 16 s, not 10: the profile ran under load, its four stable reference tests 1.5–1.9× slower than in the last profile (median 1.6). It holds 80 methods, four of them this slice's.

**Left unbuilt, and why.**
- *Nott and Kohn's best linear predictor:* the inclusion of u predicted from the current γ₋ᵤ by the history's covariance. It would act where inclusions correlate, and this measurement gives no reason to try it (above).
- *Adaptation past burn-in,* under Nott and Kohn's condition (5): only if a frozen kernel is measured to lose what continued adaptation keeps.
- *History weights in the swaps and variable changes:* only if the history birth measures well.
