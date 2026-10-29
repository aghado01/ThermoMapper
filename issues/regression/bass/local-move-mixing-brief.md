# BASS: the local moves' mixing (F-M2 design brief)

Date: 2026-10-28. Status: a draft for the user's design pass; nothing in it is decided except who ran Phase 0. Phase 0 was measured on 2026-10-28; its readings are below, and they correct the draft in three places (marked *Corrected by Phase 0*). It opens the fitting checklist's F-M2, which the user created on 2026-10-28 (the checklist's report "S10's questions decided"). The local moves' mixing is its own design item, before the efficiency pass; the run lengths and the ladder's sizing follow it.

## What S10 found

- *Untempered, the chains do not cross between structures.* On Friedman (p = 10, five inputs inert, n = 500, σ = 1), births and deaths accept about 1.5% of attempts and changes 5%. That is about twenty changes of dimension per chain in 2,000 draws.
  - Each chain holds its own band of M (10–12, 8–10, 10–12, 9–11). R̂ of M is 2.35, and σ's interval misses the truth.
  - The mixed categorical and direct functional fits reach R̂ 1.54 and 2.55.
  - R BASS 1.3.1 behaves alike (exploration only). Every lever below goes beyond R.
- *Tempered, no round trip completes.* R̂ of M falls to 1.45 and σ is covered.
  - Every chain's ladder holds one steep band, b ≈ 0.16–0.32, rejecting 0.94–1.000 even at Δb ≈ 0.04.
  - Round trips are 0 per scan against Theorem 1's ideal of 0.00144.
  - The cold rungs accept births 1.7%, deaths 1.4% and changes 5.8%.
- *On chip E's functional fixture* no pair saturated, and round trips still ran at 0.0001–0.0008 per scan against an ideal of a few hundredths.

## The frame: two factors

Syed et al.'s ideal round-trip rate assumes efficient local exploration: between swaps, each replica's state is a fresh draw at its rung. So a tempered run's rate factors as

  observed rate = (local efficiency) × (ideal rate from the ladder's rejections).

- *The ladder sets the ideal.* A pair rejecting r adds r/(1 − r) to the ideal's inefficiency, so one saturated pair holds the ideal near 0 whatever the others do. Friedman's band does this: 0.00144.
- *The local moves set the efficiency:* how far a replica's state moves between swaps. Chip E's fixture had a fair ideal and an efficiency of a few percent at most.
- *Untempered, only the local moves act.*

Levers 1–3 below act on the efficiency, and on the untempered chain directly; lever 4 acts on the ideal. Each is read on its own factor: lever 4 by the ideal rate, levers 1–3 by observed over ideal and by the untempered chain's ESS.

## Phase 0: measure first

Each measurement is a test-side harness on the existing internals, which the tests can see. No production code changes, and the defaults stay bit-identical. The fixtures are S10's Friedman, mixed categorical and direct functional fits, at S10's seeds.

**M1. Where a sweep's cost goes.**
- *What to time.* Time each part at fixed states on the cost baseline's grid (n 200, 1,000, 5,000; M 5, 20, 50; p 2, 10):
  - a birth's term draw and column evaluation;
  - the design's proposal: the new column's cross products Zᵀz, O(N·M), and the bordered Cholesky, O(M²);
  - the target's evaluation;
  - the commit on an accept (the explicit residual);
  - the τ block (its slice's evaluations);
  - the composition per retained draw;
  - tempered, a scan's swaps and hand-overs.
- *It decides lever 1.* More structure moves per sweep pay only in the share of a sweep that is not a structure move. The baseline's fixed 50–60 µs suggests they pay at small N·M; the cross products suggest they pay little at large N·M.
- *It also says whether cost, not acceptance, comes first.* For a main-effect hinge, prefix sums of each column over the variable's sorted values give Zᵀz at any knot in O(M) rather than O(N·M). That is the efficiency pass's lever, named here so that M1 can price it.

**M2. Which births fail, and why.**
- *The harness.* From fixed states (each chain's state after burn-in, on each fixture), it proposes births and scores each with the target, committing none.
- *Each proposal's outcome:*
  - an empty draw (the intercept);
  - a refusal (support, admission, numerical, cap);
  - or the Metropolis–Hastings ratio, split into its parts: the marginal likelihood's Δℓ, the count factor p(M+1)/p(M) = (M + h₁)/((M + 1)(h₂ + 1)), the birth's within-move term −ln(M + 1) − ln q(t | S) with the term prior, and the selection ratio d_{M+1}/b_M.
- *Each proposed term is tagged* by its order, whether it uses an inert input (Friedman's truth), and its knots' quantiles and signs.
- *What it reads:*
  - the share of births lost to each cause;
  - Δℓ's distribution for terms on active inputs only, against all terms;
  - the best Δℓ among the proposals, which shows how far a random birth is from a useful one.
- *The same for deaths.* Remove each term of the state in turn and record its Δℓ. At 1.4%, every term may be load-bearing, which is a local mode.
- *It decides lever 3 and bounds its gain.* If inert inputs are the births' only loss, reweighting the main effects' variable can at most double their acceptance.

**M3. Changes by knot displacement.**
- *The harness.* From the same fixed states, it groups the change move's proposals by the knot's displacement |j′ − j|, a step in the variable's distinct-value index (a quantile step), and by whether the sign flipped.
- *What it reads:* min(1, e^{Δℓ}) by group.
- *It prices lever 2 before it is built:* a local knot step accepts at the small displacements' rate.

**M4. The chains' modes.**
- *What it records:* each chain's states over the sampling phase, with their unnormalized log posterior (log marginal + log structure prior) and their terms' variables and orders.
- *What it reads:* whether the chains hold modes of comparable posterior mass, which faster local moves would cross, or whether one chain sits many nats below, stuck until a round trip frees it.

**M5. The steep band.**
- *Per rung along the frozen ladder,* the mean and spread of M and ℓ: whether the band is where the tempered posterior's M changes regime.
- *A dense fixed ladder through b ∈ [0.1, 0.4],* with no adaptation, read at Δb ≈ 0.005–0.01. Either the band's pairs fall below 0.9, a finite barrier that the sizing under-reads, so lever 4 suffices; or they stay near 1, because the band's own rungs cannot move between regimes, so levers 1–3 must act there too.
- *Under the Gaussian reading* (lever 4a), chain 0's band, two pairs at 1.000 and 0.937, holds a barrier of at least 4.2, where Λ̂ counts 1.94. *Corrected by Phase 0:* the dense ladder read it at 2.90.

## Phase 0's readings (2026-10-28)

The chip's report, with every table, is in [its brief](local-move-measurements-brief.md#reports): commits `fc88448`–`a13eb80`, report `7f76190`, reproduced here line for line on M2–M4.
- *M1, cost.*
  - At Thin = 1 the composition takes 34–94% of a sweep, above 65% wherever n ≥ 1,000: each retained draw re-evaluates every term's column for its mean (`BassStructureMoves.ColumnMean`), O(N·M).
  - The structure move takes 6–23%, and τ's block above 5% only at n = 200.
  - So k moves per sweep cost about 1 + (k − 1)·s sweeps, with s between 0.06 and 0.23, and s grows once the column means are cached.
  - Cost does not bound the levers; acceptance does. The prefix-sum cross products are worth 0–2.7% of a sweep.
- *M2, births.*
  - Births are lost to the Occam cost, not to the inputs. An admitted birth's Δℓ has its median near −½ ln N (−3.1 at N = 500), and the prior and proposal parts ask about 2 nats more, so a random birth carries a log ratio near −5.
  - Active and inert main effects accept alike: 0.0163 against 0.0141.
  - Knot placement within an active input varies acceptance about five-fold, but its best cell accepts under 5%.
- *M2, deaths.* Every term is load-bearing: 38 of Friedman's 41 terms have death α below 10⁻⁶, at a median Δℓ of −234.
- *M3, changes.*
  - A knot step of 1–3 distinct values, sign kept, accepts 0.93–0.97 on Friedman and the mixed fit, against today's redraw at 0.04.
  - On the functional fit it is coarser: 0.62 for one step of x's 60 values, 0.37 for one of t's 20.
  - Sign flips accept at most 0.02. These are acceptances, not mixing.
- *M4, modes.* The chains hold local modes, not one diffuse mode.
  - On Friedman and the functional fit, three of four chains sit 12–42 nats below the best chain's median. They hold structures that carry the surface on different interactions, and they never reach the best chain's level.
  - Leaving such a mode means replacing a load-bearing term with its substitute at once.
- *M5, the band.*
  - It is a finite barrier that the sizing under-reads, a steep but resolvable change of regime: M falls by about two terms, and ℓ's spread peaks at b ≈ 0.24.
  - At Δb = 0.01 no pair in [0.10, 0.40] rejects above 0.46, and over that window Σr reads 4.7–5.2 against the coarse ladders' 4.0–4.3.
  - With the band opened, the ideal round-trip rate rose to 0.030–0.034 per scan, about 60 round trips in 2,000 scans, and all four chains completed none. The efficiency factor is the binding one.

What the readings change in the levers below:
- lever 2 is supported (in acceptance);
- lever 1 is cheap at Thin = 1;
- lever 3a's gain is small (×1.07 pooled);
- lever 4a's increment over-reads;
- and M4 makes a move that leaves a local mode, which this draft left out of scope, a candidate (question 7).

## The levers

- *Options first.* Each lever is an option beside the current kernel, off by default until the pass (Phase 2) decides. The current kernels stay selectable whatever it decides, as per-engine defaults keep their originals.
- *Each keeps the target exact.* The ratios below are the only changes to D01's.

### Lever 1: several structure moves per sweep

- *The kernel.* A sweep runs k structure moves on D01's schedule, then the τ block: [MoveStep × k, τ] in the block scan. k = 1 is R's, and today's.
- *The ratio is unchanged.* Each move leaves the target invariant, and so does their composition.
- *k must not depend on the state.* A fixed k, or one drawn from a fixed law, keeps the sweep invariant. A k set by the state, such as k ∝ M so that each term gets about one change per sweep, does not: a state-dependent mixture of invariant kernels need not be invariant.
- *What changes.*
  - `BurnIn` and `Samples` still count sweeps, so a sweep costs k moves and retains one draw.
  - Tempered, each replica runs k moves between swaps. That is the efficiency factor's knob: more decorrelation per swap, at k times a scan's cost.
  - The ladder's adaptation rounds are counted in scans, so they lengthen in work.
- *Its measure* is ESS per second, not per sweep; M1 predicts it.

### Lever 2: local steps inside a change

- *The kernel.* A change picks a term and one of its factors uniformly, as now.
  - For a hinge factor (v, s, j), continuous or functional, with a fixed probability ρ it steps the knot index locally: j′ = j + δ, with δ from a symmetric kernel on {−w, …, w} ∖ {0}, w fixed per variable (⌈ω·n_v⌉ for a fraction ω, at least 1).
  - Otherwise, with probability 1 − ρ, it redraws the sign and knot from the prior, as now.
  - The local step keeps the sign: a flipped hinge is a different function, not a nearby one.
- *The ratio is Δℓ alone,* as for today's change:
  - the knot prior is uniform over the n_v distinct values;
  - the step is symmetric;
  - the term's and the factor's picks are the same both ways, and ρ does not depend on the state.
- *Edges and refusals.* A step outside {0, …, n_v − 1} is a rejection (a null move), which keeps the kernel symmetric at the edges. A step onto a duplicate term, or below the support count, is refused by the design, as any change is.
- *A categorical factor's local step,* taken with the same probability ρ, toggles one level: D′ = D Δ {l}, with l uniform among the n_v levels, refused when D′ is empty or holds every level.
  - The proposal is symmetric, so the ratio is Δℓ + ln π_v(D′) − ln π_v(D), the subset measure's ratio.
  - Under R's measure (the size uniform on {1, …, n_v − 1}, then uniform among the subsets of that size) the measure's term is ln C(n_v, |D|) − ln C(n_v, |D′|); under the paper's uniform over subsets it is 0.
- *Why it should pay.* The coefficients are integrated out (D04), so a knot step's acceptance is the collapsed marginal's ratio, and a nearby knot keeps most of the fit. M3 measures it.
- *Its cost is today's change's:* one column replaced.

### Lever 3: births that favour active inputs

- *The kernel.* Births draw from a mixture, q_ε(t | S) = ε·q_D01(t | S) + (1 − ε)·q_active(t | S), with ε fixed. q_active is D01's proposal except in a main-effect birth's variable (3a) or its hinge (3b). A death's reverse probability is the same mixture for the removed term under the reduced structure, as now.
- *The ratio is D01's with q_ε in place of q_D01:*
  - for a birth, the within-move term −ln(M + 1) − ln q_ε(t | S);
  - for a death, ln q_ε(dying | S ∖ dying) + ln M.
- *Two candidates for q_active, in order of cost:*
  - *3a, use-weighted main effects.* For K = 1, the variable comes from z_u ∝ w₂ + (appearances of u), the rule D01 already applies for K ≥ 2 (Nott, Kuk and Duc), in place of uniform.
    - It is exact and costs nothing: the probability is a function of the state's counts.
    - It is positive for every variable since w₂ > 0, so it may run alone, ε = 0.
    - At D01's w₂ = 5 the tilt is mild unless the counts are large; M2 bounds its gain.
    - *Corrected by Phase 0:* the gain is −7% to +24% per chain and ×1.07 pooled. this draft's M2 bounded it by "at most double", which assumed inert births never accept; M2 measured that they accept as often as active ones.
  - *3b, residual-screened hinges.* For K = 1, (v, s, j) is drawn with probability ∝ exp(κ·g(v, s, j | r_S)), where g is the single hinge's score against the current residual: the RSS reduction a lone hinge h would give, ⟨r, h⟩²/‖h‖².
    - Prefix sums over each variable's sorted values give every knot's score in O(N), so a birth costs O(p·N) more.
    - A death's reverse needs the score under the reduced structure's residual r_{S∖dying}: O(N·M) more per death.
    - It needs ε > 0.
- *The trap to avoid* (BARS, step 3.9; DiMatteo, Genovese and Kass 2001 in the docs' lineage register).
  - A proposal concentrated on good terms gives an unhelpful term a tiny reverse-birth probability, so the term is all but never removed.
  - BARS's local births left an isolated knot a death probability near e⁻³⁷. A 10% uniform component cut the knot count's autocorrelation time 12–20×, and BARS's default became the uniform base, ε = 1.
  - Here the mixture's ε bounds the trap: q_ε ≥ ε·q_D01, so a death's ratio is at most ln(1/ε) below D01's. ε = 1 restores D01 exactly.
- *The order.* 3a goes into Phase 1. 3b is built only if M2 shows the births' loss is the knots' placement within the active inputs rather than the inputs themselves, and 3a's gain is small.
- *The seam.* This lever is F-X4's first second kernel: proposal kernels behind one contract, with exact forward and reverse probabilities. KHAOS's Bernoulli-inclusion kernel (the extensions brief) plugs into the same contract, so F-X4 is built with this lever.

### Lever 4: sizing the ladder through a steep band

- *The fault.* Λ̂ = Σᵢ rᵢ counts a pair at most 1, however steep the band between its rungs: rᵢ is a censored reading of ∫λ(β) dβ over the pair. `RecommendedLevels(Λ̂)` then asks for too few rungs, and `Place` gives the band its share of an under-read barrier. Argon's QA notes found this, and the circle-back's per-pair profile (`e69d8f0`) confirmed it.
- *Three options,* each in the shared `TemperingLadder`, generic. A ladder that never resizes stays bit-identical, BARS's included.
  - *4a, an uncensored increment.* Under the Gaussian reading of a swap (chip E's report: the log swap ratio X ~ N(−σ²/2, σ²), so r = 1 − 2Φ(−σ/2)), σ grows linearly in Δβ between close rungs, so it adds across a split pair.
    - The increment h(r) = σ(r)/√(2π) = √(2/π)·Φ⁻¹((1 + r)/2) equals r to first order and grows without bound as r → 1: h(0.9) = 1.31, h(0.99) = 2.06, h(0.999) = 2.63.
    - `Size` and `Place` use h in place of r.
    - A Rao–Blackwellized r̂ of exactly 1 needs a floor, such as 1 − r̂ ≥ 1/(2·scans).
    - *Corrected by Phase 0:* h does not predict the finer ladder's Σr better than Σr does. Over the window it over-reads by 0.6–1.8 where Σr under-reads by 0.7–0.9, and a saturated pair's h is the floor's value, not the barrier's.
  - *4b, a split rule.* Any pair above a threshold (0.9, the evidence's flag) gains a rung at its barrier midpoint each round, up to `MaxLevels`.
  - *4c, thermodynamic length.* Rungs are placed by ∫ sd_β(ℓ) dβ, from each rung's own ℓ trace (`TraceRungs` keeps the slopes under the linear form), which no pair's rejection censors. But a rung holding one regime of a two-regime band reads its sd low, so this under-reads the band the other way.
- *Its measure* is the ideal round-trip rate at the settled count, against the rungs' cost. M5 says which option can work.
- *How it relates to levers 1–3.* The ideal it raises is multiplied by the efficiency the other levers set. A raised ideal with an efficiency near 0 still completes no round trip.

## After Phase 0: proposals (2026-10-29, not decided)

Talked through with the user on 2026-10-29. Each item is a proposal for the design pass; the follow-on measurements ([their brief](mode-exit-pricing-brief.md)) price the hypotheses before Phase 1 is chosen.
- *The problem is moving between representations, not the births.*
  - At a mode the residual is mostly noise, so a random birth gains about −½ ln N and is rightly rejected, and a data-informed birth would mostly find noise.
  - M4's chains hold different representations of one surface: on Friedman, chain 1 (M = 8) carries x1·x2 with two interactions beside x1 and x2 main effects, while chain 0 (M ≈ 10) carries it with five interactions and no x2 main effect.
  - The worse modes carry co-adapted terms, each load-bearing given the others.
- *Lever 1 withdrawn as an untempered default candidate.* Composition runs only on retained sweeps, so k moves per retained draw is close to `Thin` = k, with τ updated less often. It remains a tempered knob, local moves per swap, which also spreads the hand-overs' quarter of a scan.
- *Lever 3 withdrawn for sampling,* by M2 and the reframe.
  - A burn-in-only greedy variant (D10 allows adapting during burn-in) would steer chains to good modes faster.
  - But the chains would then agree because their starts stopped being dispersed, which weakens what R̂ reports. It is kept as a thought, not a lever.
- *Lever 2 kept, with a hypothesis.* The coefficients are integrated out, so sliding knots let the fit redistribute. Two co-adapted terms could drift until one is redundant, and its death then becomes cheap. Whether such a drift path exists between the modes is untested (pricing measurement P3).
- *Lever 4: resolvable, probably not the fix* (a hypothesis).
  - Likelihood tempering softens a load-bearing term only as b·Δℓ: a death costing 234 nats at b = 1 still costs about 70 at b = 0.3.
  - So the mid-ladder rungs may be frozen too, and states cannot be adapted as they travel (P4).
- *A new family: structural moves at fixed M.* With proposals matched to the prior, the first two have ratio Δℓ alone, as today's change does.
  - *Factor prune and grow within a term.* Prune drops one of the term's K factors, uniformly. Grow adds a factor: its variable uniform among the p − K variables the term does not use, its sign and knot from the prior.
    - The term prior's ratio for the prune is [C(p, K)/C(p, K − 1)]·2·n_v, where v is the dropped factor's variable.
    - The proposals' ratio is K/((p − K + 1)·2·n_v), and C(p, K)/C(p, K − 1) = (p − K + 1)/K, so they cancel.
    - Prune at K = 1 and grow at the cap are refusals; choosing grow or prune with probability ½ each keeps the selection symmetric.
    - Under several parts, the same holds within each part, the categorical part with its subset measure in place of a sign and knot.
  - *A factor's variable change.* v → v′, uniform among the p − K unused variables, with the knot from v′'s prior: the prior's n_v/n_v′ cancels the proposals' n_v′/n_v.
  - *A pair death, and its pair birth.*
    - The death picks an unordered pair uniformly, 1/C(M, 2).
    - The reverse draws two D01 births in sequence, q(t₁ | S)·q(t₂ | S + t₁) + q(t₂ | S)·q(t₁ | S + t₂), with S the reduced structure.
    - The count prior's two-step ratio enters through the target. A jointly redundant pair gains back about two terms' Occam cost, where either alone costs hundreds of nats.
  - Each is beyond R. Their value is priced by P1 and P2.
- *The composition's column means are cached first* (question 8): `9ed600a`, bit-identical. A sweep at n = 5,000, M ≈ 50 went from 694 to 148 µs, so M1's shares, taken before it, understate the structure move's share s.

## The mode-exit pricing's readings (2026-10-29)

The chip's report, with every table, is in [its brief](mode-exit-pricing-brief.md#reports): commits `17f8b5f`–`811360e`, report `f90c3cf`. It reads one fixed state per chain at fixed τ; Friedman chain 3's burned state lay below its whole sampling phase, so it was read at its most frequent structure.
- *P1, pair deaths: none is redundant.*
  - No jointly redundant pair exists in any state. The cheapest pair costs 57–354 nats, and every pair costs more than its cheaper single.
  - The co-adaptation is substitution: each term does what the other cannot do without it.
  - The pair death is withdrawn.
- *P2, structural moves at fixed M: exits on two fixtures, none on Friedman.*
  - The prior and proposal ratios cancel exactly, as derived, so each move's ratio is Δℓ alone.
  - *Mixed chain 2* (22 nats below the others): pruning x2 from x2 · x3 · z{2} costs 3.7 nats (α 0.026), and then its twin dies with Δℓ +2.5. The two steps cost 1.2 nats in all and reach the other chains' level.
  - *Functional chain 3*: changing a factor's variable, x0 → x1, costs 3.3 nats (α up to 0.038), and then two terms die together with Δℓ +7.7. The path gains 4.3 nats.
  - *Friedman's worse chains 0 and 3* have no cheap exit: their cheapest prune, grow or change costs 9.8 and 19.2 nats (α ≤ 6·10⁻⁵).
  - Each exit rests on one state.
- *P3, drift endpoints exist; their distance differs.*
  - Against each state's own profile (a sampled state is not its knots' optimum, being 1.4–17.2 nats below it), Friedman's x1 · x2 block is re-representable in every chain: 3–6 terms die within 20 nats once the remaining knots re-optimize. The x3, x4 and x5 terms stay load-bearing.
  - The endpoints are far, 137–1,237 distinct-value steps over 500 values per input. On the functional fit's worse chains they are close, 14–53 steps.
  - This is an optimization diagnostic: it bounds where a drift could end, not that a path exists.
- *P4, the band's rungs are frozen for deaths, as b·Δℓ predicted.*
  - The median term's death α stays at or below 10⁻³ on 30–31 of the 31 rungs in [0.10, 0.40], and becomes non-negligible only at b = 0.08–0.14.
  - Births (0.025–0.028), change redraws (0.15–0.17) and prunes on surplus terms loosen there. Local knot steps accept 0.83–1.00 on every rung.

What the readings say for Phase 1, as a reading (the choice is the user's):
- *Lever 2* is cheap and accepted everywhere. Its drift endpoints are near on the functional fit and far on Friedman. At 0.6% of distinct values per step, 137–1,237 steps is a long random walk.
- *The fixed-M family* (prune and grow, a factor's variable change) opens near-free two-step exits where the worse mode differs by a factor's variable or a redundant factor.
- *Neither opens Friedman's worse modes in one or two steps.* Tempering does not either: the band's rungs are frozen for deaths, and the local steps it still accepts do not remove terms.
- *Pair deaths add nothing.*

## Phase 1: build the options

- *The order:*
  1. F-X4's seam, with D01's kernel as its first implementation, bit-identical.
  2. Levers 1, 2 and 3a, as options.
  3. Lever 4's chosen option, after M5.
- *Exactness, checked against ground truth from outside the sampler, for each new kernel:*
  - under a constant likelihood, the count and order laws against F-O3's enumeration (`bass_enumeration_oracle.R`, the `laws` mode);
  - with the likelihood, on oracle A's 2,263-state fixture at a fixed τ, the structure kernel's transition probabilities between every pair of states, enumerated: π(S | τ, y)K(S, S′) = π(S′ | τ, y)K(S′, S) to rounding, which is detailed balance exactly;
  - and the chain's visit frequencies against the exact posterior, within their Monte Carlo error.
- *Bit-identical defaults:* every option off reproduces today's draws.

## Phase 2: one pre-registered comparison

- *Arms:*
  - the baseline;
  - each lever alone (1, 2, 3a);
  - all on, and all on less each one;
  - F-M1's two schedules (equal over available, prior-cancelling), on the baseline and on all on (question 1);
  - tempered against untempered at equal work (question 4).
- *Fixtures:* S10's three, plus Friedman at n = 2,000. The Occam cost per term grows as about ½ ln N, so the moves' acceptance shifts with n.
- *Seeds:* five master seeds per fixture and arm, four chains each, from exact prior draws (S7's starts).
- *What it reads:*
  - the rank-normalized bulk ESS per second of M, σ², the log posterior and the fitted values at the monitored points, and the minimum of them;
  - R̂ of each;
  - tempered, the round trips per scan and observed over ideal.
- *The decision rule, fixed before any run.* For example, a lever becomes the default when its geometric-mean ratio of the minimum ESS per second, over the fixtures, clears a threshold and no fixture falls below a floor. Otherwise it stays an option. The thresholds are the user's (question 3).
- *After the pass:*
  - the run lengths (S10's question 3), sized from the winner's ESS per sweep to a target ESS, with the R̂ stop;
  - then the ladder's other sizing question: it settles above need on the functional fixture.

## Not in scope

- *The efficiency pass's cost work,* such as prefix-sum cross products and allocations, beyond M1's reading, unless M1 shows that cost rather than acceptance bounds the levers.
- *A variational reference for tempering* (Surjanovic et al. 2022), which shortens the barrier by changing the path rather than the rungs. It is a later option if lever 4 cannot open the band.
- *An exchange move:* one term for another of the same order, at fixed dimension, with ratio Δℓ + ln π(new) − ln π(old) + ln q(old | S ∖ old) − ln q(new | S ∖ old). It is named here, but not proposed before M2. M4's local modes make it, or another move that replaces a load-bearing term at once, a candidate (question 7).
- *SMC* (stage 8) and *F-R1.*

## Open questions for the user

1. *Fold F-M1 into this pass?* I recommend it. The move schedule changes the births' and deaths' share, the quantity levers 1 and 3 change, and one pass avoids two readings that could disagree.
2. *3b in scope?* I recommend building it only on M2's evidence, given BARS's.
3. *The pass's thresholds and fixture set.*
4. *Tempering as an arm for the default?* D05 set untempered as the default until BASS's own comparison says otherwise. Tempered Friedman passed σ's criterion where untempered failed, at the cost of about 20 rungs.
5. *Lever 4's option,* after M5.
6. *Who runs it.* I recommend Phase 0 as one chip (harness only, no production change), then Phase 1 after you read its report. Decided by the user on 2026-10-28: Phase 0 runs as one chip, [its brief](local-move-measurements-brief.md).
7. *A move that leaves a local mode* (added after Phase 0). M4's chains sit in modes 12–42 nats apart whose terms are each load-bearing, so single-term moves pay tens of nats to leave them, and lever 2 only refines within them. The exchange move (Not in scope) is the simplest compound move; whether it, or tempering's round trips once the efficiency factor rises, is the route belongs to your design pass.
8. *The composition's column means* (added after Phase 0): cache them before the pass. Decided by the user on 2026-10-29: cached, bit-identical by construction. Done in `9ed600a`: five fits hashed identical against the code before it, and the Regression (828) and BASS oracle (35) suites pass. Caching changes every arm's cost alike, but lever 1's price, the structure share s, rises with it, so the pass should run on the cached cost. It is the efficiency pass's change, brought forward; whether it can be made bit-identical depends on matching the current sum's order.

9. *Price the mode exits before choosing Phase 1* (2026-10-29). Decided by the user: a chip measures P1–P4, [its brief](mode-exit-pricing-brief.md).

10. *Factor moves under several parts* (raised by the pricing chip): how the move picks its part and its factor, the pricing harness's matched selection being one choice.
11. *Pruning a part's only factor* while the term keeps another part (K = 1 → 0 in that part), the reverse of a grow into an empty part. It was the functional fit's only prune.
12. *A variable change in a part with one variable,* which has no unused variable to move to: refused, or the move unavailable there.

## Reports

Newest first.
