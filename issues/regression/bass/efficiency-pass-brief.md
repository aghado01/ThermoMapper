# BASS's efficiency pass (chip brief)

Date: 2026-11-04. Status: commissioned by the main session; the decisions below are taken. The pass acts on the cost half's ranking ([the report](performance-cost-brief.md#the-efficiency-passs-levers-ranked), 2026-10-31) lever by lever, in the recorded order, each with a before and after reading at the cost report's own points. Its report goes under "Reports" below, newest first. `docs/Maths/BayesianRegression.md`'s efficiency-pass bullet and the cost report's reading are the places it lands; the notation ledger and the shelf stay the main session's.

## Why

The cost report found where a sweep's time goes at the cap: the change proposal is the sweep (54–61% at n ≈ 1,000, 79–96% from n = 10,000), a change being a column evaluation plus that column's cross products. Four things stood out:
- columns run at 3–7× the machine's rate (1.8–2.5 ns per element pass, 5–7.5 at the categorical points), overhead-bound per element, against 0.28–0.46 ns per multiply-add for `PairwiseSum.Dot`;
- categorical membership costs 9.4–65 ns per case against about 2 for a hinge, because it is a binary search per case;
- the commit's explicit residual runs at the scalar rate (1.0–1.4 ns per N·M) and is 56% of a sweep in a 500 × 1,000 direct fit;
- τ's block and the composition are 18% of the base sweep and 36% at n = 200.

**A vectorized dot product is worth nothing** (`PairwiseSum.Dot` is at the vector rate), and this pass leaves it alone. Allocation work comes only where it falls out of these levers; the GC pause is at most 6.5% of a sweep anywhere.

Weighted fits and moving working responses landed after the report. The design re-read for this brief (2026-11-04): weights whiten a column (`W^½(z − z̄_w)`) and its Gram is built from the whitened columns, a direct functional design under weights holds stacked columns, a prior apart from the Gram adds a second centred column per term, and `ProposeFollow` re-whitens every column. Every lever below keeps all of those forms.

## The decisions

1. **Bits first.** A lever whose arithmetic is element-wise, without FMA, must leave every bit unchanged, and a before and after dump says so. A lever that reduces in a new order, or fuses a multiply-add, changes bits, and then it carries a stated tolerance against the full computation and a stated re-baseline of every pinned value it moves. Reductions keep a fixed order, as `PairwiseSum` does, so a result is the same on every machine.
2. **Vector code in the repo's pattern.** `Vector256` through the `Vector256.Add/Multiply/LoadUnsafe` surface (`PairwiseSum`, `MatrixOps`), which gives the same bits on accelerated and unaccelerated hardware. Every new primitive is checked in both modes (`DOTNET_EnableHWIntrinsic=0` for the unaccelerated one).
3. **One lever at a time, one commit per concern.** Each lever ends with its reading, its checks and a fresh-eyes review before the next starts. Tests are committed with the code they test.
4. **The chain's law.** The waiting-time chip's long chains read their results from the artifacts they write; a change here must leave a chain's law unchanged. Levers 2–5 do so by leaving the draws bit for bit. Lever 1 cannot, so it carries the law checks below.
5. **Take what the numbers support.** A lever's design is fixed by its measured break-even, not by its pricing. Where a lever does not pay at a point, that point keeps its path and the report says why.

## The points

The cost report's lever table names ten points. Each is read by its grid key through `BASS_COST_POINTS` (a tests-only addition to the cost harness: the listed points as one block "points", in the order given):

| point | key |
| --- | --- |
| base | `rich n=1000 p=10 M_max=20 K=2` |
| n 200 | `rich n=200 p=10 M_max=20 K=2` |
| n 50,000 | `rich n=50000 p=10 M_max=20 K=2` |
| M 5, n 10,000 | `rich n=10000 p=10 M_max=5 K=2` |
| M 200, n 10,000 | `rich n=10000 p=10 M_max=200 K=2` |
| p 200 | `rich n=1000 p=200 M_max=20 K=2` |
| K 5 | `rich n=1000 p=10 M_max=20 K=5` |
| L 1,024 | `categorical n=4096 p=5+2 L=1024 M_max=20 K=2` |
| direct 100 × 10 | `functional n_x=100 m=10 p=5 M_max=20 K=2` |
| direct 500 × 1,000 | `functional n_x=500 m=1000 p=5 M_max=20 K=2` |

Each point is read as the cost half read it: chain 0 burned to the declared rule, then M1's parts at its state (the median of nine batch means, round-robin with the sweep's batches). The reading per point is µs per sweep (`sweep`), the change, death and birth proposals (`change`, `death`, `birth move`), the commits, τ's block (`tau`) and the composition (`compose`), with their IQRs; the cost report's own figures for the same points stand beside them.

## The reading protocol

- **A frozen baseline.** The binary at the brief's commit is built once and copied to `artifacts\chip-efficiency-pass\baseline\`; it is never rebuilt. A lever's "after" is its commit's build in the working root.
- **Paired, in one window.** For each point, ABBA: baseline, after, after, baseline, each a separate process (`BASS_COST_BLOCK=points`, `BASS_COST_POINTS=<key>`, rows to its own folder). The figure for a side is the mean of its two medians; the A-to-A and B-to-B gaps are the noise the difference must clear. The cost report found two passes 25 minutes apart differing by a median 1.29×, so a reading from different windows is not a before and after.
- **Quiet.** A pair counts only when the harness flags none of its four runs (no other `testhost` at either end of a run), and the other processes' busy cores are recorded beside it. Anything taken under load is labelled indicative and is not in the before and after. Waiting for the window uses the Monitor tool's until-loop on `testhost`, not a sleep.
- **Below-normal priority, one thread.** The harness reads one chain on one thread; the process runs at below-normal priority. Rows append, and a pair already read is skipped on a rerun.
- **Counts need no window.** Bytes per sweep, frequencies and realized M are the same under load; bit-identity and tolerance checks run at once.
- **If no window comes** for a lever, its code and checks land, its reading is written as a command list in the report, and its before and after is marked unread. Nothing timed under load stands in.

## The levers

### 1. Prefix sums for a main-effect hinge

**What the cost report priced.** For a one-way proposal, the cross products at 16% of the base sweep (3–40% along the axes) and, with the column the same sums let a proposal skip until acceptance, 35% (14–45%).

**The identity.** For a one-way hinge `u(x) = (s(x − t))₊^q / c` and a held centred column `z`, in the cases' order on variable `v`: `Σ u z = (1/c) Σ_{s(xᵢ − t) > 0} zᵢ (s(xᵢ − t))^q`. Over the cases sorted by `x`, a table `T(z, v, s)` of that sum at every knot is a one-pass scan with anchored updates (the sum about knot k from the sum about knot k + 1: `A₀ ← A₀ + G`, `A₁ ← A₁ + δA₀`, `A₂ ← A₂ + 2δA₁ + δ²A₀`, with `G` the new group's total), which keeps every update the size of the quantity itself and avoids the cancellation of `S_q − qtS_{q−1}…` from raw prefix moments. The data-only moments (nonzero count, mass, `Σu`, `Σu²` at every knot, one table per variable and sign, positive terms only) answer the support rule and the diagonal without a column.

**What the pricing left out.** A table costs a pass over the cases per (column, variable) to build, in sorted order with a gather, against one vector dot per proposal. Whether the table pays depends on how many proposals reach one (column, variable) in the column's life. So lever 1 is built in two steps, and the second waits on the first.

*Step 1, the break-even, measured before any production change* (tests only, in a quiet window):
- R, the table's build cost over a dot's, at n = 1,000, 10,000 and 50,000, with the columns inside the cache and out of it (the cost report's dots run at 0.27–0.46 ns inside, 0.8–0.94 outside);
- H, the proposals that reach one (column, variable) in a column's life, from counters in the harness on the ten points' chains (a one-way proposal at variable v counts for every held column).

*Pre-registered rule.* A point gets the table path iff H ≥ 2R there, by the ski-rental argument: a table is built for a (column, variable) after ⌈R⌉ dot-path hits, so the worst case costs about twice the best. The gate is a hit counter per (column, variable), with ⌈R⌉ set by whether the held columns fit the cache. Where H < 2R (the reading expects the base point and p = 200 here) the dot path stays and the report says so.

*Step 2, the build,* if step 1 supports any point:
- **Where it applies:** a scalar response (its columns are held centred, or whitened under weights); a proposal whose term is one continuous hinge, from a birth, a change, or a prune to one factor; degree 1 or 2; weights of any kind that do not count rows (precision weights), the table taken over `root·z` since `Σ rootᵢ zᵢ = 0` makes the new column's centring drop out; the ridge and the g-form on the design's own Gram (the pending carries the same border and diagonal).
- **Where it cannot:** a term of order two or more or with a categorical or functional part (no one-variable sorted structure serves a product); a direct functional design (its cross products contract through the Khatri–Rao pair); degree 3 and above (the `Math.Pow` route and an expansion that loses digits with q); a prior apart from the Gram (each term also needs its prior-side column's products); weights that count rows for the support rule; and any proposal within the near-gate band below.
- **The near-gate rule.** A proposal whose bordered minimum angle is below 10⁶ × the floor, whose support sits within 1% of its minimum, or whose residual explanation passes the explicit-residual threshold (`ExplicitBelow` × 10³) takes the dot path, so the table's error never decides a boundary admission. The pending of a table proposal is flagged; **an acceptance re-derives the column, its border, diagonal and cross products by the dots and redoes the factor's bordering from them**, so the design (columns, factor, `_b`, residual) is the full computation's after every acceptance, and a snapshot's restore agrees with the live design. What the chain's state carries from the proposal (the candidate's residual and fitted sums, angles, ridge terms) holds the table's values, within the tolerance, until the next structural move.
- **Memory.** Tables are bounded by a per-design budget; past it the oldest (column, variable) tables are dropped. A restored design, a hand-over and `ProposeFollow` under moving weights start without tables.

**Bits.** Changed, and so called out: the proposal's border, diagonal, `b` entry and `zᵀr` come from the table, within the tolerance; the committed design is unchanged. A chain that never builds a table is bit for bit today's.

**Checks.**
- *Tolerance.* For every held column `z` and candidate `u`, `|T − dot| ≤ 4·K·ε·‖z‖₂‖u‖₂` with K the variable's knot count, over: the cost grid's designs at their burned states (all ten points' columns); a seeded random sweep of sizes n ∈ [3, 50,000], distributions (uniform, 99% of cases at one value, ties throughout, two clusters at the ends, integer-valued), knots at the first, last and interior ranks, both signs, both degrees, unit and precision weights, held columns including the near-collinear; each chosen case compared in the unaccelerated mode too.
- *The data-only moments* against explicit centred columns: count exactly, mass, `Σu` and `Σu²` within the same bound, the diagonal `Σu² − (Σu)²/n` within a stated cancellation factor (reported at the worst case of the sweep).
- *Admission agrees.* Over the sweep's proposals, the table path's admit/refuse equals the dot path's wherever the near-gate rule lets the table decide; the count of near-gate fallbacks is reported.
- *Acceptance restores the full computation.* After a forced acceptance of a table proposal, the design's factor, columns, `b`, coefficients and residual equal those the dot path commits, bit for bit; a snapshot taken after it restores to the same state.
- *The chain's law.* (a) On a fixture whose posterior over structures is enumerated (the `BassKernelEnumeration` fixtures), the chain with the table path forced on (gate 0) against exact π by χ², and the chain with it off by the same χ² on the same seeds; both p > 10⁻³, and the two empirical laws within the sampling band of each other. (b) Paired seeds, table path on and off, on the cost grid's base, a weighted fit and a tempered ladder: the structure sequences compared sweep by sweep (they agree until an accept decision flips at the tolerance's scale; the first difference, if any, is reported), then the retained draws' M, τ, σ² and ln marginal compared with the sampling error across seeds. (c) A long chain per path at the base point, the run lengths fixed before the run, their M, τ and ln marginal means within 4 standard errors across 8 seeds each.
- *Pinned values.* The eight `BassBitIdentityTests` hashes run with the table path at its default gate and with it forced. A hash that moves is re-baselined in the same commit, called out, with the reason (a table proposal moved a log density by less than the tolerance and a decision followed); a hash that holds at the default gate and moves when forced says exactly that.
- *Machine invariance.* The table's scan and gather use only exact per-lane operations and a fixed order; the same sweep is compared with hardware intrinsics off.

**Done:** the break-even measured and recorded; the table path built where the rule supports it, gated; every check above passing; the ten points read before and after; the points the rule leaves on the dot path named, with H and R; the near-gate fallback rate.

### 2. The column kernels

**What the cost report priced.** The column's share of a sweep, 25% at the base, 69% at M = 5 with n = 10,000, running at 3–7× the machine's rate; at most 60–85% of it is removable. Re-read after lever 1, since they overlap on the column.

**The change.** `HingeProductKernel.Apply` makes one pass per factor with a `degree` switch and `Math.Max` per element, then one division pass. Vectorize degrees 1 and 2 element-wise (`Max`-with-zero as a compare and select that returns `+0` for `−0` and passes NaN, as `Math.Max(x, 0.0)` does; `x·x` at degree 2; the factor product and the division by `c` as vector multiplies and divides), fuse the factor passes over a block of the column so the running product stays in registers, and leave the `Math.Pow` route (degree 3 and up) and the derivative route's branches scalar. The passes `BassDesign.Evaluate` makes around the kernels are the column's remaining `+4`: the support count, the mean's centring, the whitening and the prior side's centring, each element-wise (the mean and the masses are reductions through `PairwiseSum` and stay).

**Bits.** Identical. Every operation is per element and the same operation in the same order: `((1·h₁)·h₂)/c`, no FMA, `Max` with the scalar's zero sign. A reduction is not introduced.

**Checks.**
- *Bit identity of the kernel,* the new kernel against a verbatim copy of the old (in the test): degrees 1, 2 and 3; both signs; orders 1–5; knots at the ends and equal to cases' values; lengths 0–67 (every tail) and the cost grid's n; values across far ranges (`|x|` from 10⁻³⁰⁰ to 10³⁰⁰, ±0, ±∞, NaN, subnormals), comparing bit patterns; a seeded random parameter sweep of 5,000 terms; the derivative route in every variable.
- *Bit identity of the column,* `BassBasis.CaseColumn` and `OutputColumn` through the new kernels against the copy, with a categorical part, and `BassDesign.Evaluate`'s whole column (centred, whitened, with the prior side) against a copy, under unit weights, precision weights and a frozen reference.
- *The dump.* A hash of every column a design holds at the end of a chain at each of the ten points before and after is equal; the eight `BassBitIdentityTests` hashes hold; the cost harness's `Untouched` check holds.
- *Both modes.* The tests and the hashes run again with hardware intrinsics off.

**Done:** the ten points read before and after with this lever's share re-read after lever 1's outcome; the checks passing; what stays scalar (degree 3 and up, the derivative branch, the functional hinge's second pass over locations if not fused) named.

### 3. Categorical membership

**What the cost report priced.** 9.4 ns per case at L = 2 and 65 at L = 1,024 against about 2 for a hinge: `CategoricalFactor.Contains` is an `Array.BinarySearch` per case, and `SubsetDesignFactor.Value`, `Mean` and `Pair` use it too.

**The change.** The factor holds its subset as a byte table over its level codes (size the largest level + 1), built in the constructor; `Contains(level)` is `(uint)level < (uint)table.Length && table[level] != 0`, and the kernel's loop reads the table once per factor, not per case. The level list (`Levels`) stays for equality, hashing and printing.

**Bits.** Identical: membership is integer logic.

**Checks.** `Contains` against the binary search on every level in `[−2, max + 2]` for subsets of 1 to L − 1 levels at L = 2 … 1,024 and random subsets from a seeded sweep (including a level beyond the table and −1, the out-of-training code); the column through the kernel against the old loop; `SubsetDesignFactor` `Value`, `Mean`, `Pair` unchanged; the hashes hold, the categorical points' dump is equal; the membership timer of the cost harness (a replica of the old loop) is left as the before and a new timer in the same harness reads the production loop.

**Done:** the L = 1,024 point read before and after; membership per case before and after at L = 2 … 1,024.

### 4. The commit's explicit residual

**What the cost report priced.** `Refit` in `BassDesign.cs`, 4% at the base, 10% at M = 200 and 56% at the 500 × 1,000 direct point (19.7 ms per commit at 1.0–2.0 ns per N·M).

**The change.** `Fit` forms `r = ỹ − Σ_j β_j z_j` by one scalar pass over `r` per column, and a direct design's by a double loop with `ai·hj[l] − level`. Both are element-wise: the scalar form fuses blocks of columns into one pass over `r` (each element subtracts `β_j z_j` in column order, so the order of operations per element is today's) and runs on vectors; the direct form vectorizes its inner loop over locations (`bj·(ai·hj[l] − level)` as multiply, subtract, multiply, subtract). `ResidualSquares` is a reduction through `IInnerProduct`/`PairwiseSum` and stays; `Refit`'s solve and the `FittedSquares` dot stay.

**Bits.** Identical. Column-blocked fusion keeps every element's sequence of operations; there is no FMA. The one place a fused pass would change bits is accumulating `β_j z_j` before subtracting, which is not done.

**Checks.** `Fit` against a verbatim copy over a seeded sweep of N, M and shapes (scalar; direct at several n_x × m; M = 0, 1, 4, 5, 63, 200), bit-compared; the `ProposeFollow` path's residual (it calls the same `Fit`); `ExplicitResidualSquares`; the hashes hold; both modes.

**Done:** the commit part and the sweep read at the base, M 200 n 10,000, direct 100 × 10 and direct 500 × 1,000 before and after.

### 5. τ and the composition

**What the cost report priced.** 18% of the base sweep, 36% at n = 200, 14–21% at p 200, K 5 and the small direct point. τ's slice step (`BassTauBlock.Update`, `BassSweep.cs`) makes about 6 evaluations per step, each `current.WithLogTau(v)` (a state), `target.Evaluate` (an evaluation object), `Temperature.LogDensity` and the marginal; the composition is O(M²) with a fixed 1–3 µs under M ≈ 50.

**Method: find what the share is made of, then take what the numbers support.**
- Split one τ evaluation and one composition into their parts with a test-side stopwatch harness (batches of 10⁵ calls, tier 1): the closure, the state's construction, the evaluation object, `Refusal` and `LogStructurePrior`, the marginal's arithmetic (the logs and the log-gamma), and for the composition the conditional's construction, the draws and the solves.
- Decide with a pre-registered threshold: a change is built only if the split shows it removing at least 25% of τ's block or of the composition at the base point. Candidates, each tried only if the split points at it: struct-constrained generics or a direct marginal call that skips the state and the evaluation object (the runtime's dynamic PGO already devirtualizes monomorphic sites, so the gain of a generic version is measured before it is assumed); the structure-only parts of the marginal hoisted out of the slice's loop; the composition's allocations.
- Anything built keeps τ's draws bit for bit unless it evaluates the marginal in a different arithmetic order, in which case it carries a tolerance and the law check of lever 1, (b) and (c).

**Done:** the split recorded; the decision per candidate with its measured share against the threshold; what was built read before and after at the base, n 200, p 200, K 5 and the 100 × 10 direct point; what was not built, with the reading that kept it out.

## Checks and reviews common to all levers

- The eight `BassBitIdentityTests` hashes after every lever (and with intrinsics off).
- `BassDesignTests`, `BassBasisTests`, `BassDesignKernelTests`, `BassRunResumeTests`, `BassMovingResponseTests` and the other BASS fixtures pass in Release after every lever; a failure that predates the change is named, not fixed.
- The cost harness's `Untouched` check passes at the ten points.
- A seeded random parameter sweep, far ranges included, for every new numeric primitive, compared with its scalar form; chosen cases beside it.
- **Fresh-eyes review at the end of each lever,** for correctness (a tail, a zero sign, a weight form, an unaccelerated path, a path the chosen cases miss) and for coverage gaps (a design form the tests do not drive, a point the reading skips), fixed before the lever is reported.

## Sequencing and what this pass moves

- **G-7, G-8 and the statistical half** ([the statistical half](performance-statistical-brief.md)) wait on a quiet window and read better after this pass: a rerun saved. The report names which of their numbers this pass would move (the cost per sweep, and so error per unit time, at every point they read).
- **The waiting-time chip** ([its brief](waiting-time-brief.md)) runs long chains from its own artifacts folder. Levers 2–5 leave its draws bit for bit; lever 1 is a numerical change at the tolerance, which the chip is told of through the shared notes folder before lever 1 lands.
- Other sessions edit under `src/maths/regression`: `git status` and `git log` on a lever's files precede each lever, and commits stay small.

## Working rules

- **`main`,** no worktree, no subagent, nothing pushed. Commit by explicit path (`git commit -- paths`), leave nothing staged (the index is shared), end each message with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>` (the attribution the session's harness sets; the commissioning note named Opus 5.5).
- **Build and test** with `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-efficiency-pass\` on each `dotnet` command; delete that folder at the end.
- **Machine:** shared. Heavy runs go at below-normal priority with capped threads; timing blocks wait for a quiet window; anything timed under load is labelled indicative.
- **Other sessions' uncommitted work stays untouched.** At this brief's start the tree holds the user's untracked `issues/` notes and a deleted `opus-bass-completion-remainder.md`; none is mine.
- **Line endings:** briefs and C# are LF.
- **Files this pass owns:** `src/maths/regression/spline/bass/` (`DesignKernels.cs`, `DesignKernel.cs`, `BassTerm.cs`, `BassDesign.cs`, `BassSweep.cs`, `BassTarget.cs` and what lever 1 and 5 add beside them), the tests it adds under `tests/regression/`, the cost harness's point list, this brief, the efficiency-pass bullet in `docs/Maths/BayesianRegression.md`, and an appended reading in the cost report.

## Report

Under "Reports" below, newest first:
- per lever: what changed, the bits class and the tolerance or the equality it was checked to, the checks with their numbers, the before and after at the points (µs per sweep and the parts, the paired runs and their spread), what was left unbuilt and why;
- lever 1's break-even (R, H) and where the table pays;
- the pinned values re-baselined, if any, each with its reason;
- which of G-7, G-8 and the statistical half's numbers this pass would move;
- the machine, the windows, the labelled-indicative readings and any deviation from this brief.

## Reports

Newest first.

### The efficiency pass built (chip, 2026-11-04)

**What landed** (on `main`, nothing pushed):
- `55e5ac0`: this brief, and the cost harness's `BASS_COST_POINTS` (the lever table's ten points by their grid keys).
- `a9f3211`: levers 2, 3 and 4, bit for bit: `ColumnOps.cs`, `Hinge.Scale`, `HingeProductKernel.Multiply`, `CategoricalFactor`'s bitset, and `BassDesign`'s column and `Fit` passes.
- `4e62cf6`: lever 1, prefix tables (`HingeTables.cs` and `BassDesign`), with its tests and the break-even harness. `693689f` corrects it (below): tables only on a plain chain's design.
- `05ea2b5`: lever 5, `Sigma2PriorTerms`, the count prior's collection weight kept by M, and the composition's sort.
- `ee36f95`: the cost harness and this pass's timers read the thread's own cycle counter when `BASS_COST_CLOCK=cycles`; `b637c17`: the `Untouched` check's restored design asks for tables as the chain's does; `980f258`: the new numeric kernels compile fully optimized from their first call (a vector loop in tier-0 code is slower than the scalar loop it replaces; no bit moves); `29000f9`: the in-process old-against-new reading of lever 5's memos; `651cb90`: a chain whose model moves the working response between sweeps (a latent block) keeps the dots.

**The bits, by lever.**

| lever | bits | what was checked, and how far |
| --- | --- | --- |
| 2. column kernels | identical | `HingeKernel_IsTheScalarForm_BitForBit`: 12,800 random terms (degrees 1–4, orders 1–4, lengths 0–1,001, ties, far ranges `|x|` 10⁻³⁰⁰–10³⁰⁰, ±0, ±∞, NaN, subnormals), and `ColumnOps` against its scalar loops at 17 lengths; both with `DOTNET_EnableHWIntrinsic=0`; the eight pinned hashes with tables off |
| 3. categorical membership | identical | `Contains` against the binary search at every level in `[−3, L + 126]` and at `int.MinValue`, `int.MaxValue` for 660 subsets at L = 2…5,000; `ZeroOutside` against the per-case loop |
| 4. commit residual | identical | `SubtractScaledColumns` against column-by-column at 0–63 columns and 0–1,000 rows; `SubtractScaledProduct`; the direct design's `Fit` runs the same loop |
| 5. τ and the composition | identical | 59,797 marginal evaluations (50,787 proper, 9,010 improper, 16,954 under a moving scale) and the collection weight at M = 0…1,099 and past the table, bit for bit against the direct forms; the eight hashes with tables off |
| 1. prefix tables | tolerance | below |

**Lever 1.**
- *The break-even, measured first* (`BassPrefixBreakEvenTests`, tests only). One-way proposals' reach per (column, variable) in 20,000-sweep replays of burned chains, a table built at a pair's first read:

  | point | M | one-way proposals per sweep | pair reads | tables built | reads per table | tables live at once |
  | --- | --- | --- | --- | --- | --- | --- |
  | base | 20 | 0.466 | 177,023 | 3,767 | 47.0 | 159 |
  | n 200 | 13 | 0.429 | 84,648 | 9,892 | 8.6 | 152 |
  | n 50,000 (4,000 sweeps) | 20 | 0.343 | 26,046 | 573 | 45.5 | 155 |
  | M 5, n 10,000 | 5 | 0.199 | 15,896 | 329 | 48.3 | 8 |
  | M 200, n 10,000 (6,000 sweeps) | 200 | 0.312 | 372,296 | 3,841 | 96.9 | 2,000 |
  | p 200 | 9 | 0.445 | 73,800 | 32,163 | 2.3 | 843 |
  | K 5 | 20 | 0.398 | 151,088 | 3,563 | 42.4 | 156 |
  | L 1,024 | 20 | 0.318 | 60,411 | 1,128 | 53.6 | 98 |
  | S10, p 5 | 20 | 0.299 | 113,651 | 372 | 305.5 | 97 |

  R, a table's build (both signs) over one dot, read 11–17 at n = 10,000 and 50,000 (indicative, under load; the n = 1,000 reading was in tier-0 code and is discarded). The rule H ≥ 2R puts the break-even at about 22–34 reads per table: the base, K 5, L 1,024, S10, n 50,000 and M 5 and 200 pass it; n 200 (8.6) and p 200 (2.3) fail. The scan's error against the explicit dot read 1.6·10⁻¹⁷–1.1·10⁻¹⁶ of ‖z‖‖u‖ in the break-even harness and at most 93 ε in the sweep.
- *What was built, against the brief.* A pure cache. The brief's gate (build at the (r + 1)th read) makes a design's numbers depend on its history: a restored or handed-over design has no tables and would take the dots where the live one took tables, so a resumed chain would differ from the uninterrupted one in the low bits. So the choice of path is fixed by the design's size and nothing the cache holds: a plain `BassChain`'s design reads tables when the fit is scalar with at most 16 continuous predictors and at least 500 cases (the measured break-even), for the leading variables whose tables fit 64 MiB at the cap; every other design keeps the dots. The first paired readings then showed the second condition: a design pays about R dots a table before it reads one back, so a restored 100-sweep batch, a particle's mutation step, a ladder's handed-over replica and a draw's design lose by it; `693689f` turns tables on at a `BassChain`'s construction and restore only.
- *Where it applies:* a scalar response, degree 1 or 2, a proposal of one continuous hinge, unit or precision weights with the prior in the Gram, the g-form and the normalized ridge, categorical terms among the held columns. *Where it cannot:* a term of order two or more, or with a categorical or functional part; a direct functional design; a prior apart from the Gram (`WeightedPriorGram.Unweighted`, the weighted default); weights that count rows; degree 3 and up; a ladder, SMC or a short-lived design; a model whose latent block moves the working response between sweeps (every move re-whitens the columns and discards the tables); a variable past the budget (at n = 50,000 with M_max = 20, 4 of 10 variables; at M_max = 200 and n = 10,000, 2 of 10).
- *Deviations from the brief,* each decided before the code that carries it: no hit-count gate (above); the column is still evaluated at the proposal, so the diagonal, ỹ and residual cross products and the support rule are the dots' (after lever 2 the column costs about a pass per factor and the data-only moment tables the brief sketched would not repay their size); the near-gate fallback is 10⁴ floors, not 10⁶ (the tables' error is 10⁻¹³ of the norms, which decides nothing above 10⁻⁶); an acceptance does not re-derive the border (a restore copies the exact factor, and a table entry is stored once, so the error does not accumulate); the bound is 4·K·ε + 8ε·log₂(n + 1) of ‖z‖‖u‖, the second term being the dot's own.
- *Checks.* `BassPrefixTableTests`, 19 tests:
  - the scan against the dot: 5,448 sums over n = 3…20,000, five x distributions (uniform, a point mass of 99%, ties, two end clusters, integers), three column kinds, unit and precision weights, both signs, both degrees; the worst error 92.9 ε of ‖z‖‖u‖, inside the bound everywhere;
  - a design's proposals with tables on against off, 4,000 scripted births, deaths and changes each, in five forms (unit, weighted, the ridge, categorical terms among the columns, weighted and categorical): admission equal at every proposal, ‖r‖² within 3.2·10⁻¹¹, Q within 5.5·10⁻¹¹, angles within 1.1·10⁻⁹ (relative), every Gram entry within 2·10⁻¹⁵ of the columns' norms after about 1,900 acceptances, 0–3 proposals in 4,000 re-decided by the dots (2,118–2,217 in the enumerated fixtures, whose K is 3);
  - the chain's law: the enumerated fixtures' chains (unit, weighted) with tables forced on against exact π by M's law and the ten largest states (4 chains × 30,000 sweeps), both within the existing tolerance; paired seeds, 6,000 sweeps, unit and weighted: structure sequences identical throughout, ln densities within 6.3·10⁻¹³ and 3.4·10⁻¹², τ identical; a 4-rung ladder, 500 scans: identical; **the base point, 8 chains × 20,000 sweeps, tables on against off: identical structure sequences throughout, ln densities within 5.9·10⁻¹³, 1,349,518 table entries read, 76 proposals re-decided by the dots.** The brief's check (c), two laws compared by their means, was not needed: no pair parted;
  - purity: a chain restored from a snapshot with no tables continues 800 sweeps bit for bit beside the live chain that had built tables;
  - a design not asked reads no tables (a bare design, a ladder's replica, a small fit, a wide one, a model with a moving working response); a prior apart from the Gram keeps the dots bit for bit.
- *Pinned values.* The two untempered Friedman hashes moved, as the brief expected, and are re-baselined in `BassBitIdentityTests` with the reason: equal `ee7d60a2738b56d8` → `07d9e2ca9822490c`, cancelling `414e30ea867f9058` → `9e1e7fb5e46ca40b`. With `BASS_TABLES=off` all eight former values hold on the final tree (and the tempered Friedman and the mixed and functional fixtures hold at their former values by default). Draws, not law, moved: a chain checkpointed by an older build and resumed by this one continues in the same law and not the same bits.

**Levers 2 to 5 in turn.**
- *Lever 2, the column kernels.* A term's factors run in one pass each in four lanes (degrees 1 and 2), the last also dividing by c; `Math.Max`'s zero sign and NaN are kept by a compare and a select; the library power (degree 3 and up) and the derivative route stay scalar, as does the tail. The passes `Evaluate` makes around the kernels (support count, centring, whitening, the prior side's centring, a stacked direct column's outer product) are lanes too, and `ColumnOps` carries `AggressiveOptimization`, because a vector loop in tier-0 code is far slower than the scalar loop it replaces. Re-read after lever 1: the column's share of a change falls from the report's 10.9 of 25.5 µs to about a third of that, and the cross products become the larger part of the change until lever 1 (below).
- *Lever 3, membership.* A bitset over the level codes, 64 to a word, built with the factor; `ZeroOutside` reads it once per factor. The memory is L/8 bytes, less than the sorted level list's.
- *Lever 4, the commit's residual.* Element-wise, so bit for bit: a scalar design's columns subtract in blocks of four (each entry meets its columns in column order, as before), and a direct design's inner loop over locations is in lanes. Not element-wise, and unchanged: the solve, `ResidualSquares` (a `PairwiseSum` reduction through `IInnerProduct`), and `FittedSquares`.
- *Lever 5, τ and the composition.* The split came first (`BassTauSplitTests`, one τ evaluation and the composition, each part timed apart). At the base (loaded, so read for its shares): a τ evaluation is a state (`WithLogTau`, 57–99 ns, 3–5%), an evaluation object (`Evaluate`, about half of what remains) and the marginal (`Marginal` then `LogMarginal`, a third to two-fifths of an evaluation), beside the slice's own stepping (0.3–0.4 µs of 6–16, 3%). The two `LnΓ` of the marginal (`lnΓ(g₁ + ν/2)`, `lnΓ(g₁)`) and the two of `Evaluate` (the count prior's `lnΓ(h₁ + M) − lnΓ(h₁)`) depend on neither τ nor the structure but on M, so they are kept: `Sigma2PriorTerms` (eight slots by ν, each entry the expression the marginal evaluates) and `LogCollectionWeight` by M. Each is the old expression in the old order, so the bits are the old bits. In one process, interleaved (loaded, cycle clock): the marginal 324 → 72 ns (×4.5), the collection weight 264 → 8 ns (×32). The split before and after, in two processes: the τ block 10.98 → 6.16 µs at the base, 15.92 → 6.67 at n = 200, 9.47 → 6.30 at K 5. The pre-registered threshold (a change is built only if it removes at least 25% of the block) was met by these two memos and not by the others: the state's construction and the evaluation object's allocation are 3–5% of an evaluation, the slice's stepping 3%, so **struct-constrained generics and a direct marginal path were not built**, and the runtime's dynamic PGO is left where it was. The composition's 6.7–15 µs (loaded) is spread over its block (0.4), the conditional (0.6–1.0), the draws and the interning and column-mean lookups by M; the one change it took was its sort (`Array.Sort` on the ids in place of two LINQ passes and an array), under 25% of it, so the rest was left.

**The before and after, at the cost report's points.** *No quiet window came.* Two long test hosts of other sessions (4.4 and 2.3 hours of CPU by the end of the readings), and a third for part of the time, ran from before the pass began, and the machine kept 7–10 cores busy throughout. The record's pairs (wall clock, ABBA in one quiet window) are therefore **unread**, and the commands to read them are below. What the pass has instead is **indicative**, taken on the harness's thread-cycle clock (`BASS_COST_CLOCK=cycles`, `ee36f95`): a thread starved by the load keeps its cycles, so one build can be set beside another in one window. It is not the wall clock (the cycles inflate with the cache and memory traffic of the processes around them: the base sweep reads about 180–490 µs where the cost report's quiet machine read 22), and two runs of the same binary differ by up to 2× (the base before: 487 and 253; direct 100 × 10 before: 246 and 512). Read the ratios to ±40% and the parts that move by a factor of 2 or more. Four builds, each its commit's archive with the harness over it and the kernels compiled fully optimized: before (`55e5ac0`), levers 2–4 (`a9f3211`), those with lever 1 (`4e62cf6` and its correction) and all five (`980f258`). Each point ran A, B, B, A (and, at the base and M 200, the stages in a palindrome); sweeps and parts in µs on the cycle clock, the mean of the runs, the harness's warm tables.

| point | stage | sweeps (each run) | sweep | change | death | commit (change) | τ | compose | sweep, ratio to before |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base (n 1,000, p 10, M 20, K 2) | before | 487, 253 | 370 | 413 | 57 | 135 | 25.5 | 50.1 | 1.00 |
| | levers 2–4 | 343, 272 | 308 | 340 | 49 | 49.8 | 19.4 | 31.2 | 0.83 |
| | + lever 1 | 251, 234 | 242 | 205 | 57 | 51.1 | 22.9 | 38.9 | 0.66 |
| | all five | 173, 182 | 177 | 173 | 47 | 40.9 | 10.1 | 25.6 | 0.48 |
| n 200 | before | 118, 179 | 148 | 80.6 | 20.0 | 34.6 | 24.5 | 17.9 | 1.00 |
| | all five | 115, 103 | 109 | 50.8 | 17.8 | 8.4 | 12.2 | 9.8 | 0.74 |
| M 200, n 10,000 | before | 2431, 2568 | 2500 | 3945 | 297 | 5542 | 4.5 | 102 | 1.00 |
| | levers 2–4 | 2148, 2406 | 2277 | 3868 | 305 | 2023 | 5.1 | 101 | 0.91 |
| | + lever 1 | 1829, 2485 | 2157 | 4662 | 313 | 2091 | 5.2 | 110 | 0.86 |
| | all five | 2365, 3423 | 2894 | 6198 | 441 | 2872 | 3.9 | 128 | 1.16 |
| p 200 (M 8) | before | 101, 125 | 113 | 59.8 | 5.4 | 55.4 | 7.5 | 22.3 | 1.00 |
| | all five | 143, 151 | 147 | 66.6 | 18.8 | 18.7 | 9.7 | 18.7 | 1.30 |
| L 1,024 | before | 304, 195 | 249 | 365 | 11.0 | 281 | 5.6 | 5.5 | 1.00 |
| | all five | 132, 135 | 133 | 232 | 12.3 | 84.5 | 3.8 | 8.5 | 0.53 |
| direct 100 × 10 | before | 246, 512 | 379 | 368 | 75.2 | 170 | 29.8 | 45.1 | 1.00 |
| | all five | 252, 296 | 274 | 294 | 62.1 | 176 | 11.8 | 29.3 | 0.72 |
| direct 500 × 1,000 | before | 1702, 1228 | 1465 | 1997 | 10.8 | 36857 | 4.7 | 5.6 | 1.00 |
| | all five | 764, 966 | 865 | 1271 | 7.2 | 9415 | 2.2 | 3.8 | 0.59 |

- *What the table shows.* At the base the sweep falls to about half (0.48), by stage 0.83, 0.66, 0.48; at L 1,024 to 0.53, at direct 500 × 1,000 to 0.59, at direct 100 × 10 and n 200 to 0.72 and 0.74. **Two points the noise swallows:** M 200 at n = 10,000 (DRAM-bound, 16 MB of columns per proposal, and the tables there serve 2 of 10 variables) reads 1.16 with stage ratios 0.91 and 0.86 that the spread of 1.0 between two F runs does not separate; and p 200, where tables are off by rule and levers 2–5 are the whole change, reads 1.30 with a 20% gap between the two A runs and a death proposal (which no lever touches) at 5.4 before and 18.8 after, so it is noise, not a regression found.
- *The parts that moved far beyond the spread* (the report's lever, by part): the commit (lever 4) at 0.24–0.30 of its time wherever N·M is large or the column count is (base 135 → 41, L 1,024 281 → 85, n 200 34.6 → 8.4, direct 500 × 1,000 36,857 → 9,415; M 200: 5,542 → 2,023 at levers 2–4); τ's block (lever 5) at 0.40–0.50 (base 25.5 → 10.1, n 200 24.5 → 12.2, direct 100 × 10 29.8 → 11.8, direct 500 × 1,000 4.7 → 2.2); the composition at 0.55 (base 50 → 26, n 200 17.9 → 9.8, direct 100 × 10 45 → 29), which the pass changed little, so part of it is the load; the change proposal at the base 413 → 173 (levers 2–4: 340; lever 1: 205; lever 5 has no part in it, so the last step is noise).
- *In one process, interleaved (loaded, cycle clock),* which cancels the load better than the pairs: a column's evaluation, support count and centring 12–24 → 2.7–6.6 ns per case (×2.1–6.3 over n = 200…50,000, degrees 1 and 2, orders 1 and 2); membership 2.5× at L = 2 (17.6 → 7.0 ns), 7.8–22.8× at L = 4…256 (to 1.3–1.8 ns) and 44.8× at L = 1,024 (91.8 → 2.05 ns); the commit's `Fit` ×7.0 at N = 1,000, M = 20 (1.90 → 0.27 ns per N·M), ×4.7 and ×3.0 at the M 5 and M 200 shapes, ×2.9 at n = 50,000, ×4.7 and ×4.4 for a direct design at 100 × 200 and 500 × 1,000 (3.5 → 0.76 and 3.0 → 0.69), and **×0.9 at direct 100 × 10, where the inner loop is ten locations, shorter than the lanes pay for**; the marginal ×4.5, the collection weight ×32.
- *Commands for the record's pairs, in a quiet window* (no other `testhost`, the busy cores at most 1.5 as the statistical brief states): `git archive 55e5ac0 src projects tests Directory.Build.props | tar -x -C <before>` and the same for the head (`980f258` or later) into `<after>`; `dotnet build <dir>\tests\Maths.Regression.Tests\Maths.Regression.Tests.csproj -c Release -p:ArtifactsRoot=<dir>\out\` in each; then, for each key in the points' table, A, B, B, A, each its own process: `BASS_QUALIFICATION=1 BASS_COST_BLOCK=points BASS_COST_POINTS=<key> BASS_COST_OUT=<folder> dotnet test <dll> --filter FullyQualifiedName~BassCostQualificationTests.Timed_Blocks`, below-normal priority, no `BASS_COST_CLOCK`. The head's harness warms a chain's tables; the before has none to warm.
- *Rows.* The reach, break-even and indicative rows are in [efficiency-pass-results.csv](efficiency-pass-results.csv).

**Left unbuilt, and why.**
- *Lever 1:* the column's skip until acceptance (35% in the cost report's pricing) and the data-only moment tables it needs: after lever 2 the column costs about a third of what it did, and the tables answer only the M border dots; so the diagonal, ỹ and residual dots and the support rule stay the dots'. Tables for M ≥ 100 at large n beyond the 64 MiB budget (an anchored, subsampled table would cut the memory by the sampling interval and finish a lookup over at most that many cases; M 200 at n = 10,000 needs 320 MB for all ten variables). Tables for ladders and SMC, which a hand-over or a fork would have to carry (a table is immutable once built, so sharing it by reference is possible). Interactions and categorical terms (no one-variable order serves a product). Reach between p = 10 and p = 200 was not read: 16 is a rule set below the measured failure at 200 (2.3 reads per table) and above the measured passes, not a measured edge.
- *Lever 2:* degree 3 and up (`Math.Pow`), the derivative route's branches, and the tail after the last full lane group stay scalar.
- *Lever 3:* the membership loop is a table read and a branch per case; a gather-based lane form was not tried (1.3–2 ns per case against the hinge's 3).
- *Lever 4:* a direct design with few locations (direct 100 × 10, ×0.9) needs a loop blocked over cases, which was not written.
- *Lever 5:* the closure and state allocation, the slice's stepping, and the composition's allocations, each under the 25% threshold.

**Which numbers of G-7, G-8 and the statistical half this pass would move.**
- *G-7* ([the statistical half](performance-statistical-brief.md)): its cost blocks, µs per work unit at each (fixture, backend), fall: the untempered chain by the sweep's ratio above (about half at the base, indicative), tempering and SMC by levers 2–5 only (they keep the dots); so every error-per-unit-time figure rises by the same factor, since the draws are the same (bit for bit under lever 2–5, in law under lever 1). The grid's error per unit work stands. Rerun the cost blocks only, after this pass.
- *G-8*, `BassMoments` alone: untouched; the subset factors' `Mean` and `Pair` call `Contains`, now a bit test, a change that does not show beside `PiecewisePolynomial`'s quadrature.
- *The cost report's own tables* (per-part shares, the lever table): the change proposal, the commit, τ and the composition rows; its analytic counts stand (the harness's `cross` replica still times the dots).
- *The waiting-time chip* ([its brief](waiting-time-brief.md)): its chains' law is unchanged. Draws are bit for bit under levers 2–5; under lever 1 a plain chain of a scalar fit with at most 16 continuous predictors and 500 cases or more differs in the last digits (the paired chains kept every structure over 160,000 sweeps), and `BASS_TABLES=off` restores the old bits. Artifacts it wrote with an older build resume in the same law, not the same bits. The shared notes folder carries this.

**Pre-existing failures, found by the regression runs and not mine.** `BassMixingBirthTests.BirthsAndDeaths_FromBurnedStates` (three cases: a fit's `Working` is null where the test expects the design's) and `BassQualificationRecoveryTests.Friedman_ScalarContinuous` (σ's interval [1.0127, 1.1609] misses 1) fail on the frozen baseline binary, `55e5ac0`, with the same values.

**Regression runs.** The eight pinned hashes, accelerated and with hardware intrinsics off (levers 2–4), with tables off on the final tree, and by default with the two re-baselined; `BassDesignTests`, `BassBasisTests`, `BassDesignKernelTests`, `BassRunResumeTests` on levers 2–4. On the head less its last commit (`980f258` and the tests after it), 1,228 tests of the BASS and shared-regression classes passed and 3 failed, the `BassMixingBirthTests` cases above; the run was stopped at its tail, where `BassKernelSpectrumTests`'s R̂ grid (5 minutes a case under this load) was still going, on fixtures of fewer than 500 cases that lever 1 does not touch. A first run on levers 2–4's build passed 967 tests with the same three failures and the Friedman recovery one, all four also failing on the baseline. After the last commit (`651cb90`), 220 tests passed: the prefix-table, moving-response, observation-model, weighted-core and resume classes and the eight hashes. The chip's own tests: `BassColumnBitsTests` 9, `BassMarginalTermsTests` 4, `BassPrefixTableTests` 19 (the break-even and primitive harnesses run behind `BASS_QUALIFICATION`).

**The machine, the windows, the deviations.**
- Intel Core Ultra 5 125H, 18 logical processors, Release, .NET 10.0.5, the default GC and tiered settings. Heavy runs at below-normal priority, one thread per chain; builds into `artifacts\chip-efficiency-pass\` (deleted at the end).
- Deviations from the brief, in the order met: (1) the cycle clock, a harness addition, because no quiet window came; (2) lever 1's gate and its consequences (above); (3) `NearGate` 10⁴ floors; (4) the tolerance's second term; (5) tables only on a plain chain; (6) the commit messages end with Claude Sonnet 5.5, the attribution the session's harness set (the commissioning note named Opus 5.5); (7) a stray build output at `D:\aghado01\ThermoMapperbin` and `…obj` (a script's parameter clash with `$root`) was deleted.
