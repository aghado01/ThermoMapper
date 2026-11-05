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
