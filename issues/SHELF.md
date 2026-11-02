# The Shelf

A place to jot down action items for later when they come up organically in the course of doing other work

## Dummy Example

Each entry can house a given concept and may be appended or amended over time as an idea incubates before it is promoted

## RNG one-die migration (from the seed-derivation audit, 2026-07-17)

The seed-derivation aliasing audit (see `project-rng-hygiene` memory) fixed the three arithmetic
derivation sites (DistributedSpred `+1009·block`, BicKSweep `+7919·restart`, ChainEnsemble additive
readout salt — all now SeedTree/SplitMix64). What remains is the **one-die** half of the theme:
`System.Random` still serves as the die in the GMM stack (`GaussianMixtureModel.RandomInitialize`,
`BicKSweep`), `kmeans`, `MapperGMM`, `ICA`, `MatrixOps`, and the synthetic generators. Migrating
them to `Xoshiro256PlusPlus` is API churn (signatures take `Random`) with cross-version stream
stability as the payoff — `System.Random`'s algorithm is not pinned across .NET versions, Xoshiro's
is ours. Do it as one deliberate pass, regenerating any pinned fixture expectations.

## EyeTorusToy noise-stream derivation (audit note, 2026-07-17)

`EyeTorusToy.cs` derives its noise stream as `seed ^ 0x9E3779B9` — the XOR cousin of the fixed
additive aliasing (main stream of seed `s ^ C` ≡ noise stream of seed `s`). Benign in practice:
it's a fixture generator and the two streams serve different draws. Fold a `SeedTree.Derive(seed, 2)`
pair into the next moment the fixture's generated clouds are allowed to churn (any pinned
expectations regenerate with it); not worth a churn cycle on its own.

## Left over from closing the PH loose ends (2026-11-02)

None of these blocks anything. Context: `issues/ph/ph-loose-ends-brief.md`, its report.

- **Compatible cycles along a zigzag bar.** `Bar.Cycle` names a bar's class where it begins. Under `Flag` the chain of a hole
  that opens never fits its first frame, and a cycle can move before a bar ends. Dey, Hou and Morozov (arXiv:2410.20565)
  compute the sequence; the coned run's V may give it directly (derived, unverified). Pick up when a consumer asks for the
  cycle at a frame.
- **Involuted representatives on negative values.** `PersistentInvolutedHomology` with `representatives` hangs on an essential
  H₁ bar when an edge value is negative: `CycleReconstruction` runs a shortest path weighted by those values. Rips never gets
  there. Pick up before anything non-Rips asks that engine for cycles.
- **The relative reading's timing on a quiet machine.** Today's ran at 80–100% load, so `ComputeBothReadings` against `ComputeBarcode` has a sign at
  42,000 triangles and no reliable size, and no sign at n≤200. Set `RELATIVE_READING_TIMING` to an output file and run `BothReadings_Timing` in Release,
  about 25 minutes under load at the default three passes. Pick up on a quiet machine, or if a caller wants both readings and the internal method's cost matters.
- **`PersistentCohomology` on a `CellFiltration`.** Throws, since its H₀ path reads `GetVertices`. Pick up if a cell complex
  ever wants the clearing path.
