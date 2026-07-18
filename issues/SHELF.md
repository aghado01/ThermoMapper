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
