# Extension seams chip — hand-over (2026-10-31)

All seven seams landed on `main`, one commit each; the full report is prepended under "Reports" in
`issues/regression/bass/extensions-preparation-brief.md` (`530f42e`).

| Seam | Commit | Contract |
|---|---|---|
| F-X3 + O-X4 | `9b0c189` | `IDesignKernel<TPoints>` per term type (column, derivative, factors/scale); `BassBasis` and `BassMoments` compose from `DesignKernels`; `BassMoments.FactorDerivative`, `BassBasis.CaseDerivative/OutputDerivative` |
| F-X1 | `ebe2308` | `IBassObservationModel.Working` → `WorkingResponse` (y − offset, weights); `BassModel.Create(config, data, observation, form)`; the core refuses weights |
| F-X2 | `27b86b9` | `BassSweepSlot`/`BassSweepSlots` (D04's order); chain, ladder replicas and SMC kernel assemble from it; `IBassObservationModel.Blocks(slot)` |
| F-X6 | `0131c78` | `IInnerProduct` (Dot over stacked vectors, ln|Q|), `KhatriRao.Gram/Border/CrossProducts(…, output)`, `WorkingResponse.Output`; design refuses Q ≠ I |
| N-X1 | `df00275` | `KhatriRao.FaceSplittingGram/CrossProducts/HalfBandwidth`, `KroneckerGram` |
| N-X2 | `aa90524` | `WideBandCholesky`, `SpdFactors.UsesWideBand/ForBand` |
| N-X3 | `83d3d09` | `GigParameters`, GIG forms of `ConjugateUpdates.NormalVariance/NormalPrecision` |
| review / cleanup | `cbacbdb`, `d62f596` | ladder-slot test; chain and SMC descriptions |

No draw changed: the eight `BassBitIdentityTests` hashes held at every step and on the final tree.
Final cumulative suites: LinAlg 139/139, Samplers 277/277, Regression 1,383/1,383, Oracle 106/108
(the two failures are `DiagramMetricsTransportParityTests`, a 120 s R-oracle timeout that fails at
the baseline too). Baseline Oracle runs time out whenever several suites share the machine.

For the main session to verify: the working tree still holds the user's notes and the
flat-histogram chip's `BassFlatHistogramQualificationTests.cs` edit, unstaged. `artifacts\chip-ext\`
is deleted. No `docs/Maths` or checklist edits were made.
