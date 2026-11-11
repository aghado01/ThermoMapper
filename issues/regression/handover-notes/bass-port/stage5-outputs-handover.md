# Stage 5 (O1–O4) hand-over — 2026-10-29

Commits on main: 2fec1ee (O1), 29e7e23 (O3 direct), bb065c6 (O2 + O3 basis), 502fb0f (O4), dd533be (O-O1 parity), ac4730e (review).
Full report: issues/regression/bass/stage5-outputs-brief.md, "Reports".

What O7 (basis sensitivity) can read:
- `BassBasisFit` (bass/output/BassComponents.cs): `Decomposition.Map` is a `LinearReconstruction` (y_m `Offsets`, y_s `Scales`,
  K via `Basis(l, k)` / `BasisMatrix()`, `Physical(C)` = diag(y_s)·C·diag(y_s)); `Components[k]` (a `BassFit` per component,
  each with its own `Terms` table over the same predictors and knots); `ComponentData(k)`; internal `Model(k)`.
- `JointDraws.Draw(k, j)`: component k's pooled draw index in joint draw j. Joint draws weigh equally; `ChainLengths` when
  paired within chains, null after resampling.
- All components share one `BassInputTransform` (same predictors), so cross-component hinge moments read one set of knots.
- `BassFit.Predict(data, cases, options)` takes the training data, as PointwiseLogDensity does (terms index knots).
