# Grassmann tight-cluster probe (2026-11-08, Grassmann discussion session)

Probe of `GeometricMedian` / `GeometricMean` (IRLS, `IrlsOptions.Default`) on Gr(2,10) and Gr(3,20) clusters with spreads 1e-2..1e-9, against truth from outside the solver, with a Euclidean control and a probe-only sine-based Grassmann. Results and conclusions: `issues/regression/discussion/grassmann-stiefel-subspaces-brief.md`, Report, "probe" entry.

To rerun: `git archive <rev> Directory.Build.props projects/Maths.Geometry projects/Maths.LinAlg projects/Maths.Distance projects/Maths.Distributions projects/Maths.Rng src/maths | tar -x -C <dir>/snap`, put these files in `<dir>/probe`, `dotnet build -c Release` the csproj (it references `..\snap\projects\Maths.Geometry`), run the exe. `run.log` is the run at d571dff (2 s).
