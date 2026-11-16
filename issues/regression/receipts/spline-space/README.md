# Spline-space probe receipts

Tracked copies of the probe sources and result summaries behind the
[spline-space probes report](../../spline-space-probes-brief.md), at the paths they have under
`artifacts/spline-space/`.

- **Source pin.** The probe builds against the latent-recovery pilot's snapshot,
  `git archive be979dd src projects Directory.Build.props` unpacked into `artifacts/latent-recovery/source/`.
- **Build.** From `artifacts/spline-space/probe/`:
  `dotnet build -c Release -p:ArtifactsRoot=<repo>/artifacts/spline-space/build/`.
- **Run.** `SplineSpaceProbe <command> <results dir>` with `DOTNET_PROCESSOR_COUNT=6`; the process lowers its
  priority and caps its thread pool. Commands: `p1`, `p1-summary`, `p2`, `p2-summary`, `p2-diagnostic`,
  `p3`, `p3-sweep`, `p3-uniform`, `p3-summary`. Every seed derives from the gate's master seed by `SeedTree`.
- **What stays local.** Console logs, the build, `p1-designs.csv` and `p1-products.csv` (one row per design,
  1.9 MB together) and `p2-runs.csv`; the summaries hold every cell they report.
