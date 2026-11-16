# Latent-recovery pilot receipts

Tracked copies of the pilot's probe sources, settings and result summaries
([brief and report](../../latent-recovery-pilot-brief.md)). Each file sits at the path it has under
`artifacts/latent-recovery/`.

- **Source pin.** The probe builds against `git archive be979dd src projects Directory.Build.props`
  unpacked into `artifacts/latent-recovery/source/`, never the shared working tree.
- **Build.** From `artifacts/latent-recovery/probe/`:
  `dotnet build -c Release -p:ArtifactsRoot=<repo>/artifacts/latent-recovery/build/`.
- **Run.** `LatentRecoveryProbe <command> <results dir> [args]` with `DOTNET_PROCESSOR_COUNT=6`; the process
  lowers its own priority to below normal and caps its thread pool. Commands: `b1`, `b2`, `b1-gate`,
  `b2-gate`, `b2-supplement`, `s1-peek`, `s1-run <k>`, `s1-gate`, `s2-peek`, `s2-run <k>`, `s2-gate`;
  `rerun` after `<k>` runs a cell's fresh replicates. Every seed derives from the settings' master seed
  by `SeedTree`, so a rerun on the pinned snapshot reproduces the rows.
- **What stays local.** Console logs, the build and the replicate rows (`*-replicates.csv`, 1.2 MB for
  Stage 1). The gate summaries (`s1-gate.txt`, `s2-readings.txt`) hold every cell; rerunning a stage with
  the tracked settings on the pinned snapshot reproduces the rows. `results/mixing/s1-mixing-summary.csv`
  keeps the post-gate mixing check's per-replicate convergence and chain spread.
