# SW retention chip: hand-over note (2026-10-29)

The full report is in `D:\aghado01\ThermoMapper\issues\spc\sw-retention-chip-brief.md`, under "Reports".

What other sessions sharing this checkout should know:
- **SPCX is now v11.** `AccumulatorSerializer` refuses v8–v10 files. `SpcxDiskFrameSink.TryLoad` does not catch the refusal, so a sweep pointed at an old checkpoint directory throws `InvalidDataException`. Delete old checkpoint directories.
- **Clustering.Graphical.SPC now has `InternalsVisibleTo VizCore.Tests`.**
- **`SweepProfile.AdditionalChannels` gains `BinderCumulant`,** which adds a column to the tabular export.
- **New `AccumulationSpec` options, both off by default:** `BurnInTrace` and `ExactClusterSizeHistogram`. The CLI token parser (`SpcCommand.ParseAccumulationSpec`) doesn't know them yet.
- **Nothing in `src/maths` or `issues/regression` was touched.**
