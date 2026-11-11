# LMP rank comparison: raw results (2026-11-01)

Built from commit `9abcf23` (harness `716e1f7`); the write-up is in `issues/spc-BARS/graph-compiler-brief.md`.

- `runs.tsv`, `graphs.tsv`: the raw rows. **Their headers name only the first row's (iris, one-level) columns**, so BlattHierarchy's three-level rows do not align under it. The harness was fixed afterwards to pad one schema; `analyze.py` reads these files with the per-dataset schema.
- `analysis.txt`: `analyze.py`'s tables (grid, six seeds, cold-edge rerun, compile-only readings).
