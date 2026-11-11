# Conformal Lean chip: hand-over, 2026-11-10

Landed: slice 10 grid reading (bb5df90), slice 13 (f998131), slice 14 (27a6662; 340 declarations, standard axioms, full checks green in 1,820 s).

Ready in scratch, not spliced. Scratch = C:\Users\azrie\AppData\Local\Temp\claude\6b380889b4cf6b3a\93f88779-bc7a-494c-a9d6-45d614c0ca5f\scratchpad
- Slice 15 (PCP):
  - `s15.lean` compiles clean against HEAD; 9 declarations. `s15b.lean` is an older draft; ignore it.
  - `s15_block.py`: the checks pass. It needs slice 14's `w_refused` and CB-4's `prod`.
  - `s15_ledger.md`, `s15_report.md`.
- Slice 16 (batches):
  - `s16.lean` compiles clean (via `S16Test.lean`). It is two-sided:
    - `batch_interval`: both sides on one rank event;
    - `batch_coverage`: the 1 − γ fraction;
    - `card_rankSeqs`: |H| = C(n+m, m);
    - `cal_lt_orderEmbOfFin` / `orderEmbOfFin_lt_cal`: the order-statistic sandwich.
  - `s16_block.py` was updated for the lower bound but not rerun since that edit.
  - No ledger section or report has been drafted yet.

Per item:
1. Splice the Lean before `end Conformal`.
2. Splice the checks before the final `print()` / `print("All assertions hold.")`.
3. Append the ledger section, then update the status line: count with the broad declaration regex, 340 + new.
4. Prepend the report under "## Reports\n\nNewest first.\n\n".
5. Gates: module and barrel build (`run-lake.ps1`), axioms probe (pattern in `AxiomsProbe10.lean`), `lean/scripts/meta-ci.ps1 -NoBuild`, D: free space, full checks in the background (~30 min, below-normal priority).
6. Commit the four files by pathspec.

Slice 18 (Snell & Griffiths 2502.13228; Gressmann et al. 1801.00753) not started.
