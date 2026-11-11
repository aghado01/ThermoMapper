# F-R1 hand-over (deterministic profile)

The report is in `issues/regression/bass/deterministic-profile-brief.md`, under Reports, 2026-10-31. The per-fit results are in `deterministic-profile-results.csv`, beside the brief.

**Commits:**
- `8349663`: the harness.
- `fb35f51`: the review fixes (the CSV quoting and report wording; no reading changes).
- `c9e6191`: the report and the CSV.

**The verdict: no arm clears.** All 144 fits read a predictive R̂ above 1.05 (min 1.23, median 2.14 noise-free and 1.74 noisy). The least ESS is about 5 of 8,000 draws, and D10's verdict is N everywhere. So the gate leaves every benchmark unscored, and F-R1 closes as read. F-M2 read the same mixing at S10's lengths.

**Context only, not a result:** without the gate, A6 (ridge, τ sampled, λ (1, 1)) would score 0.598 and A4 (g, Zellner–Siow, λ (1, 1)) 0.825, and both would clear. A6's lead rests most on the functional benchmark's design 0. On the noisy controls, ungated, no arm reaches 0.9.

The files in this folder:
- `fr1-deterministic-profile-report.md`: the harness's own report.
- `fr1-deterministic-profile-results.csv`: the same CSV, with the verdict re-quoted.
- `fr1_rule_check.py`: the independent recomputation of the rule. Run it as `python -I fr1_rule_check.py <csv> 3`.

`artifacts\chip-fr1\` and the run directories are deleted.
