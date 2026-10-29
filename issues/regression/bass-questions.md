Gathered the hub and checklists; now checking your three untracked BASS question files for open items.

BASS's fitting, outputs and presentation data are built and verified. The local-move comparison pass is running, and calibration (stage 7) and SMC (stage 8) haven't started. I went through the hub, the five checklists, the segment briefs, the F-M2 records and your three untracked question files. Your files hold nothing open: argon's notes were all handled in the 2026-10-28 circle-back, and the S10 questions were decided.

## Where BASS stands

| Stage | Status |
|---|---|
| 0–3: baseline, spec, numerics, sampler infrastructure | Done. Numerics checklist 14/17; the 3 open items are extension prep. |
| 4: fitting (MCMC, tempering, restart) | Done and qualified by S10. Fitting checklist 28/37; open: F-M2, F-M1 (now inside F-M2's pass), F-R1, and seams F-X1–X6. |
| F-M2: local-move mixing | Phase 1 options landed; **Phase 2 pass running**. |
| 5–6: basis/PCA, prediction, sensitivity, effects | Done. Outputs checklist 19/20; only seam O-X4 is open. |
| 7: calibration | **Not started.** It waits on D09, which is still a draft. Calibration checklist 2/10. |
| 8: SMC (fitting adapter, then cut-calibration SMC) | **Not started.** |
| 9: presentation | The data is done. Rendering waits for a consumer above maths and for viz v2's contracts. |
| Extensions E1–E13 | All candidates; none scheduled. |

## Open questions

**Waiting on your decision**
1. **D09, the calibration design.** It blocks stage 7. The draft's choices:
   - the cut posterior as the target, against plug-in or fully joint calibration;
   - the functional and multivariate covariance formulations as declared types;
   - variance-group priors, defaulting to R's IG(ν/2, ν·ŝ²/2), with R's heavy-tailed kernel chosen per group;
   - **the discrepancy's prior.** R fixes its precision at 1; the draft says to declare it but names no default.
   - categorical calibration parameters, which R's API lacks;
   - a nested reference sampler that measures the error from finite inner runs.
2. **D10's stopping semantics.** Still open. It decides what "converged" means beyond today's opt-in R̂ stop, and it gates the run-length defaults.
3. **F-R1, the deterministic-simulator profile.** Drafted 2026-10-26. Its coefficient prior (ridge against g-prior for noise-free simulators) is to be set by a pre-registered comparison that hasn't been designed.

**Answered by the pass, or following it**

4. **Which mixing option becomes the default.** This is decided by the fixed rule. I'll verify the verdict and carry it out.
5. **Run lengths.** They're still R's provisional 9,000 burn-in + 1,000 draws, to be sized from the winning option's ESS per sweep.
6. **The ladder's sizing.** On the functional fixture it settles above need, and on Friedman a steep band remains.
7. **Friedman's worse modes.** No option opens them, and neither does tempering. If the pass confirms this, the candidates are:
   - several moves per swap on tempered runs (lever 1, kept as a knob);
   - the variational reference for tempering;
   - a compound move that replaces a load-bearing term at once.
8. **The efficiency pass,** the cost work held until F-M2 closes.

**Defaults chips set as reversible, never put to you**

9. A choice of mine: the factor moves share the change's slot in thirds.
10. From stage 6: the output mean y_m is placed in the location effect for basis sensitivity.
11. From stage 9:
    - Sobol indices are summarized over the defined draws;
    - `NoiseSd` = E[σ];
    - the trace axis counts a ladder's scan as one sweep;
    - basis training-fit intervals take one predictive draw per joint draw.

**For stage 8**

12. The shared weighted quantile jumps between exactly equal and nearly equal weights. Once SMC brings unequal weights, a weighted type 7 through Kish's effective sample size is the fix.
13. OQ2: an exact reference rung on long ladders is still open. The 4- and 8-rung study showed no payoff.

**Parked prep and extensions**

14. **Seams:** F-X1–X6, C-X1/X2, N-X1–X3, O-X4. C-X1 (a declared measurement-error model) and C-X2 (emulator draws behind an interface) fit naturally into stage 7 itself.
15. **Extensions E1–E13:** robust noise, ABMRS, KHAOS, active subspaces, cross-engine evidence and the rest. None has a trigger yet.

**Stale labels in the records**

16. Two section headers still read "Open (2026-10-23)" although the work has landed:
    - OQ17 in the fitting brief;
    - OQ18 in the numerics brief.
17. The hub's decision table is headed "Status (2026-10-27)" but carries updates through 10-29.

My recommendation: take D09 next while the pass runs, since it unblocks stage 7 and doesn't touch the pass's files. D10 can follow once the pass's ESS numbers give the stopping rule something concrete to target.