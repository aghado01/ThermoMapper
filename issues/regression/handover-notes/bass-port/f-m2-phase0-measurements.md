# F-M2 Phase 0: the measurements (chip report pointer)

Date: 2026-10-28. The full report is in `issues/regression/bass/local-move-measurements-brief.md`, section "Reports" ("Phase 0 measured: M1–M5").

Commits on `main`:
- `fc88448`: the helper and M2.
- `6b365ff`: M3. `98381db`: M4. `11f1204`: M5. `a13eb80`: M1.
- `7f76190`: the report.

There is no production change and no internal accessor. Regression is 828/828 with the switch off.

The readings, in short (detail and strength in the report):
- *M1:* at Thin = 1 the composition is 34–94% of a sweep, the structure move 6–23% and τ 0.1–30%. A second structure move per sweep adds only the structure share. The machine was loaded, so the shares hold to about ±20% of a sweep.
- *M2:* births fail on the Occam cost (Δℓ ≈ −½ ln N against a threshold of about 2 nats), not on the inert inputs. Active and inert main effects accept alike (1.6% against 1.4%), and lever 3a moves acceptance by −7% to +24%. Every held term is load-bearing: deaths cost tens to hundreds of nats.
- *M3:* a kept-sign knot step of 1–3 values accepts 0.93–0.97 on Friedman and the mixed fit, against the redraw's 0.04. The functional fit is much steeper: 0.62 and 0.37 at one step.
- *M4:* the untempered chains sit in different modes 12–42 nats apart in median log posterior (Friedman, functional). On the mixed fit, one chain sits 22 nats below the other three.
- *M5:* at Δb = 0.01 the steep band's pairs reject at most 0.28–0.46, so it is a finite barrier that the sizing under-reads. With the band opened, the ideal round trips reach 0.03 per scan, yet none is observed: the local moves bound the run. The Gaussian Σh over-reads the band.
