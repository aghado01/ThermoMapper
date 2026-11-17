The coordinator has read the covariate-shift docstring change (0fede4d) and recorded it with the Lean proof it cites (deaf0da), which it verified from a clean copy. Every guarantee the conformal layer reports is now proved in Lean and verified. Nothing is in flight here.

What remains is follow-ups, all yours to pick or leave:
- **Mondrian:** run the pooled and per-stratum monitors side by side, splitting the false-alarm budget.
- **Strata in the other forms:** add strata to the full and covariate-shift forms.
- **Sleepy Jumper:** let ε drift slowly as well as jump.
- **Full form:** prune draws by weight to cut the cost on long streams.
- **Estimated likelihood ratio:** a bound on what one costs in validity.
- **The WATCH paper:** put it on the shelf so the covariate-shift method can be checked against its source.
- **Slice 15's factorization:** a stop that reads R² instead of the residual.
- **Spline density ratio:** model interactions between inputs, and decide whether `Knots` becomes the score mixture's default.