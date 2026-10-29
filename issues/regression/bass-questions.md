Chip E reports done. I'm checking its commits and report and rerunning the suites myself, then continuing S10 with the tempered recovery, the cost baseline and the docs.

Chip E's report confirms my Friedman finding: round trips are ~100x below ideal because cold-rung local moves, not the ladder, bottleneck mixing. While the suites run, I'll check whether the sampled-g leave-one-out run is still alive.

The suites are running in my artifacts root, so the tempered Friedman run waits until they finish. Meanwhile I'm reading the existing BASS docs so the new sections match their style.

I'll draft the docs once the last numbers are in: the sampled-g leave-one-out run, the tempered Friedman run and the cost baseline. The suites are still running.

The suites still hold my artifacts root, so I'm running the tempered Friedman recovery from a third root, in parallel.

Three runs are going: the suites on chip E's work, the sampled-g leave-one-out run, and Friedman tempered and untempered. I'll carry on with the cost baseline and the docs as they finish.

Tempering roughly halves the disagreement on Friedman, though it doesn't cure it:

| | untempered | tempered (self-sized) |
|---|---|---|
| R̂ of M | 2.35 | 1.45 |
| R̂ of fitted values | 1.3–1.96 | 1.03–1.10 |
| RMSE | 0.33 | 0.28 |
| coverage | 91.6% | 97.2% |
| σ's interval | misses 1 | [0.99, 1.14], covers 1 |
| time | 0.9 s | 55 s |

The tempered run settles at 20–23 rungs, but one pair still rejects every swap, and no complete round trip happened. As chip E found, the local moves limit mixing, not the ladder. Still waiting on the suites and the leave-one-out run.

All of Chip E's test suites pass (Regression 812/812, Samplers 186/186, Distributions 129/129, LinAlg 120/120, BASS/Wallenius 37/37), though regression ran slower than before likely due to parallel load—I'll recheck that on a quiet machine. Next I'll establish the cost baseline on a single thread, with only the leave-one-out run on another core.

The forced cap on the count prior backfired—concentrating starts at M=50 makes the support rule nearly impossible to satisfy on 200 points, so the sampler correctly bails after 10⁶ draws. I'll revert to the default count prior and instead let a rich, low-noise surface naturally call for more terms, checking the resulting M against each cap.

Cost scales roughly linearly in N·M at large n (e.g., 186/465/1,184 µs for M=5/20/50 at n=5,000), with ~50–60 µs fixed overhead at small n, and the realized M stays under the cap with 10 inputs. Next I'll pull the exact text of the config and status rows to place the docs pass.

I have all results except the sampled-g leave-one-out run, so I'll proceed with the docs pass: adding four new BASS sections plus updates to the intro, API, config, status, source-map, and samplers docs.

The docs pass is in. Now the S10 report in the fitting checklist and the hub's status line, then the commits. The sampled-g leave-one-out verdict goes in once its run finishes.

Both commits are in: the qualification runs (`1e3d3f5`) and the docs with S10's report (`bdf248f`). Now the project memory.

S10 is done and chip E's self-sizing ladder is checked: both are on `main`, with one result still running.

**Chip E checked.** I read its commits (`ec565f6`…`c388461`) and reran the suites: Regression 812, Samplers 186, Distributions 129, LinAlg 120 and the BASS and Wallenius oracle tests 37, all passing. On the functional fixture where 8 fixed levels had failed, the ladder now settles at 11–16 rungs and the two forms' evidence estimates agree. The between-chain error now shows what the old within-chain error hid.

**What S10 built:**
- **The residual check (`3b5fcc0`).** You declare an index on the call: an output variable, an explicit order, or Moran's I over nearest neighbours. With several output variables the default stays null until you declare one. The posterior mean is now computed per term, not per draw.
- **Qualification runs (`1e3d3f5`)** in `tests/regression/`, behind `BASS_QUALIFICATION=1`.
- **Docs (`bdf248f`).** Four new BASS sections in [BayesianRegression.md](docs/Maths/BayesianRegression.md), chip E's shared members in [Samplers.md](docs/Maths/Samplers.md), and S10's report in the fitting checklist.

**What the qualification found:**
- **Leave-one-out check.** I scored PSIS against the exact leave-one-out density from enumerating every structure (fixed g). The retained content's k̂ can be trusted: no unit it passed was wrong. The collapsed content's k̂ is too pessimistic: it flagged 72 of 81 units, and all 72 estimates were accurate. That matches BARS. The same check under sampled g, the default, is still running, and I'll add its verdict to the report when it lands.
- **Recovery, and the main finding.** BASS mixes slowly between structures under its defaults.
  - Births and deaths are accepted about 1.5% of the time, so each chain stays in a narrow band of M. R̂ of M is 1.5–2.5 on all three surfaces.
  - R's BASS on its own Friedman sample behaves the same, so this is the method, not the port.
  - Untempered Friedman fails the σ criterion I registered: its interval is [1.013, 1.161] against a true 1. I kept the criterion as registered rather than move it.
  - Tempered, Friedman passes (σ covered, R̂ 1.45), but still completes no round trip. Chip E saw the same: the local moves, not the ladder, limit mixing.
  - The categorical and functional surfaces meet their criteria (coverage 98–99%), though their R̂ is high too.
- **Column mass vs count.** The two support rules give the same predictive spread at the edges, at both noise levels, so the count default stands.
- **Cost.** A sweep costs about linearly in N·M, from 68 µs (n = 200, M = 5) to 1,184 µs (n = 5,000, M = 50). That's roughly 6 times faster than R; my earlier "13 times" came from wall time with chains in parallel.

**Five open questions, with my recommendations:**
1. **Friedman's failed σ criterion:** keep it recorded as a failure of the untempered defaults.
2. **Slow local moves:** open this as its own design item before the efficiency pass. Improving the moves goes beyond what R does.
3. **The provisional run lengths:** leave them to the efficiency pass.
4. **Chip E's six constants:** take them as built.
5. **The ladder settling above what it needs:** it errs on the safe side, so leave it to the efficiency pass.

The full Regression suite took 5 minutes, against 1.5 earlier, while two benchmarks were also running. I'll re-time it once the machine is quiet.