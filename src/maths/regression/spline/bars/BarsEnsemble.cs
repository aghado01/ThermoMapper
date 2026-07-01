using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The curve's peak (argmax) as a posterior, not a point estimate — the downstream T_c readout. Each draw
/// contributes its exact closed-form peak (reduce-per-draw, never argmax-of-the-mean), pooled into a location
/// with a 95% credible interval and its own R̂/ESS.
/// </summary>
public sealed record PeakPosterior(
    double LocationMean, double LocationLo, double LocationHi,
    double HeightMean, double LocationRHat, double LocationEss);

/// <summary>
/// One transition resolved from the pooled fields — the per-mode generalization of <see cref="PeakPosterior"/>
/// (the dominant mode, k = 1). <see cref="Location"/> and its credible interval come from the peak-intensity λ(T)
/// mass over the mode's support (epistemic — where the transition is); <see cref="SpanLeft"/>/<see cref="SpanRight"/>
/// come from the coverage π(T) (structural — how wide the feature is); <see cref="Mass"/> = ∫λ over the support
/// (≈ posterior probability of a transition there). <see cref="LocationRHat"/> is the cross-chain R̂ on the mode
/// location, computed matching-free from the per-chain λ histograms (each chain's λ over the support IS its own
/// peak-location distribution there) — ≈1 when the chains agree, &gt;1 a disagreement / multimodality flag, NaN if
/// fewer than two chains resolved the mode. Modes are ordered by mass, dominant first.
/// </summary>
public sealed record PeakMode(
    double Location, double LocationLo, double LocationHi,
    double SpanLeft, double SpanRight, double Mass, double LocationRHat);

/// <summary>
/// The pooled posterior fit of a BARS run plus its convergence diagnostics, all on the reported functionals
/// (knot count, fitted grid values) — the quantities you read, monitored for convergence.
/// <see cref="RHatFit"/> is the per-grid-point consensus map: where it is elevated the chains disagree on the
/// fit (multimodality / under-determination), where it is ≈1 the fit is settled.
/// <see cref="PeakIntensity"/> is the grid-aligned peak-intensity λ(T): the matching-free posterior density
/// of significant transitions (∑ = <see cref="PeakCountMean"/>), its modes the persistent transitions.
/// <see cref="SpanCoverage"/> is the grid-aligned span-coverage π(T): the posterior probability each T lies
/// within a peak's FWHM-esque span — its plateaus the structural feature widths, its complement the stable regimes.
/// <see cref="PeakModes"/> resolves each transition from those fields — location + credible interval from λ, span
/// from π — with <see cref="Peak"/> the dominant mode (the k = 1 case).
/// </summary>
public sealed record BarsResult(
    double[] Grid,
    double[] Fit,
    double[] RHatFit,
    double MeanKnots,
    double PeakCountMean,
    double[] PeakIntensity,
    double[] SpanCoverage,
    double RHatKnots,
    double RHatFitMax,
    double EssKnots,
    PeakPosterior Peak,
    IReadOnlyList<PeakMode> PeakModes,
    int SamplesUsed,
    double AcceptanceRate,
    int[] ChainSeeds,
    IReadOnlyList<MoveStat> MoveStats);

/// <summary>
/// The BARS client of the kernel-agnostic ensemble shell (<see cref="ChainEnsemble"/>): it knows how to start one
/// reversible-jump replica (overdispersed start, optional adaptive-τ kernel or parallel-tempering ladder, optional
/// scale-mixture weight resampler) and how to reduce each draw to the readout functionals (knot count + fitted
/// grid values). The shell owns the fan-out, the burn → rounds → R̂ early-stop, and the pooled R̂/ESS; this owns
/// the BARS semantics and assembles <see cref="BarsResult"/> from the pooled functionals + its own pooled sinks
/// (peak-intensity λ(T), span-coverage π(T), the per-draw argmax cloud). Parallelism is across chains, never data.
/// </summary>
public sealed class BarsEnsemble : IEnsembleModel<KnotConfig>
{
    private readonly IBasis _basis;
    private readonly IObservationModel _model;
    private readonly IComplexityPrior _prior;
    private readonly IKnotKernel _kernel;
    private readonly IWeightResampler? _resampler;
    private readonly double? _dmgkConstant;     // null ⇒ equal birth/death/relocate weights
    private readonly double? _adaptKernelTau;    // set ⇒ each chain tunes an AdaptiveLocalBetaKernel in burn-in
    private readonly int? _temperLevels;        // set ⇒ each chain is a parallel-tempering ladder (cold = the sample)
    private readonly double _temperBetaMin;

    // Per-run context, set at the top of Run and read (only) by StartChain/the chain runs during that call. Run
    // drives the shell synchronously, so these are valid for the whole run; the fan-out reads them, never writes.
    private double[] _x = null!;
    private double[] _y = null!;
    private double[] _grid = null!;
    private int _weightEvery;
    private int _startDispersion;
    private double _peakProminence;
    private double _spanDropFraction;

    /// <param name="dmgkConstant">When set (DMGK uses 0.4), birth/death/relocate use the prior-aware DMGK
    /// schedule (<see cref="DmgkSchedule"/>) instead of equal weights — better dimension mixing.</param>
    /// <param name="adaptKernelTau">When set, each chain proposes with an <see cref="AdaptiveLocalBetaKernel"/>
    /// started at this τ and tunes it to the relocate target during burn-in (then freezes) — robust to a
    /// mis-scaled proposal. Overrides <paramref name="kernel"/> for proposals.</param>
    /// <param name="temperLevels">When set (≥ 2), each "chain" becomes a <see cref="ParallelTempering{T}"/> ladder of
    /// this many replicas, and the cold (β = 1) replica supplies the draws — so the readout is unchanged but the
    /// chain can cross modal barriers (a multi-peak T_c posterior). Mutually exclusive with the resampler and with
    /// adaptive-τ in this version (a plain fixed-kernel tempered ladder).</param>
    /// <param name="temperBetaMin">The hot end of the geometric β-ladder (default 0.05).</param>
    public BarsEnsemble(IBasis basis, IObservationModel model, IComplexityPrior prior, IKnotKernel kernel,
                        IWeightResampler? resampler = null, double? dmgkConstant = null, double? adaptKernelTau = null,
                        int? temperLevels = null, double temperBetaMin = 0.05)
    {
        ArgumentNullException.ThrowIfNull(basis);
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(prior);
        ArgumentNullException.ThrowIfNull(kernel);
        if (temperLevels is int lv)
        {
            if (lv < 2) throw new ArgumentOutOfRangeException(nameof(temperLevels), "Tempering needs ≥ 2 levels.");
            if (resampler is not null || adaptKernelTau is not null)
                throw new ArgumentException("Tempering is mutually exclusive with the resampler and adaptive-τ.", nameof(temperLevels));
        }
        _basis = basis;
        _model = model;
        _prior = prior;
        _kernel = kernel;
        _resampler = resampler;
        _dmgkConstant = dmgkConstant;
        _adaptKernelTau = adaptKernelTau;
        _temperLevels = temperLevels;
        _temperBetaMin = temperBetaMin;
    }

    /// <summary>Functionals R̂ runs on: knot count (index 0) + each fitted grid value (the per-grid consensus map).</summary>
    public int FunctionalDim => 1 + _grid.Length;

    /// <summary>Only the knot count (index 0) gets ESS — the per-grid sequences would be a needless O(n²)·g blow-up.</summary>
    public int EssDim => 1;

    /// <summary>
    /// Start one reversible-jump replica from the shell's seeds: <paramref name="chainSeed"/> drives the
    /// overdispersed start, the chain (or the tempering sub-ladder, which re-derives from this same seed), and the
    /// resampler RNG; <paramref name="readoutSeed"/> drives the per-sample coefficient draws, decoupled from the chain.
    /// </summary>
    public IChainRun<KnotConfig> StartChain(int chainSeed, int readoutSeed)
    {
        var rng = new Xoshiro256PlusPlus(chainSeed);
        var run = new BarsChainRun(this, rng, OverdispersedStart(rng, _startDispersion), _grid.Length)
        {
            ReadoutRng = new Xoshiro256PlusPlus(readoutSeed),
            Kernel = _adaptKernelTau is double t0 ? new AdaptiveLocalBetaKernel(t0) : _kernel,
        };
        if (_temperLevels is int levels)
        {
            run.Ladder = NewLadder(run, levels, chainSeed);
            run.SetChain(run.Ladder.ColdChain);   // the cold (β = 1) replica is the sampler; its draws are the posterior
        }
        else
        {
            run.SetChain(NewChain(run));
        }
        return run;
    }

    /// <remarks>
    /// <paramref name="rHatTarget"/> &gt; 0 stops sampling early once R̂ on knot count and every grid point is ≤
    /// it (checked each round of <paramref name="batchSize"/> samples; a non-positive batch is a single round of
    /// <paramref name="samples"/>).
    /// </remarks>
    public BarsResult Run(double[] x, double[] y, double[] grid,
                          int chains = 4, int masterSeed = 0, int burn = 2000, int samples = 2000,
                          int startDispersion = 4, int weightEvery = 25, double rHatTarget = 0.0, int batchSize = 0,
                          double peakProminence = 0.1, double spanDropFraction = 0.5)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(grid);
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 1) throw new ArgumentOutOfRangeException(nameof(samples));

        _x = x;
        _y = y;
        _grid = grid;
        _weightEvery = weightEvery;
        _startDispersion = startDispersion;
        _peakProminence = peakProminence;
        _spanDropFraction = spanDropFraction;

        (EnsembleRun summary, IReadOnlyList<IChainRun<KnotConfig>> chainRuns) =
            new ChainEnsemble().Run<KnotConfig>(this, chains, masterSeed, burn, samples, rHatTarget, batchSize);

        int g = grid.Length;
        int collected = summary.SamplesUsed;
        double denom = chains * (double)collected;

        var bars = new BarsChainRun[chains];
        for (int ci = 0; ci < chains; ci++) bars[ci] = (BarsChainRun)chainRuns[ci];

        // Pooled functionals: index 0 = knot count, indices 1..g = fitted grid values.
        var fit = new double[g];
        for (int gi = 0; gi < g; gi++) fit[gi] = summary.FunctionalMean[1 + gi];
        double meanK = summary.FunctionalMean[0];

        // Consensus map: clamp a degenerate-variance grid point's R̂ to 0 (the historical FitRHat display rule);
        // the knot R̂ is reported raw, as before.
        var rHatFit = new double[g];
        double rHatFitMax = 0.0;
        for (int gi = 0; gi < g; gi++)
        {
            double r = summary.FunctionalRHat[1 + gi];
            double rv = double.IsNaN(r) || double.IsInfinity(r) ? 0.0 : r;
            rHatFit[gi] = rv;
            if (rv > rHatFitMax) rHatFitMax = rv;
        }
        double rHatK = summary.FunctionalRHat[0];
        double essK = summary.FunctionalEss[0];

        // BARS's own pooled sinks (the orchestrator does not know what these mean).
        var peakIntensity = new double[g];
        var spanCoverage = new double[g];
        for (int gi = 0; gi < g; gi++)
        {
            double si = 0.0, sc = 0.0;
            for (int ci = 0; ci < chains; ci++) { si += bars[ci].PeakIntensitySum[gi]; sc += bars[ci].SpanCoverageSum[gi]; }
            peakIntensity[gi] = si / denom;
            spanCoverage[gi] = sc / denom;
        }

        double peakCountMean = 0.0;
        for (int ci = 0; ci < chains; ci++) peakCountMean += bars[ci].PeakCountSum;
        peakCountMean /= denom;

        var moveAtt = new Dictionary<string, long>();
        var moveAcc = new Dictionary<string, long>();
        for (int ci = 0; ci < chains; ci++)
        {
            BankMoveStats(bars[ci]);   // fold the final (post-last-rebuild) segment
            foreach (var kv in bars[ci].MoveAtt) moveAtt[kv.Key] = moveAtt.GetValueOrDefault(kv.Key) + kv.Value;
            foreach (var kv in bars[ci].MoveAcc) moveAcc[kv.Key] = moveAcc.GetValueOrDefault(kv.Key) + kv.Value;
        }
        var moveStats = new List<MoveStat>(moveAtt.Count);
        foreach (string key in moveAtt.Keys)
            moveStats.Add(new MoveStat(key, moveAtt[key], moveAcc.GetValueOrDefault(key)));

        PeakPosterior peak = SummarisePeak(bars, collected, denom);
        var chainLambda = new double[chains][];
        for (int ci = 0; ci < chains; ci++) chainLambda[ci] = bars[ci].PeakIntensitySum;
        var peakModes = ExtractModes(grid, peakIntensity, spanCoverage, chainLambda);
        return new BarsResult(grid, fit, rHatFit, meanK, peakCountMean, peakIntensity, spanCoverage, rHatK, rHatFitMax, essK, peak, peakModes, collected, summary.AcceptanceRate, summary.ChainSeeds, moveStats);
    }

    private static PeakPosterior SummarisePeak(BarsChainRun[] states, int n, double denom)
    {
        int c = states.Length;
        var sums = new double[c];
        var sqs = new double[c];
        var peakChains = new double[c][];
        var all = new List<double>();
        double mean = 0.0, height = 0.0;
        for (int i = 0; i < c; i++)
        {
            sums[i] = states[i].PeakLocSum;
            sqs[i] = states[i].PeakLocSumSq;
            peakChains[i] = states[i].PeakDraws.ToArray();
            all.AddRange(states[i].PeakDraws);
            mean += states[i].PeakLocSum;
            height += states[i].PeakHeightSum;
        }
        mean /= denom;
        height /= denom;
        all.Sort();
        double lo = Quantile(all, 0.025);
        double hi = Quantile(all, 0.975);
        double rHat = ChainDiagnostics.RHat(sums, sqs, n);
        double ess = ChainDiagnostics.Ess(peakChains);
        return new PeakPosterior(mean, lo, hi, height, rHat, ess);
    }

    private static double Quantile(List<double> sorted, double q)
    {
        int n = sorted.Count;
        if (n == 0) return double.NaN;
        if (n == 1) return sorted[0];
        double pos = q * (n - 1);
        int lo = (int)Math.Floor(pos);
        double frac = pos - lo;
        return lo + 1 < n ? sorted[lo] * (1.0 - frac) + sorted[lo + 1] * frac : sorted[lo];
    }

    // Advance the chain one step (the IChain.Step body): step the lone chain or the tempered ladder (returning the
    // cold state), then run the optional scale-mixture weight resampler — every weightEvery steps it redraws w and
    // rebuilds the RJ knot chain | w (banking the retired chain's acceptance + per-move counts first).
    private KnotConfig StepOnce(BarsChainRun st)
    {
        st.Config = st.Ladder is not null ? st.Ladder.Step() : st.RjChain.Step();

        if (_resampler is not null)
        {
            st.SinceWeight++;
            if (st.SinceWeight >= _weightEvery)
            {
                st.W = _resampler.Resample(st.Config, st.W, _x, _y, _basis, _model, st.Rng);
                st.BankedAccepted += st.RjChain.Accepted;
                st.BankedAttempts += st.RjChain.Attempts;
                BankMoveStats(st);
                st.SetChain(NewChain(st));
                st.SinceWeight = 0;
            }
        }
        return st.Config;
    }

    // Reduce one post-burn draw: fit the grid values into the functionals (index 0 = knot count), then deposit the
    // BARS-specific sinks — the per-draw closed-form argmax cloud, and per significant peak its location into the
    // grid-aligned intensity accumulator (pooled → λ(T), ∑_grid λ = PeakCountMean) and its FWHM-esque span into the
    // coverage accumulator (pooled → π(T)).
    private void AccumulateDraw(BarsChainRun st, KnotConfig draw, Span<double> functionals)
    {
        int g = _grid.Length;
        double[,] zTrain = _basis.Design(draw, _x);
        double[] coef = _model.SampleCoefficients(zTrain, _y, st.W, st.ReadoutRng);
        double[,] zGrid = _basis.Design(draw, _grid);
        for (int gi = 0; gi < g; gi++)
        {
            double f = 0.0;
            for (int j = 0; j < coef.Length; j++) f += zGrid[gi, j] * coef[j];
            functionals[1 + gi] = f;
        }
        functionals[0] = draw.Count;

        if (_basis is SplineBasis splineBasis)
        {
            (double peakX, double peakF) = SplineExtrema.Argmax(draw, coef, splineBasis);
            st.PeakLocSum += peakX;
            st.PeakLocSumSq += peakX * peakX;
            st.PeakHeightSum += peakF;
            st.PeakDraws.Add(peakX);

            var spans = SplineExtrema.SignificantPeakSpans(draw, coef, splineBasis, _peakProminence, _spanDropFraction);
            st.PeakCountSum += spans.Count;
            for (int p = 0; p < spans.Count; p++)
            {
                PeakSpan sp = spans[p];
                st.PeakIntensitySum[NearestGridIndex(_grid, sp.Location)] += 1.0;
                for (int gi = 0; gi < g; gi++)
                    if (_grid[gi] >= sp.Left && _grid[gi] <= sp.Right) st.SpanCoverageSum[gi] += 1.0;
            }
        }
    }

    private ReversibleJumpChain<KnotConfig> NewChain(BarsChainRun st)
        => new(Moves(st.Kernel), new SplineTarget(_basis, _model, _prior, _x, _y, st.W), st.Config, st.Rng);

    private ParallelTempering<KnotConfig> NewLadder(BarsChainRun st, int levels, int seed)
    {
        var target = new SplineTarget(_basis, _model, _prior, _x, _y, st.W);
        double[] betas = ParallelTempering<KnotConfig>.GeometricLadder(levels, _temperBetaMin);
        return new ParallelTempering<KnotConfig>(Moves(st.Kernel), target, st.Config, betas, seed);
    }

    private IRjMove<KnotConfig>[] Moves(IKnotKernel kernel)
    {
        var birth = new KnotBirthMove(kernel);
        var death = new KnotDeathMove(kernel);
        var relocate = new KnotRelocateMove(kernel);
        return _dmgkConstant is double c
            ? DmgkSchedule.Wrap(birth, death, relocate, _prior, c)
            : new IRjMove<KnotConfig>[] { birth, death, relocate };
    }

    // Burn-in with optional adaptive proposal-scale tuning: when the chain's kernel is adaptive, nudge τ toward
    // the relocate target every adaptEvery steps (vanishing Robbins–Monro step), then leave it frozen for the
    // sampling phase (the diminishing-adaptation condition for ergodicity).
    private void AdaptiveBurn(BarsChainRun st, int burn)
    {
        if (st.Kernel is not AdaptiveLocalBetaKernel ak)
        {
            for (int k = 0; k < burn; k++) StepOnce(st);
            return;
        }

        const int adaptEvery = 200;
        long prevAtt = 0, prevAcc = 0;
        int done = 0, round = 0;
        while (done < burn)
        {
            int step = Math.Min(adaptEvery, burn - done);
            for (int k = 0; k < step; k++) StepOnce(st);
            done += step;
            (long att, long acc) = RelocateStat(st.RjChain);
            long dA = att - prevAtt, dC = acc - prevAcc;   // dA ≤ 0 ⇒ chain rebuilt (resampler) → skip this round
            prevAtt = att; prevAcc = acc;
            if (dA > 0) ak.Adapt((double)dC / dA, stepSize: 3.0 / Math.Sqrt(round + 1));
            round++;
        }
    }

    private static (long Attempts, long Accepted) RelocateStat(ReversibleJumpChain<KnotConfig> chain)
    {
        foreach (MoveStat m in chain.MoveStats())
            if (m.Key == "relocate") return (m.Attempts, m.Accepted);
        return (0, 0);
    }

    // Fold the chain's current per-move counts into the chain-state banks (the chain resets on resampler rebuilds).
    private static void BankMoveStats(BarsChainRun st)
    {
        foreach (MoveStat m in st.RjChain.MoveStats())
        {
            st.MoveAtt[m.Key] = st.MoveAtt.GetValueOrDefault(m.Key) + m.Attempts;
            st.MoveAcc[m.Key] = st.MoveAcc.GetValueOrDefault(m.Key) + m.Accepted;
        }
    }

    // Nearest grid index to x for the peak-intensity deposit (grid ascending — it is the evaluation grid).
    private static int NearestGridIndex(double[] grid, double x)
    {
        int idx = Array.BinarySearch(grid, x);
        if (idx >= 0) return idx;
        int ins = ~idx;                              // first index with grid[index] > x
        if (ins <= 0) return 0;
        if (ins >= grid.Length) return grid.Length - 1;
        return (x - grid[ins - 1]) <= (grid[ins] - x) ? ins - 1 : ins;
    }

    // Resolve transitions from the pooled fields: each local maximum of λ(T) is a mode; its location + credible
    // interval come from λ's mass over the mode's support, its span from the π ≥ ½ region around it. The dominant
    // (highest-mass) mode is the k = 1 case carried by PeakPosterior; this generalizes it to every transition.
    private static List<PeakMode> ExtractModes(double[] grid, double[] lambda, double[] pi, double[][] chainLambda)
    {
        int g = grid.Length;
        var modes = new List<PeakMode>();
        for (int i = 1; i < g - 1; i++)
        {
            if (lambda[i] <= 0.0 || !(lambda[i] > lambda[i - 1] && lambda[i] >= lambda[i + 1])) continue;

            int lo = i;
            while (lo > 0 && lambda[lo - 1] > 0.0 && lambda[lo - 1] <= lambda[lo]) lo--;
            int hi = i;
            while (hi < g - 1 && lambda[hi + 1] > 0.0 && lambda[hi + 1] <= lambda[hi]) hi++;

            double mass = 0.0, wsum = 0.0;
            for (int j = lo; j <= hi; j++) { mass += lambda[j]; wsum += lambda[j] * grid[j]; }
            if (mass <= 0.0) continue;

            double location = wsum / mass;
            double loCi = MassQuantile(grid, lambda, lo, hi, mass, 0.025);
            double hiCi = MassQuantile(grid, lambda, lo, hi, mass, 0.975);

            int sl = i, sr = i;
            while (sl > 0 && pi[sl - 1] >= 0.5) sl--;
            while (sr < g - 1 && pi[sr + 1] >= 0.5) sr++;
            bool covered = pi[i] >= 0.5;

            double rhat = ModeRHat(chainLambda, grid, lo, hi);
            modes.Add(new PeakMode(location, loCi, hiCi,
                covered ? grid[sl] : grid[i], covered ? grid[sr] : grid[i], mass, rhat));
        }
        modes.Sort((u, v) => v.Mass.CompareTo(u.Mass));   // dominant first
        return modes;
    }

    // Matching-free between-chain R̂ on a mode's location: each chain's λ over the support [lo,hi] is that chain's
    // own peak-location histogram there, giving its mean and within-chain variance directly. Gelman–Rubin on those
    // (B = Var of the chain means, W = mean within-variance, n = mean per-chain peak count) → ≈1 when chains agree.
    private static double ModeRHat(double[][] chainLambda, double[] grid, int lo, int hi)
    {
        int m = chainLambda.Length;
        var mean = new double[m];
        var varc = new double[m];
        double nbar = 0.0;
        int seen = 0;
        for (int c = 0; c < m; c++)
        {
            double cnt = 0.0, s = 0.0, s2 = 0.0;
            for (int j = lo; j <= hi; j++)
            {
                double w = chainLambda[c][j];
                cnt += w; s += w * grid[j]; s2 += w * grid[j] * grid[j];
            }
            if (cnt <= 0.0) continue;
            mean[seen] = s / cnt;
            varc[seen] = Math.Max(0.0, s2 / cnt - mean[seen] * mean[seen]);
            nbar += cnt;
            seen++;
        }
        if (seen < 2) return double.NaN;        // can't assess convergence from one chain
        nbar = Math.Max(1.0, nbar / seen);

        double gm = 0.0;
        for (int c = 0; c < seen; c++) gm += mean[c];
        gm /= seen;
        double between = 0.0, within = 0.0;
        for (int c = 0; c < seen; c++) { between += (mean[c] - gm) * (mean[c] - gm); within += varc[c]; }
        between /= seen - 1;                     // Var of the chain means (= B/n in Gelman–Rubin)
        within /= seen;
        if (within <= 0.0) return double.NaN;
        double varPlus = (nbar - 1.0) / nbar * within + between;
        return Math.Sqrt(varPlus / within);
    }

    // Step-CDF quantile of λ-mass over the support [lo, hi].
    private static double MassQuantile(double[] grid, double[] w, int lo, int hi, double total, double q)
    {
        double target = q * total, cum = 0.0;
        for (int j = lo; j <= hi; j++) { cum += w[j]; if (cum >= target) return grid[j]; }
        return grid[hi];
    }

    private static KnotConfig OverdispersedStart(Xoshiro256PlusPlus rng, int dispersion)
    {
        int k = dispersion <= 0 ? 0 : rng.NextInt(dispersion + 1);
        var knots = new double[k];
        for (int i = 0; i < k; i++) knots[i] = rng.NextDouble();
        Array.Sort(knots);
        return new KnotConfig(knots);
    }

    /// <summary>
    /// One BARS reversible-jump replica + its thread-local accumulators. As an <see cref="IChain{KnotConfig}"/> its
    /// <see cref="Step"/> advances the chain (or tempered ladder) and runs the resampler block, and its acceptance
    /// is banked-inclusive (it survives the resampler's chain rebuilds); as an <see cref="IChainRun{KnotConfig}"/>
    /// it burns adaptively, reduces a draw to functionals, and carries the BARS sinks (peak cloud, λ/π fields).
    /// </summary>
    private sealed class BarsChainRun : IChain<KnotConfig>, IChainRun<KnotConfig>
    {
        private readonly BarsEnsemble _owner;
        public readonly Xoshiro256PlusPlus Rng;
        public readonly double[] PeakIntensitySum;
        public readonly double[] SpanCoverageSum;
        public readonly List<double> PeakDraws = new();
        public KnotConfig Config;
        public double[]? W;
        public IKnotKernel Kernel = null!;
        public Xoshiro256PlusPlus ReadoutRng = null!;   // per-draw coefficient draws, decoupled from the chain RNG
        public ParallelTempering<KnotConfig>? Ladder;   // set ⇒ the RJ chain is this ladder's cold replica
        public readonly Dictionary<string, long> MoveAtt = new();   // per-move attempts, banked across rebuilds
        public readonly Dictionary<string, long> MoveAcc = new();   // per-move accepts
        public double PeakLocSum;
        public double PeakLocSumSq;
        public double PeakHeightSum;
        public double PeakCountSum;
        public long BankedAccepted;
        public long BankedAttempts;
        public int SinceWeight;

        private ReversibleJumpChain<KnotConfig> _rjChain = null!;

        public BarsChainRun(BarsEnsemble owner, Xoshiro256PlusPlus rng, KnotConfig start, int gridSize)
        {
            _owner = owner;
            Rng = rng;
            Config = start;
            PeakIntensitySum = new double[gridSize];
            SpanCoverageSum = new double[gridSize];
        }

        /// <summary>The live RJ chain — for the tempered case, the ladder's cold replica (reassigned on resampler rebuilds).</summary>
        public ReversibleJumpChain<KnotConfig> RjChain => _rjChain;
        public void SetChain(ReversibleJumpChain<KnotConfig> chain) => _rjChain = chain;

        // IChain — banked-inclusive so the resampler's chain rebuilds don't drop acceptance history.
        public KnotConfig Step() => _owner.StepOnce(this);
        public long Accepted => BankedAccepted + _rjChain.Accepted;
        public long Attempts => BankedAttempts + _rjChain.Attempts;

        // IChainRun
        public IChain<KnotConfig> Chain => this;
        public void Burn(int steps) => _owner.AdaptiveBurn(this, steps);
        public void Accumulate(in KnotConfig draw, Span<double> functionals) => _owner.AccumulateDraw(this, draw, functionals);
    }
}
