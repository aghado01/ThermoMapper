using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's fit by the flat-histogram sampler (N17 Phase 2): <see cref="MulticanonicalRun"/> over D04's sweep on π_w
/// (<see cref="BassMulticanonicalModel"/>), then assembled into a <see cref="BassFit"/> whose draws are the production
/// chains', weighted at b = 1.
/// <list type="bullet">
/// <item><b>The pilot</b> (<see cref="Pilot"/>), when the run reads one: <see cref="BassFlatHistogramPilot.Chains"/> ladders
/// on the linear form, each from an exact draw of the prior on F, burned as a tempered fit burns them (adapted after rounds
/// of 1, 2, 4, … scans, a self-sizing ladder sized after every round but the last), then frozen and traced. Each ladder is
/// read by <see cref="TemperingPilot.FromLadders"/>, since self-sizing ladders freeze their own rungs, and the rungs pool
/// by b: a rung several ladders share (b = 1, b = 0) holds their draws, its inefficiency their draws over their summed
/// effective draws.</item>
/// <item><b>Streams.</b> The master seed derives the run's root and the pilot's (<see cref="SeedTree"/>); the run's root
/// derives each chain's two seeds, as the ensemble shell does, and the pilot's root each ladder's start and stream.</item>
/// <item><b>Draws.</b> Each chain's production draws, composed as an MCMC fit's retained draws are, so each keeps its full
/// state (D04), in one store per chain; the terms merge in chain order. <see cref="BassFit.DrawWeights"/> are
/// e^{ℓ₁}/w(ℓ₁) self-normalized over every chain's draws (<see cref="Reweighting.Weigh"/>), chain-major in draw order,
/// the stores' pooled order, so every view reads them as it reads an SMC fit's.</item>
/// <item><b>The evidence.</b> The ratio estimator's ln Z(b) on a grid from b_lo to 1 (<see cref="BassFlatHistogramEvidence"/>),
/// under proper σ² and intercept priors; the fit's evidence ln Z(1) only from b_lo = 0.</item>
/// <item><b>The budget.</b> A call takes at most <see cref="BassFlatHistogramSpec.MaxIterations"/> iterations; a fit short
/// of the lengths is flagged, reports no evidence or residual check, and its checkpoint continues it. A complete run's
/// checkpoint is null: its lengths are totals.</item>
/// <item><b>Restart state.</b> The engine's checkpoint, which carries the learning, both phases' diagnostics and every
/// chain's own state (<see cref="BassMulticanonicalChainState"/>), fills the ensemble checkpoint's one engine slot, since
/// the chains share the learned weights.</item>
/// </list>
/// </summary>
internal static class BassFlatHistogramFit
{
    /// <summary>The fit of <paramref name="model"/>, whose config sets <see cref="BassSampling.FlatHistogram"/>, or the
    /// continuation of <paramref name="resumeFrom"/>, which must be a flat-histogram fit's checkpoint on the same config
    /// and data.</summary>
    public static BassFit Run(BassModel model, string fingerprint, BassRunOptions options, EnsembleCheckpoint? resumeFrom)
    {
        BassSampling sampling = model.Config.Sampling;
        BassFlatHistogramSpec spec = sampling.FlatHistogram ?? throw new ArgumentException("The config does not run the flat histogram.", nameof(model));
        FlatHistogramSpec run = spec.Run;
        BassFlatHistogramState? restart = null;
        if (resumeFrom is not null)
        {
            if (resumeFrom.Chains.Count != 1 || resumeFrom.Chains[0].Engine is not BassFlatHistogramState state)
                throw new ArgumentException("Not a BASS flat-histogram fit's checkpoint.", "checkpoint");
            restart = state;
        }
        int[] roots = SeedTree.Derive(sampling.MasterSeed, 2);   // the run's, the pilot's
        bool piloted = restart is null && run.NeedsPilot;
        long pilotScans = piloted ? (long)spec.Pilot!.Chains * (spec.Pilot.BurnIn + spec.Pilot.TracedScans) : 0;
        var sweeps = new BassSweepCount(options, pilotScans + (long)run.Chains * (run.LearningIterations + run.ProductionDraws) * sampling.Thin);
        TemperingPilot? pilot = null;
        BassFlatHistogramPilotReadout? pilotReadout = restart?.Pilot;
        MulticanonicalWeights? shape = null;
        if (piloted)
        {
            (pilot, BassTemperingReadout ladders) = Pilot(model, spec.Pilot!, roots[1], options, sweeps);
            shape = StartShape(run, pilot);
            pilotReadout = new BassFlatHistogramPilotReadout(ladders, [.. pilot.InverseTemperatures], [.. pilot.Inefficiencies], pilot.Draws,
                                                             NodeDraws(pilot, shape));
        }
        var client = new BassMulticanonicalModel(model, fingerprint, run.LearningIterations, sweeps);
        var parallel = new ParallelOptions { MaxDegreeOfParallelism = options.MaxDegreeOfParallelism, CancellationToken = options.Cancellation };
        MulticanonicalResult<BassState> result;
        try
        {
            result = MulticanonicalRun.Run(run, model.Target, client, roots[0], pilot, restart?.Run, spec.MaxIterations, options.Checkpoint, parallel);
        }
        catch (InvalidOperationException e) when (shape is not null && run.Start == FlatHistogramStart.Wham)
        {
            throw new InvalidOperationException(
                $"The flat-histogram run failed from its WHAM start over the window [{shape.Low:G6}, {shape.High:G6}], K = {shape.Nodes}, b_lo = " +
                $"{shape.LowSlope}, the pilot having drawn nothing under nodes [{string.Join(", ", pilotReadout!.EmptyNodes)}]: {e.Message}", e);
        }
        sweeps.End();
        MulticanonicalWeights start = result.WeightTrace[0].Weights;
        if (shape is not null && (start.Low, start.High, start.Nodes, start.LowSlope) != (shape.Low, shape.High, shape.Nodes, shape.LowSlope))
            throw new InvalidOperationException("The run resolved another window, K or b_lo than the pilot's readout.");
        return Assemble(model, fingerprint, spec, roots, result, client, pilotReadout, options.Checkpoint);
    }

    /// <summary>
    /// The tempering pilot of <paramref name="model"/> under <paramref name="spec"/> from <paramref name="seed"/>: each ladder on
    /// the linear form from an exact draw of the prior on F, burned, frozen and traced, read by
    /// <see cref="TemperingPilot.FromLadders"/> and pooled by b (<see cref="Pooled"/>); with the ladders' readout over their
    /// traced scans. The ladders run in parallel within the run's bound; nothing depends on it.
    /// </summary>
    internal static (TemperingPilot Pilot, BassTemperingReadout Ladders) Pilot(BassModel model, BassFlatHistogramPilot spec, int seed,
                                                                            BassRunOptions options, BassSweepCount sweeps)
    {
        BassModel linear = BassModel.Create(model.Config, model.Data, TemperedForm.Linear);
        int[] seeds = SeedTree.Derive(seed, spec.Chains);
        var ladders = new BassLadder[spec.Chains];
        Parallel.For(0, spec.Chains, new ParallelOptions { MaxDegreeOfParallelism = options.MaxDegreeOfParallelism, CancellationToken = options.Cancellation }, c =>
        {
            int[] streams = SeedTree.Derive(seeds[c], 2);   // the start's draws, the ladder's
            BassState start = new BassPriorSampler(linear).DrawState(new Xoshiro256PlusPlus(streams[0])).State;
            var ladder = new BassLadder(linear, spec.Tempering, start.Terms, start.LogTau, streams[1]);
            void Scan()
            {
                options.Cancellation.ThrowIfCancellationRequested();
                ladder.Step();
                sweeps.Add(1);
            }
            int done = 0;
            if (ladder.Spec.Adaptive)
                for (int round = 1; done + round <= spec.BurnIn; round *= 2)
                {
                    for (int k = 0; k < round; k++) Scan();
                    done += round;
                    if (ladder.SizesItself && done + 2L * round <= spec.BurnIn) ladder.Resize();
                    else ladder.Adapt();
                }
            for (; done < spec.BurnIn; done++) Scan();
            ladder.EndBurnIn(traceRungs: true);
            for (int s = 0; s < spec.TracedScans; s++) Scan();
            ladders[c] = ladder;
        });
        TemperingPilot pooled = Pooled([.. ladders.Select(l => TemperingPilot.FromLadders([l.Tempering]))]);
        return (pooled, BassTemperingReadout.Of(ladders));
    }

    /// <summary>
    /// <paramref name="pilots"/>' rungs pooled by b, in order of first appearance: a rung several pilots share holds their
    /// draws in pilot order, and its inefficiency is their draws over their summed effective draws, Σn / Σ(n/g).
    /// </summary>
    internal static TemperingPilot Pooled(IReadOnlyList<TemperingPilot> pilots)
    {
        var b = new List<double>();
        var levels = new List<List<double>>();
        var effective = new List<double>();
        var index = new Dictionary<double, int>();
        foreach (TemperingPilot p in pilots)
            for (int i = 0; i < p.Rungs; i++)
            {
                double rung = p.InverseTemperatures[i];
                if (!index.TryGetValue(rung, out int k))
                {
                    index.Add(rung, k = b.Count);
                    b.Add(rung);
                    levels.Add([]);
                    effective.Add(0.0);
                }
                ReadOnlySpan<double> draws = p.Levels(i);
                foreach (double l in draws) levels[k].Add(l);
                effective[k] += draws.Length / p.Inefficiencies[i];
            }
        var inefficiencies = new double[b.Count];
        for (int k = 0; k < inefficiencies.Length; k++) inefficiencies[k] = Math.Max(1.0, levels[k].Count / effective[k]);
        return new TemperingPilot(b, [.. levels.Select(l => l.ToArray())], inefficiencies);
    }

    // The start's window, b_lo and K as the run resolves them: the spec's, or the pilot's where it leaves them.
    private static MulticanonicalWeights StartShape(FlatHistogramSpec run, TemperingPilot pilot)
    {
        FlatHistogramWindow window = run.Window ?? Wham.Window(pilot, run.WindowQuantile, run.WindowMargin);
        int nodes = run.Nodes ?? Wham.Nodes(pilot, window);
        return MulticanonicalWeights.Flat(window.Low, window.High, run.LowSlope ?? pilot.HotEnd, nodes);
    }

    // The pilot's draws under each node's hat over `shape`'s window, every rung's: where a node holds none, WHAM's ĝ has no
    // mass and the warm start interpolates its weight.
    private static double[] NodeDraws(TemperingPilot pilot, MulticanonicalWeights shape)
    {
        var counts = new double[shape.Nodes];
        for (int i = 0; i < pilot.Rungs; i++)
            foreach (double level in pilot.Levels(i))
            {
                if (!shape.TryLocate(level, out int segment, out double fraction)) continue;
                counts[segment] += 1.0 - fraction;
                counts[segment + 1] += fraction;
            }
        return counts;
    }

    private static BassFit Assemble(BassModel model, string fingerprint, BassFlatHistogramSpec spec, int[] roots, MulticanonicalResult<BassState> result,
                                    BassMulticanonicalModel client, BassFlatHistogramPilotReadout? pilot, bool checkpoint)
    {
        BassConfig config = model.Config;
        BassSampling sampling = config.Sampling;
        int chains = spec.Run.Chains;
        var runs = new BassMulticanonicalChain[chains];
        for (int c = 0; c < chains; c++) runs[c] = (BassMulticanonicalChain)result.Chains[c];

        // The draws, each chain's table merged in chain order, and their weights at b = 1.
        var terms = new BassTermTable();
        var stores = new DrawStore<BassDraw>[chains];
        for (int c = 0; c < chains; c++)
        {
            BassTermTable own = runs[c].Terms;
            var map = new int[own.Count];
            for (int id = 0; id < map.Length; id++) map[id] = terms.Intern(own[id]);
            stores[c] = new DrawStore<BassDraw>();
            foreach (BassDraw d in runs[c].Draws) stores[c].Add(BassEnsemble.Remap(d, map));
        }
        MulticanonicalDraws? draws = result.Draws;
        int perChain = draws?.Draws ?? 0;
        for (int c = 0; c < chains; c++)
            if (stores[c].Count != perChain) throw new InvalidOperationException($"Chain {c} kept {stores[c].Count} draws, production {perChain}.");
        WeightedDraws? posterior = draws is null ? null : Reweighting.Weigh(draws, 1.0);
        double[] weights = posterior?.Weights ?? [];

        // The accounts over production (all the run while it learns), each chain's weighted share at M_max.
        var accounts = new BassChainAccount[chains];
        var refusals = new Tally<BassRefusal>();
        var support = new BassSupportTally(model.Prior);
        var kinds = new BassMoveTally();
        long emptyBirths = 0;
        double massAtCap = 0.0;
        MoveStat[]? moves = null;
        for (int c = 0; c < chains; c++)
        {
            BassMulticanonicalChain run = runs[c];
            BassStructureMoves structure = run.Chain.Structure;
            IReadOnlyList<MoveStat> now = run.Chain.Moves.MoveStats();
            IReadOnlyList<MoveStat>? from = run.LearningMoves;
            var counted = new MoveStat[now.Count];
            for (int i = 0; i < now.Count; i++)
                counted[i] = from is null ? now[i] : new MoveStat(now[i].Key, now[i].Attempts - from[i].Attempts, now[i].Accepted - from[i].Accepted);
            moves ??= new MoveStat[now.Count];
            for (int i = 0; i < now.Count; i++)
                moves[i] = new MoveStat(counted[i].Key, moves[i].Attempts + counted[i].Attempts, moves[i].Accepted + counted[i].Accepted);

            double atCap = 0.0, mass = 0.0;
            for (int s = 0; s < perChain; s++)
            {
                double w = weights[c * perChain + s];
                mass += w;
                if (stores[c][s].Structure.TermIds.Length == model.MaxBasis) atCap += w;
            }
            massAtCap += atCap;
            accounts[c] = new BassChainAccount(run.Start, counted, structure.Kinds.Copy(), structure.Refusals.Copy(), structure.EmptyBirths,
                                               structure.Support.Copy(), mass > 0.0 ? atCap / mass : double.NaN);
            kinds.MergeFrom(structure.Kinds);
            refusals.MergeFrom(structure.Refusals);
            emptyBirths += structure.EmptyBirths;
            support.MergeFrom(structure.Support);
        }

        // π_w's convergence diagnostics and D10's verdict over production; unavailable while learning.
        EnsembleRun? production = result.Production;
        ConvergenceDiagnostic none = new(RHatDiagnostic.Unavailable, RHatDiagnostic.Unavailable, RHatDiagnostic.Unavailable,
                                         EssDiagnostic.Unavailable, EssDiagnostic.Unavailable);
        int f = 0;
        ConvergenceDiagnostic Next() => production?.Monitored[f++] ?? none;
        ConvergenceDiagnostic count = Next();
        ConvergenceDiagnostic? tau = model.Target.SamplesTau ? Next() : null;
        ConvergenceDiagnostic logPosterior = Next();
        ConvergenceDiagnostic sigma2 = Next();
        IReadOnlyList<BassMonitoredPoint> points = client.PointInputs.Points;
        var fitted = new BassMonitoredValue[points.Count];
        for (int p = 0; p < fitted.Length; p++) fitted[p] = new BassMonitoredValue(points[p], Next());
        ConvergenceVerdict verdict = production?.Verdict
            ?? new ConvergenceVerdict(ConvergenceOutcome.Undetermined, ChainEnsemble.ReferenceRHatTarget, sampling.MinEssPerChain,
                [.. new[] { BassDiagnostics.PredictiveGroup, BassDiagnostics.LevelGroup, BassDiagnostics.StructuralGroup }
                    .Select(g => new GroupVerdict(g, ConvergenceOutcome.Undetermined, false, []))]);
        double acceptance = production?.AcceptanceRate ?? (double)moves!.Sum(m => m.Accepted) / moves!.Sum(m => m.Attempts);

        // The flat-histogram readout, with the evidence curve of a complete run under proper σ² and intercept priors.
        bool complete = result.Complete;
        bool defined = sampling.ReportEvidence && model.Target.Sigma2 is not null && model.Target.Intercept != GPriorCentre.Flat;
        BassFlatHistogramEvidence? curve = complete && defined && draws is not null ? BassFlatHistogramEvidence.Of(draws) : null;
        var readout = new BassFlatHistogramReadout(result.Weights, result.Learning.GainScale, result.Learning.GainStart, result.Learning.Iterations,
                                                   result.LearningDiagnostics, result.WeightTrace, result.ProductionDiagnostics,
                                                   draws is null ? [] : PerChain(draws, result.ProductionDiagnostics!),
                                                   posterior?.EffectiveSampleSize ?? double.NaN, posterior?.ChainEffectiveSampleSizes ?? [], draws, curve,
                                                   pilot, complete);

        var diagnostics = new BassDiagnostics(count, tau, logPosterior, sigma2, fitted, moves!, acceptance, kinds, refusals, emptyBirths, support,
                                              massAtCap, perChain, false, verdict, Separation(draws, client, points.Count), accounts, null, null, readout);
        BassEvidence? evidence = curve is { IsEvidence: true } && model.Target.Scale.IsProper ? BassEvidence.OfFlatHistogram(model, curve) : null;
        EnsembleCheckpoint? restart = checkpoint && !complete
            ? new EnsembleCheckpoint(new RunIdentity(sampling.MasterSeed, chains, fingerprint), (int[])roots.Clone(), perChain,
                                     [new ChainCheckpoint(new BassFlatHistogramState(result.Checkpoint!, pilot), [], new ChainTrace(0), 0, 0)])
            : null;
        return new BassFit(config, model.Resolved, model.Transform, new RunIdentity(sampling.MasterSeed, chains, fingerprint), terms, stores, weights,
                           diagnostics, evidence, complete ? BassResidualDependence.Of(model, terms, stores, weights: weights) : null, restart);
    }

    // Each chain's production visits: its occupancy over the window from its own levels, its round trips the run's.
    private static FlatHistogramDiagnostics[] PerChain(MulticanonicalDraws draws, FlatHistogramDiagnostics merged)
    {
        MulticanonicalWeights w = draws.Weights;
        var chains = new FlatHistogramDiagnostics[draws.Chains];
        for (int c = 0; c < chains.Length; c++)
        {
            var occupancy = new double[w.Nodes];
            long inside = 0, below = 0, above = 0;
            foreach (double level in draws.Levels(c))
            {
                if (w.TryLocate(level, out int segment, out double fraction))
                {
                    occupancy[segment] += 1.0 - fraction;
                    occupancy[segment + 1] += fraction;
                    inside++;
                }
                else if (level < w.Low) below++;
                else if (level > w.High) above++;
            }
            chains[c] = new FlatHistogramDiagnostics(occupancy, w.FlatOccupancy(), inside, below, above, [merged.RoundTrips[c]], [merged.LastEnd[c]]);
        }
        return chains;
    }

    // The mode-separation readings from each chain's reweighted means at b = 1: the log posterior and the fitted values.
    private static BassModeSeparation Separation(MulticanonicalDraws? draws, BassMulticanonicalModel client, int points)
    {
        if (draws is null) return new BassModeSeparation(double.NaN, double.NaN);
        int level = client.LevelIndex, first = level + 2;
        double[] logPosterior = Reweighting.Expectation(draws, 1.0, level).ChainMeans;
        var fitted = new double[draws.Chains][];
        for (int c = 0; c < fitted.Length; c++) fitted[c] = new double[points];
        for (int p = 0; p < points; p++)
        {
            double[] means = Reweighting.Expectation(draws, 1.0, first + p).ChainMeans;
            for (int c = 0; c < fitted.Length; c++) fitted[c][p] = means[c];
        }
        return BassModeSeparation.Of(logPosterior, fitted);
    }
}

/// <summary>A flat-histogram fit's restart state: the engine's checkpoint, the learning's and both phases' state with
/// every chain's, and the pilot's readout (null without a pilot), which a resumed run reports without rerunning it.</summary>
internal sealed record BassFlatHistogramState(MulticanonicalCheckpoint Run, BassFlatHistogramPilotReadout? Pilot) : IChainState;
