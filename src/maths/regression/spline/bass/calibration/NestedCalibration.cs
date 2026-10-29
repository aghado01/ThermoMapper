using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The nested reference sampler's settings (D09, C3): how many emulator draws the outer loop takes, and each inner run's
/// ensemble, which runs the conditional sampler (<see cref="CalibrationSampling"/>'s burn-in, thinning, decorrelation
/// and ladder) until N14's verdict reads <see cref="ConvergenceOutcome.Converged"/> or its budget runs out. Values only.
/// </summary>
public sealed record ReferenceSampling
{
    /// <summary>The emulator draws S the outer loop takes, without replacement from the retained draws (by weight, with
    /// replacement, when the draws are weighted); null takes every retained draw once.</summary>
    public int? OuterDraws { get; init; }

    /// <summary>Each inner run's independent chains (R̂ needs at least 2).</summary>
    public int Chains { get; init; } = 4;

    /// <summary>Each inner chain's budget of kept draws: the run stops earlier when its verdict reads converged.</summary>
    public int MaxSamples { get; init; } = 4000;

    /// <summary>Kept draws per chain between the verdict's readings.</summary>
    public int BatchSize { get; init; } = 250;

    /// <summary>The R̂ target the verdict reads; N14's reference, 1.01, by default.</summary>
    public double RHatTarget { get; init; } = ChainEnsemble.ReferenceRHatTarget;

    /// <summary>The bulk and tail ESS per chain every monitored functional needs; N14's default, 100.</summary>
    public double MinEssPerChain { get; init; } = ChainEnsemble.MinEssPerChain;

    /// <summary>Checks the fields, naming each in the reason.</summary>
    public void Validate()
    {
        if (OuterDraws is < 1) throw Bad("OuterDraws", "must be at least 1");
        if (Chains < 2) throw Bad("Chains", "must be at least 2, for R̂");
        if (MaxSamples < 1) throw Bad("MaxSamples", "must be at least 1");
        if (BatchSize < 1) throw Bad("BatchSize", "must be at least 1");
        if (!(RHatTarget > 1.0 && double.IsFinite(RHatTarget))) throw Bad("RHatTarget", "must be finite and above 1");
        if (!(MinEssPerChain >= 0.0 && double.IsFinite(MinEssPerChain))) throw Bad("MinEssPerChain", "must be finite and not negative");
    }

    private static ArgumentException Bad(string field, string reason) => new($"CalibrationConfig.Reference.{field}: {reason}");
}

/// <summary>
/// One inner run of the nested sampler: the emulator draw e_s it conditions on, its ensemble's verdict (N14), the kept
/// draws of each chain in order, and the run's lengths: burn-in and thinning per chain, kept draws per chain.
/// </summary>
public sealed record InnerRun(int EmulatorDraw, ConvergenceVerdict Verdict, IReadOnlyList<IReadOnlyList<CalibrationDraw>> Chains,
                              int BurnIn, int Thin, int SamplesPerChain, double AcceptanceRate)
{
    /// <summary>Whether the verdict reads converged, the condition for pooling.</summary>
    public bool Converged => Verdict.Outcome == ConvergenceOutcome.Converged;
}

/// <summary>
/// The nested reference sampler's result (D09's outputs): every inner run with its verdict, and the pooled posterior over
/// the converged ones. A run that did not converge is kept and reported (<see cref="NotConverged"/>), never pooled
/// unless asked.
/// </summary>
public sealed record NestedCalibration(IReadOnlyList<CalibrationParameter> Parameters, IReadOnlyList<InnerRun> Runs)
{
    /// <summary>The runs whose verdict did not read converged.</summary>
    public IReadOnlyList<InnerRun> NotConverged => [.. Runs.Where(r => !r.Converged)];

    /// <summary>
    /// The cut posterior's draws, p_cut(θ | y) = (1/S′)·Σ_s p(θ | y, e_s) over the S′ pooled runs: each run's draws, all its
    /// chains together, weighted 1/(S′·n_s), so every emulator draw counts once whatever its run's length. Pools the
    /// converged runs, or every run when <paramref name="includeNotConverged"/> says so.
    /// </summary>
    public (IReadOnlyList<CalibrationDraw> Draws, double[] Weights) Pooled(bool includeNotConverged = false)
    {
        InnerRun[] runs = [.. Runs.Where(r => includeNotConverged || r.Converged)];
        var draws = new List<CalibrationDraw>();
        var weights = new List<double>();
        foreach (InnerRun run in runs)
        {
            int n = run.Chains.Sum(c => c.Count);
            foreach (IReadOnlyList<CalibrationDraw> chain in run.Chains)
                foreach (CalibrationDraw draw in chain)
                {
                    draws.Add(draw);
                    weights.Add(1.0 / (runs.Length * (double)n));
                }
        }
        return (draws, [.. weights]);
    }
}

/// <summary>
/// The nested conditional reference sampler (D09, C3): for s = 1, …, S it takes an emulator draw e_s from the retained
/// draws and runs the conditional sampler on p(θ, σ², v | y, e_s) as a <see cref="ChainEnsemble"/> from prior draws,
/// in rounds, until N14's verdict reads converged or the budget runs out, keeping every inner run's draws and verdict.
/// The observations never update e: that is the cut. R refreshes e every iteration and takes a few θ steps under each,
/// which does not equilibrate under any e; this sampler does, and the error left by finite inner runs is measured by
/// varying their length.
/// </summary>
public static class NestedCalibrationSampler
{
    /// <summary>Runs the sampler on <paramref name="model"/> under its config's <see cref="CalibrationConfig.Sampling"/>
    /// and <see cref="CalibrationConfig.Reference"/>.</summary>
    public static NestedCalibration Run(CalibrationModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        CalibrationSampling sampling = model.Config.Sampling;
        ReferenceSampling reference = model.Config.Reference;
        sampling.Validate();
        reference.Validate();
        int[] seeds = SeedTree.Derive(sampling.MasterSeed, 2);
        int[] outer = OuterDraws(model.Emulator, reference.OuterDraws, new Xoshiro256PlusPlus(seeds[0]));
        int[] innerSeeds = SeedTree.Derive(seeds[1], outer.Length);
        var runs = new InnerRun[outer.Length];
        var ensemble = new ChainEnsemble();
        for (int s = 0; s < outer.Length; s++)
        {
            var inner = new ConditionalEnsemble(model.Target(outer[s]), sampling);
            (EnsembleRun summary, IReadOnlyList<IChainRun<CalibrationState>> chains) = ensemble.Run(
                inner, reference.Chains, innerSeeds[s], sampling.BurnIn, reference.MaxSamples, reference.RHatTarget,
                Math.Min(reference.BatchSize, reference.MaxSamples), reference.MinEssPerChain);
            runs[s] = new InnerRun(outer[s], summary.Verdict, [.. chains.Select(c => (IReadOnlyList<CalibrationDraw>)((ConditionalChainRun)c).Draws)],
                                   sampling.BurnIn, sampling.Thin, summary.SamplesUsed, summary.AcceptanceRate);
        }
        return new NestedCalibration(model.Config.Parameters, runs);
    }

    // The outer loop's emulator draws: every retained draw once; a sample without replacement of S of them; or, for
    // weighted draws, S draws by weight with replacement.
    internal static int[] OuterDraws(ICalibrationEmulator emulator, int? count, Xoshiro256PlusPlus rng)
    {
        int total = emulator.DrawCount;
        if (emulator.DrawWeights is { } weights)
        {
            int s = count ?? total;
            double sum = weights.Sum();
            if (!(sum > 0.0) || weights.Any(w => !(w >= 0.0) || double.IsPositiveInfinity(w)))
                throw new ArgumentException("The emulator's draw weights must be finite, not negative and not all zero.", nameof(emulator));
            var drawn = new int[s];
            for (int i = 0; i < s; i++)
            {
                double u = rng.NextDouble() * sum, cumulative = 0.0;
                drawn[i] = total - 1;
                for (int d = 0; d < total; d++)
                    if ((cumulative += weights[d]) > u) { drawn[i] = d; break; }
            }
            return drawn;
        }
        if (count is not { } k || k >= total) return [.. Enumerable.Range(0, total)];
        int[] all = [.. Enumerable.Range(0, total)];
        for (int i = 0; i < k; i++)   // a partial Fisher–Yates: the first k are a uniform sample without replacement
        {
            int j = i + rng.NextInt(total - i);
            (all[i], all[j]) = (all[j], all[i]);
        }
        return all[..k];
    }
}

/// <summary>
/// One emulator draw's conditional as an ensemble model (N14's shell): each chain a <see cref="CalibrationChain"/> from
/// a prior draw, overdispersed against the posterior, thinned by <see cref="CalibrationSampling.Thin"/> scans per kept
/// draw, adapting during burn-in only. Its functionals, all monitored, in four groups the verdict reads together: θ (each
/// continuous coordinate in [0, 1], and an indicator per categorical level but the last), the variances (ln σ²_g), the
/// discrepancy coefficients v, and the fit (the log-likelihood).
/// </summary>
internal sealed class ConditionalEnsemble : IEnsembleModel<CalibrationState>
{
    private readonly CalibrationTarget _target;
    private readonly CalibrationSampling _sampling;
    private readonly int _continuous, _groups, _coefficients, _levels;

    public ConditionalEnsemble(CalibrationTarget target, CalibrationSampling sampling)
    {
        _target = target;
        _sampling = sampling;
        CalibrationModel model = target.Model;
        _continuous = model.Layout.Continuous.Count;
        _groups = model.Measurement.Groups;
        _coefficients = model.Coefficients;
        _levels = model.Layout.Categorical.Sum(c => c.Levels.Length - 1);
        var names = new List<string>();
        foreach (ContinuousParameter p in model.Layout.Continuous) names.Add(p.Name);
        foreach (CategoricalParameter c in model.Layout.Categorical)
            for (int l = 0; l + 1 < c.Levels.Length; l++) names.Add($"{c.Name}={c.Levels[l]}");
        for (int g = 0; g < _groups; g++) names.Add($"ln σ²[{g}]");
        for (int k = 0; k < _coefficients; k++) names.Add($"v[{k}]");
        names.Add("log-likelihood");
        FunctionalNames = names;
        int theta = _continuous + _levels;
        var groups = new List<FunctionalGroup>();
        if (theta > 0) groups.Add(new FunctionalGroup("theta", [.. Enumerable.Range(0, theta)]));
        groups.Add(new FunctionalGroup("variances", [.. Enumerable.Range(theta, _groups)]));
        if (_coefficients > 0) groups.Add(new FunctionalGroup("discrepancy", [.. Enumerable.Range(theta + _groups, _coefficients)]));
        groups.Add(new FunctionalGroup("fit", [theta + _groups + _coefficients]));
        Groups = groups;
    }

    public int FunctionalDim => _continuous + _levels + _groups + _coefficients + 1;
    public int MonitoredDim => FunctionalDim;
    public IReadOnlyList<string> FunctionalNames { get; }
    public IReadOnlyList<FunctionalGroup> Groups { get; }

    public IChainRun<CalibrationState> StartChain(int chainSeed, int readoutSeed)
    {
        CalibrationState start = _target.Model.Start(new Xoshiro256PlusPlus(readoutSeed));
        return new ConditionalChainRun(new CalibrationChain(_target, _sampling, start, chainSeed), _sampling.Thin, this);
    }

    public void Reduce(CalibrationChain chain, Span<double> functionals)
    {
        CalibrationState state = chain.Current;
        int f = 0;
        for (int i = 0; i < _continuous; i++) functionals[f++] = state.Unit[i];
        ParameterLayout layout = _target.Model.Layout;
        for (int j = 0; j < layout.Categorical.Count; j++)
            for (int l = 0; l + 1 < layout.Categorical[j].Levels.Length; l++) functionals[f++] = state.Levels[j] == l ? 1.0 : 0.0;
        for (int g = 0; g < _groups; g++) functionals[f++] = Math.Log(state.Variances[g]);
        for (int k = 0; k < _coefficients; k++) functionals[f++] = state.Coefficients[k];
        functionals[f] = ((CalibrationEvaluation)chain.Cold.CurrentEvaluation).LogLikelihoodValue;
    }

    public CalibrationDraw Draw(CalibrationChain chain) => CalibrationDraw.Of(_target, chain.Current, chain.Cold.CurrentEvaluation);
}

/// <summary>One inner chain: the calibration chain thinned to one kept draw per <c>thin</c> scans, adapting during burn-in,
/// keeping its draws.</summary>
internal sealed class ConditionalChainRun : IChainRun<CalibrationState>, IChain<CalibrationState>
{
    private readonly CalibrationChain _chain;
    private readonly int _thin;
    private readonly ConditionalEnsemble _model;
    private readonly List<CalibrationDraw> _draws = [];

    public ConditionalChainRun(CalibrationChain chain, int thin, ConditionalEnsemble model)
    {
        _chain = chain;
        _thin = thin;
        _model = model;
    }

    public IChain<CalibrationState> Chain => this;
    public IReadOnlyList<CalibrationDraw> Draws => _draws;
    public long Accepted => _chain.Accepted;
    public long Attempts => _chain.Attempts;

    public CalibrationState Step()
    {
        for (int t = 0; t < _thin; t++) _chain.Step();
        return _chain.Current;
    }

    public void Burn(int steps) => _chain.Burn(steps);

    public void Accumulate(in CalibrationState draw, Span<double> functionals)
    {
        _model.Reduce(_chain, functionals);
        _draws.Add(_model.Draw(_chain));
    }
}
