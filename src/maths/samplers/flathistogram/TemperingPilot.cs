using System;
using System.Collections.Generic;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// A tempering pilot's rung draws for a flat-histogram warm start: each rung's inverse temperature b_i and its draws of
/// the untempered log likelihood ℓ₁, sampled from <c>prior·e^{b_i·ℓ₁}</c>, the linear likelihood path. What
/// <see cref="Wham"/> combines into ln g, and what a window, b_lo and K are read from. A rung at b = 1 is required, as
/// the window's top is read there. Each rung carries its statistical inefficiency g_i = 1 + 2τ_i, the draws per effective
/// draw, by which WHAM weighs it (Ferrenberg and Swendsen 1989).
/// </summary>
public sealed class TemperingPilot
{
    private readonly double[] _inverseTemperatures;
    private readonly double[][] _levels;
    private readonly double[] _inefficiencies;

    /// <param name="inverseTemperatures">The rungs' b, distinct, in [0, 1], one of them exactly 1.</param>
    /// <param name="levels">Each rung's ℓ₁ draws, finite, at least two per rung; copied.</param>
    /// <param name="inefficiencies">Each rung's statistical inefficiency g_i ≥ 1, finite; null takes 1 for every rung,
    /// the value for independent draws.</param>
    public TemperingPilot(IReadOnlyList<double> inverseTemperatures, IReadOnlyList<double[]> levels,
                          IReadOnlyList<double>? inefficiencies = null)
    {
        ArgumentNullException.ThrowIfNull(inverseTemperatures);
        ArgumentNullException.ThrowIfNull(levels);
        if (inverseTemperatures.Count < 1 || levels.Count != inverseTemperatures.Count)
            throw new ArgumentException("One array of draws per rung is required.", nameof(levels));
        _inverseTemperatures = [.. inverseTemperatures];
        _levels = new double[levels.Count][];
        var seen = new HashSet<double>();
        for (int i = 0; i < _inverseTemperatures.Length; i++)
        {
            double b = _inverseTemperatures[i];
            if (!(b >= 0.0 && b <= 1.0)) throw new ArgumentException($"Rung {i}'s b = {b} is outside [0, 1].", nameof(inverseTemperatures));
            if (!seen.Add(b)) throw new ArgumentException($"Two rungs share b = {b}.", nameof(inverseTemperatures));
            if (levels[i] is not { Length: >= 2 }) throw new ArgumentException($"Rung {i} needs at least two draws.", nameof(levels));
            foreach (double l in levels[i])
                if (!double.IsFinite(l)) throw new ArgumentException($"Rung {i} holds a draw ℓ₁ = {l}.", nameof(levels));
            _levels[i] = (double[])levels[i].Clone();
            Draws += levels[i].Length;
        }
        _inefficiencies = new double[_inverseTemperatures.Length];
        if (inefficiencies is not null && inefficiencies.Count != _inverseTemperatures.Length)
            throw new ArgumentException("One inefficiency per rung is required.", nameof(inefficiencies));
        for (int i = 0; i < _inefficiencies.Length; i++)
        {
            double g = inefficiencies?[i] ?? 1.0;
            if (!(g >= 1.0 && double.IsFinite(g)))
                throw new ArgumentException($"Rung {i}'s inefficiency {g} is not a finite number of at least 1.", nameof(inefficiencies));
            _inefficiencies[i] = g;
        }
        ColdRung = Array.IndexOf(_inverseTemperatures, 1.0);
        if (ColdRung < 0) throw new ArgumentException("The pilot needs a rung at b = 1.", nameof(inverseTemperatures));
        HotRung = 0;
        for (int i = 1; i < _inverseTemperatures.Length; i++) if (_inverseTemperatures[i] < _inverseTemperatures[HotRung]) HotRung = i;
    }

    /// <summary>The rungs' b.</summary>
    public IReadOnlyList<double> InverseTemperatures => _inverseTemperatures;

    /// <summary>The rungs.</summary>
    public int Rungs => _inverseTemperatures.Length;

    /// <summary>The draws over every rung.</summary>
    public long Draws { get; }

    /// <summary>The rung at b = 1.</summary>
    public int ColdRung { get; }

    /// <summary>The rung of least b.</summary>
    public int HotRung { get; }

    /// <summary>The least b, the pilot's hot end.</summary>
    public double HotEnd => _inverseTemperatures[HotRung];

    /// <summary>Rung <paramref name="rung"/>'s ℓ₁ draws.</summary>
    public ReadOnlySpan<double> Levels(int rung) => _levels[rung];

    /// <summary>The rungs' statistical inefficiencies g_i.</summary>
    public IReadOnlyList<double> Inefficiencies => _inefficiencies;

    /// <summary>
    /// The pilot from tempering ladders' rung traces (<see cref="ParallelTempering{TState}.TraceRungs"/>), each ladder on
    /// the same rungs, on the likelihood path, tracing slopes: a rung's slope ∂ℓ_b/∂b is ℓ₁ where ℓ_b is linear in b,
    /// which the draws are checked for. Every hotter rung's log ratio to the next colder one must be
    /// <c>(b_c − b)·slope</c>, to rounding; a collapsed likelihood, not linear in b, is refused, since WHAM needs the
    /// rungs' laws to be <c>prior·e^{b·ℓ₁}</c>. The ladders are traced for the same number of scans, so each rung's draws
    /// are one chain per ladder, and its inefficiency is <c>g = max(1, n/ESS)</c>, ESS the effective sample size of the
    /// mean of its ℓ₁ chains (<see cref="ChainDiagnostics.MeanEss"/>); 1 where every draw is the same.
    /// </summary>
    public static TemperingPilot FromLadders<TState>(IReadOnlyList<ParallelTempering<TState>> ladders)
    {
        ArgumentNullException.ThrowIfNull(ladders);
        if (ladders.Count == 0) throw new ArgumentException("At least one ladder is required.", nameof(ladders));
        double[] b = ladders[0].InverseTemperatures;
        var chains = new double[b.Length][][];
        for (int i = 0; i < b.Length; i++) chains[i] = new double[ladders.Count][];
        int scans = -1;
        int l = 0;
        foreach (ParallelTempering<TState> ladder in ladders)
        {
            ArgumentNullException.ThrowIfNull(ladder);
            if (ladder.Path != TemperingPath.Likelihood) throw new ArgumentException("A pilot runs on the likelihood path.", nameof(ladders));
            if (ladder.RungTraces is not IReadOnlyList<ChainTrace> traces || !ladder.TracesSlopes)
                throw new ArgumentException("Every ladder must trace its rungs with slopes.", nameof(ladders));
            double[] rungs = ladder.InverseTemperatures;
            if (rungs.Length != b.Length) throw new ArgumentException("The ladders hold different rungs.", nameof(ladders));
            for (int i = 0; i < b.Length; i++)
                if (rungs[i] != b[i]) throw new ArgumentException("The ladders hold different rungs.", nameof(ladders));
            if (scans >= 0 && traces[0].Count != scans) throw new ArgumentException("The ladders were traced for different numbers of scans.", nameof(ladders));
            scans = traces[0].Count;
            for (int i = 0; i < b.Length; i++)
            {
                double[] slopes = traces[i].Draws(RungFunctional.Slope);
                if (i > 0)
                {
                    double[] ratios = traces[i].Draws(RungFunctional.LogRatio);
                    double step = b[i - 1] - b[i];
                    for (int s = 0; s < slopes.Length; s++)
                    {
                        double residual = ratios[s] - step * slopes[s];
                        if (!(Math.Abs(residual) <= 1e-9 * (1.0 + Math.Abs(ratios[s]))))
                            throw new ArgumentException(
                                $"Rung {i}'s log ratio departs from (b_c − b)·ℓ₁ by {residual:G3}: the target's ℓ_b is not linear in b.", nameof(ladders));
                    }
                }
                chains[i][l] = slopes;
            }
            l++;
        }
        var levels = new double[b.Length][];
        var inefficiencies = new double[b.Length];
        for (int i = 0; i < b.Length; i++)
        {
            var pooled = new List<double>(ladders.Count * scans);
            foreach (double[] chain in chains[i]) pooled.AddRange(chain);
            levels[i] = [.. pooled];
            EssDiagnostic ess = ChainDiagnostics.MeanEss(chains[i]);
            inefficiencies[i] = ess.State == DiagnosticState.Estimated ? Math.Max(1.0, pooled.Count / ess.Value) : 1.0;
        }
        return new TemperingPilot(b, levels, inefficiencies);
    }
}
