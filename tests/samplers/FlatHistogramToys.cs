using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;

namespace Maths.Samplers.Tests;

/// <summary>
/// The flat-histogram engine's toy targets with exact references (N17 Phase 1). Toy 1 is a double well: x uniform on
/// [−2, 2], ℓ₁(x) = −c·(x² − 1)², whose density of states the prior's midpoint quadrature gives to any accuracy, and so
/// SAMC's limit θ* (<see cref="SamcLearning.FixedPoint"/>), the tempered laws and ln Z(b). Toy 2 is a q-state Potts model
/// on a periodic lattice, prior uniform over the qᴺ configurations and ℓ₁ = J·(agreeing neighbour pairs), whose
/// density of states enumerates exactly.
/// </summary>
internal static class FlatHistogramToys
{
    // ---- Toy 1: the double well ------------------------------------------------------------------------------------

    public const double WellHalfWidth = 2.0;

    /// <summary>The double well's target at depth <paramref name="c"/>.</summary>
    public static ISamplingTarget<double> DoubleWell(double c) => SamplingTarget.Separable<double>(
        x => Math.Abs(x) <= WellHalfWidth ? -Math.Log(2.0 * WellHalfWidth) : double.NegativeInfinity,
        x => WellLevel(x, c));

    public static double WellLevel(double x, double c)
    {
        double d = x * x - 1.0;
        return -c * d * d;
    }

    /// <summary>The prior's midpoint quadrature on <paramref name="cells"/> cells: each cell's centre and ℓ₁, its log
    /// prior mass ln(1/cells).</summary>
    public static (double[] X, double[] Levels, double[] LogMasses) WellQuadrature(double c, int cells)
    {
        var x = new double[cells];
        var levels = new double[cells];
        var masses = new double[cells];
        double h = 2.0 * WellHalfWidth / cells, logMass = -Math.Log(cells);
        for (int j = 0; j < cells; j++)
        {
            x[j] = -WellHalfWidth + (j + 0.5) * h;
            levels[j] = WellLevel(x[j], c);
            masses[j] = logMass;
        }
        return (x, levels, masses);
    }

    /// <summary>
    /// ln g(ℓ), the prior's density of states, in closed form: with r = √(−ℓ/c), the levels x² = 1 + r and, for r &lt; 1,
    /// x² = 1 − r, each at ±x, contribute p(x)/|ℓ′(x)| = (1/4)/(4c·|x|·r). −∞ outside (−9c, 0).
    /// </summary>
    public static double WellLogDensityOfStates(double level, double c)
    {
        if (!(level < 0.0 && level > -9.0 * c)) return double.NegativeInfinity;
        double r = Math.Sqrt(-level / c);
        double sum = 1.0 / Math.Sqrt(1.0 + r) + (r < 1.0 ? 1.0 / Math.Sqrt(1.0 - r) : 0.0);
        return Math.Log(sum / (8.0 * c * r));
    }

    /// <summary>The hats' occupancy inside the window, E_w[φ_k | in], under <c>prior·w</c>, by the quadrature: what a chain
    /// on these weights visits.</summary>
    public static double[] WellOccupancy(MulticanonicalWeights weights, double c, int cells = 400_000)
    {
        (_, double[] levels, _) = WellQuadrature(c, cells);
        var occupancy = new double[weights.Nodes];
        double max = double.NegativeInfinity;
        foreach (double l in levels) if (weights.TryLocate(l, out _, out _)) max = Math.Max(max, weights.LogWeight(l));
        double total = 0.0;
        foreach (double l in levels)
        {
            if (!weights.TryLocate(l, out int s, out double u)) continue;
            double p = Math.Exp(weights.LogWeight(l) - max);
            occupancy[s] += p * (1.0 - u);
            occupancy[s + 1] += p * u;
            total += p;
        }
        for (int k = 0; k < occupancy.Length; k++) occupancy[k] /= total;
        return occupancy;
    }

    /// <summary>ln Z(b) = ln E_prior[e^{b·ℓ₁}] by the quadrature.</summary>
    public static double WellLogEvidence(double c, double b, int cells = 400_000)
    {
        (_, double[] levels, double[] masses) = WellQuadrature(c, cells);
        var terms = new double[cells];
        for (int j = 0; j < cells; j++) terms[j] = masses[j] + b * levels[j];
        return LogSumExp(terms);
    }

    /// <summary>E_b[f(x)] under the tempered posterior prior·e^{b·ℓ₁}, by the quadrature.</summary>
    public static double WellExpectation(double c, double b, Func<double, double> f, int cells = 400_000)
    {
        (double[] x, double[] levels, _) = WellQuadrature(c, cells);
        double max = double.NegativeInfinity;
        foreach (double l in levels) max = Math.Max(max, b * l);
        double num = 0.0, den = 0.0;
        for (int j = 0; j < cells; j++)
        {
            double w = Math.Exp(b * levels[j] - max);
            num += w * f(x[j]);
            den += w;
        }
        return num / den;
    }

    /// <summary>The double well's model: one random-walk move per step; functionals ℓ₁, x and the indicator x &gt; 0.</summary>
    public sealed class WellModel(double c, double step, double start = 1.0, int? stuckSeed = null) : IMulticanonicalModel<double>
    {
        public double Depth => c;
        public int FunctionalDim => 3;
        public int MonitoredDim => 2;
        public IReadOnlyList<string> FunctionalNames => ["l1", "x", "right"];
        public string Fingerprint => new RunFingerprint().Add(c).Add(step).Add(start).Add(stuckSeed ?? 0).ToString();

        public IMulticanonicalChain<double> StartChain(MulticanonicalTarget<double> target, int chainSeed, int readoutSeed)
        {
            var chain = new ReversibleJumpChain<double>(Moves(chainSeed), target, start, new Xoshiro256PlusPlus(chainSeed));
            return new MulticanonicalChain<double>(chain, target, Reduce);
        }

        // The stuck chain is never resumed, so a restored chain takes the common step.
        public IMulticanonicalChain<double> RestoreChain(MulticanonicalTarget<double> target, IChainState state)
            => MulticanonicalChain<double>.Restore([new GaussianStep(step)], target, state, Reduce);

        // A chain whose seed is stuckSeed takes steps too small to cross the window.
        private IRjMove<double>[] Moves(int chainSeed) => [new GaussianStep(stuckSeed == chainSeed ? 1e-7 : step)];

        private void Reduce(in double x, Span<double> functionals)
        {
            functionals[0] = WellLevel(x, c);
            functionals[1] = x;
            functionals[2] = x > 0.0 ? 1.0 : 0.0;
        }
    }

    /// <summary>A Gaussian random walk on a scalar: a symmetric proposal.</summary>
    public sealed class GaussianStep(double scale) : IRjMove<double>
    {
        public string Key => "walk";
        public string ReverseKey => "walk";
        public double Weight(double state) => 1.0;
        public Proposal<double>? Propose(double current, Xoshiro256PlusPlus rng) => new Proposal<double>(current + scale * Normal.SampleStandard(rng), 0.0, 0.0);
    }

    /// <summary>ln w's departure from a reference, up to a constant: the RMS over the nodes <paramref name="use"/> marks
    /// of (a − b) less its mean there.</summary>
    public static double CentredRms(IReadOnlyList<double> a, IReadOnlyList<double> b, bool[]? use = null)
    {
        double mean = 0.0;
        int n = 0;
        for (int k = 0; k < a.Count; k++)
        {
            if (use is not null && !use[k]) continue;
            mean += a[k] - b[k];
            n++;
        }
        mean /= n;
        double ss = 0.0;
        for (int k = 0; k < a.Count; k++)
        {
            if (use is not null && !use[k]) continue;
            double d = a[k] - b[k] - mean;
            ss += d * d;
        }
        return Math.Sqrt(ss / n);
    }

    public static double LogSumExp(ReadOnlySpan<double> values)
    {
        double max = double.NegativeInfinity;
        foreach (double v in values) max = Math.Max(max, v);
        if (double.IsNegativeInfinity(max)) return max;
        double s = 0.0;
        foreach (double v in values) s += Math.Exp(v - max);
        return max + Math.Log(s);
    }
}
