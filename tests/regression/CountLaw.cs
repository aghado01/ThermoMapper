using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Tests;

/// <summary>
/// The constant-likelihood count-law check. With the likelihood held constant, a reversible-jump sampler's
/// stationary dimension law must be its count prior exactly; an error in the location measure, a Hastings factor
/// or the move-selection ratio shows up as a different law. Generic over the target, move set and carrier, so
/// every trans-dimensional engine (BARS, and BASS's D01 check) runs the same oracle. The standard errors are batch
/// means over one long chain.
/// </summary>
internal static class CountLaw
{
    /// <summary>
    /// Observed and expected frequencies for dimensions 0..<c>bins − 2</c>, the last bin pooling everything above,
    /// with the batch-means standard error of each observed frequency.
    /// </summary>
    public sealed record Result(double[] Observed, double[] Expected, double[] StandardError)
    {
        /// <summary>True when every bin satisfies |observed − expected| ≤ z·SE + floor.</summary>
        public bool Agrees(double z, double floor)
        {
            for (int j = 0; j < Observed.Length; j++)
                if (Math.Abs(Observed[j] - Expected[j]) > z * StandardError[j] + floor) return false;
            return true;
        }

        public string Describe()
        {
            var sb = new StringBuilder();
            for (int j = 0; j < Observed.Length; j++)
            {
                string label = j == Observed.Length - 1 ? $"≥{j}" : j.ToString(CultureInfo.InvariantCulture);
                sb.Append(CultureInfo.InvariantCulture, $" k={label}: {Observed[j]:F4} vs {Expected[j]:F4} (se {StandardError[j]:F4});");
            }
            return sb.ToString();
        }
    }

    /// <param name="moves">The move set, as the engine will run it.</param>
    /// <param name="target">A target whose likelihood is constant, so only the prior and the measure remain.</param>
    /// <param name="start">The chain's start.</param>
    /// <param name="dimension">The dimension of a state (the knot or basis-function count).</param>
    /// <param name="logPrior">The intended count prior, up to a constant; −∞ off its support.</param>
    /// <param name="bins">Number of reported bins, the last pooling the tail.</param>
    /// <param name="burn">Steps discarded before counting.</param>
    /// <param name="steps">Steps counted.</param>
    /// <param name="batches">Batches for the standard errors; must divide <paramref name="steps"/>.</param>
    /// <param name="seed">The chain's seed.</param>
    public static Result Sample<TState>(
        IReadOnlyList<IRjMove<TState>> moves, IRjTarget<TState> target, TState start, Func<TState, int> dimension,
        Func<int, double> logPrior, int bins, int burn, int steps, int batches, int seed)
    {
        if (bins < 2) throw new ArgumentOutOfRangeException(nameof(bins));
        if (batches < 2 || steps % batches != 0) throw new ArgumentOutOfRangeException(nameof(batches));

        double[] expected = ExpectedLaw(logPrior, bins);

        var chain = new ReversibleJumpChain<TState>(moves, target, start, new Xoshiro256PlusPlus(seed));
        for (int i = 0; i < burn; i++) chain.Step();

        int batchLength = steps / batches;
        var batchFreq = new double[batches, bins];
        for (int b = 0; b < batches; b++)
        {
            for (int i = 0; i < batchLength; i++)
                batchFreq[b, Math.Min(dimension(chain.Step()), bins - 1)] += 1.0;
            for (int j = 0; j < bins; j++) batchFreq[b, j] /= batchLength;
        }

        var observed = new double[bins];
        var se = new double[bins];
        for (int j = 0; j < bins; j++)
        {
            double mean = 0.0;
            for (int b = 0; b < batches; b++) mean += batchFreq[b, j];
            mean /= batches;
            double ss = 0.0;
            for (int b = 0; b < batches; b++) ss += (batchFreq[b, j] - mean) * (batchFreq[b, j] - mean);
            observed[j] = mean;
            se[j] = Math.Sqrt(ss / (batches - 1) / batches);
        }
        return new Result(observed, expected, se);
    }

    // The prior normalized over its support. Terms are summed until they fall 60 nats below the largest seen,
    // past the reported bins; a prior that never decays that far is not proper and fails here.
    private static double[] ExpectedLaw(Func<int, double> logPrior, int bins)
    {
        var logs = new List<double>();
        double max = double.NegativeInfinity;
        for (int k = 0; ; k++)
        {
            double lp = logPrior(k);
            logs.Add(lp);
            if (lp > max) max = lp;
            if (k >= bins && (double.IsNegativeInfinity(lp) || lp < max - 60.0)) break;
            if (k > 100_000) throw new InvalidOperationException("The count prior does not decay; it is not proper.");
        }

        var expected = new double[bins];
        double total = 0.0;
        for (int k = 0; k < logs.Count; k++)
        {
            double p = Math.Exp(logs[k] - max);
            total += p;
            expected[Math.Min(k, bins - 1)] += p;
        }
        for (int j = 0; j < bins; j++) expected[j] /= total;
        return expected;
    }
}
