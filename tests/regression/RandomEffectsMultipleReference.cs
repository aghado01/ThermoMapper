using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;

namespace Maths.Regression.Tests;

// Small finite structures and a full case-covariance oracle, independent of the production whitening.
internal static class RandomEffectsMultipleReference
{
    internal static (BassConfig Config, BassData Data) Fixture(string kind, bool drawn = true, bool tauDrawn = false)
    {
        var example = RandomEffectsWhiteningTests.Example(kind);
        int n = example.Groups[0].Length;
        double[] x = [.. Enumerable.Range(0, n).Select(i => ((7 * i) % n + 0.25) / n)];
        double[] y = [.. x.Select((v, i) => 0.6 + 0.5 * (2 * v - 1) + 0.45 * Math.Sin(1.7 * i) + 0.3 * Math.Cos(0.6 * i))];
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxBasis = 1, MaxInteraction = 1, Legendre = new LegendreBasis { MaxDegree = 2 } },
            Coefficients = BassCoefficientPrior.RawRidge(tauDrawn ? TauPrior.Gamma(2, 2) : TauPrior.Held(0.7)),
            Count = new BassCountPrior(1.5, 1.0), Intercept = GPriorCentre.Zero, Sigma2 = Sigma2Prior.InverseGamma(3, 2),
            RandomIntercepts = [.. example.Groups.Select((g, k) => new BassRandomIntercept(g.Name, example.Scales[k],
                example.Parents[k] is null ? RandomInterceptConvention.ObservationCentred : RandomInterceptConvention.WithinParent,
                drawn ? new BassGroupScalePrior(k == 0 ? 0.8 : 1.2) : null, example.Parents[k]?.Name))],
            Sampling = new BassSampling { Chains = 1, BurnIn = 1000, Samples = 1000, MasterSeed = 21029 },
        };
        return (config, new BassData([new ContinuousPredictor("x", x, new InputDomain(0, 1))], new ScalarResponse(y), example.Groups));
    }

    internal static BassTerm[][] Structures(BassModel model) => [[], .. BassTermSpace.Terms(model).Select(t => new[] { t })];

    internal static double[,] Covariance(BassModel model, ReadOnlySpan<double> scales)
    {
        var blocks = model.RandomEffects!;
        var parents = blocks.Specifications.Select(s => s.ParentGrouping is { } name ? model.Data.Groupings.Single(g => g.Name == name) : null).ToArray();
        var fixture = new RandomEffectsWhiteningTests.Fixture(blocks.Groupings.ToArray(), parents, scales.ToArray());
        return RandomEffectsWhiteningTests.Reference(fixture).Covariance;
    }

    internal static double DenseLogMarginal(BassModel model, BassTerm[] terms, ReadOnlySpan<double> scales, double tau)
    {
        double[,] covariance = Covariance(model, scales);
        int n = covariance.GetLength(0), p = terms.Length + 1;
        var columns = new double[p][];
        columns[0] = [.. Enumerable.Repeat(1.0, n)];
        for (int j = 1; j < p; j++)
        {
            columns[j] = new double[n];
            model.Basis.CaseColumn(terms[j - 1], model.Inputs.Cases, columns[j]);
        }
        // Raw ridge with a zero-centred intercept: coefficient covariance is I/τ.
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                foreach (double[] column in columns) covariance[i, j] += column[i] * column[j] / tau;
        var factor = new DenseCholesky(n);
        factor.Decompose(covariance);
        if (!factor.Status.PositiveDefinite) throw new InvalidOperationException("The dense reference covariance is not positive definite.");
        double[] y = model.Data.Response.Values.ToArray(), solution = factor.Solve(y);
        double quadratic = y.Zip(solution, (a, b) => a * b).Sum();
        const double shape = 3, rate = 2;
        return SpecialFunctions.LogGamma(shape + 0.5 * n) - SpecialFunctions.LogGamma(shape) + shape * Math.Log(rate)
            - (shape + 0.5 * n) * Math.Log(rate + 0.5 * quadratic) - 0.5 * n * Math.Log(2 * Math.PI) - 0.5 * factor.LogDet;
    }

    internal sealed record Law(double LogZ, double[] Structures, double[] LogScales);

    // Integrate each half-Cauchy by its own probability coordinate. The full (0,∞) range is covered,
    // without a cut-off in ln s; Gauss nodes avoid both endpoints. Refine the rule before using its result.
    internal static Law Enumerate(BassModel model, int nodes)
    {
        BassTerm[][] structures = Structures(model);
        double[] prior = [.. structures.Select(t => model.Prior.LogCollectionWeight(t.Length) + t.Sum(model.Prior.LogPrior))];
        double priorMass = prior.Select(Math.Exp).Sum();
        GaussLegendre rule = GaussLegendre.Of(nodes);
        double[] mass = new double[structures.Length], moments = new double[2];
        double top = double.NegativeInfinity, total = 0;
        for (int a = 0; a < nodes; a++)
            for (int b = 0; b < nodes; b++)
            {
                int[] indices = [a, b];
                double[] scales = [.. indices.Select((v, k) => model.RandomEffects!.Specifications[k].ScalePrior!.Scale
                    * Math.Tan(0.25 * Math.PI * (1 + rule.Nodes[v])))];
                double cell = Math.Log(rule.Weights[a] * rule.Weights[b] / 4);
                for (int j = 0; j < structures.Length; j++)
                {
                    double log = prior[j] + cell + DenseLogMarginal(model, structures[j], scales, model.Target.HeldTau!.Value);
                    if (log > top)
                    {
                        double rescale = Math.Exp(top - log);
                        total *= rescale;
                        for (int k = 0; k < mass.Length; k++) mass[k] *= rescale;
                        for (int k = 0; k < moments.Length; k++) moments[k] *= rescale;
                        top = log;
                    }
                    double weight = Math.Exp(log - top);
                    total += weight;
                    mass[j] += weight;
                    for (int k = 0; k < 2; k++) moments[k] += weight * Math.Log(scales[k]);
                }
            }
        return new Law(top + Math.Log(total) - Math.Log(priorMass), [.. mass.Select(v => v / total)], [.. moments.Select(v => v / total)]);
    }

    // Direct level-space prior simulation, projected separately inside each declared parent; no production contrast basis.
    internal static double[][] PriorEffects(BassModel model, double[] scales, double sigma2, Xoshiro256PlusPlus rng)
    {
        var result = new double[scales.Length][];
        for (int k = 0; k < result.Length; k++)
        {
            var grouping = model.RandomEffects!.Groupings[k];
            int[] codes = grouping.Codes.ToArray();
            int[] parent = model.RandomEffects.Specifications[k].ParentGrouping is { } name
                ? model.Data.Groupings.Single(g => g.Name == name).Codes.ToArray() : new int[codes.Length];
            int[] levels = [.. codes.Distinct().Order()];
            double[] counts = [.. levels.Select(g => (double)codes.Count(c => c == g))];
            int[] families = [.. levels.Select(g => parent[Array.FindIndex(codes, c => c == g)])];
            double[] z = [.. levels.Select(_ => Normal.SampleStandard(rng))];
            result[k] = new double[levels.Length];
            for (int g = 0; g < levels.Length; g++)
            {
                double norm = 0, dot = 0;
                for (int h = 0; h < levels.Length; h++)
                    if (families[h] == families[g]) { norm += counts[h] * counts[h]; dot += counts[h] * z[h]; }
                result[k][g] = Math.Sqrt(sigma2) * scales[k] * (z[g] - counts[g] * dot / norm);
            }
        }
        return result;
    }
}
