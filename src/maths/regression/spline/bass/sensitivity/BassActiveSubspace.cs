using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Subspaces;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// E8's readout (<see cref="BassFit.ActiveSubspace"/>, the active-subspace brief, slice 2): C per draw on one of two routes, read in
/// a declared metric, then the posterior's summaries (<see cref="SubspacePosterior"/>) and decision 6's checks.
/// <list type="bullet">
/// <item><b>Exact</b> (decision 1): the moment service's derivative moments under the product law, for terms that have factors, which
/// every term type built today has. <see cref="ActiveSubspaceRoute.Auto"/> takes it.</item>
/// <item><b>Sampled:</b> the gradient readout (<see cref="BassFit.Gradient"/>) at N points drawn once from the law and shared by every
/// draw, so the draws' spread is the posterior's and not the Monte Carlo error's (decision 5). Each draw's subspace comes from the thin
/// SVD of its N × p gradient matrix over √N, not from C, and each entry of C carries its standard error.</item>
/// </list>
/// The metric (decision 2) is a map L from metric to box coordinates, C = LᵀC_boxL. Under a product law L is diagonal: 1 for the box,
/// one over the box's width for raw units, the law's standard deviation in box units for whitening.
/// </summary>
internal static class BassActiveSubspaces
{
    private const int Batch = 64;   // draws read from the gradient readout at once

    public static ActiveSubspaceAnalysis Of(BassFit fit, BassData data, ActiveSubspaceSpec spec)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(spec);
        string fingerprint = Bass.Fingerprint(fit.Config, data);
        if (fingerprint != fit.Identity.Fingerprint)
            throw new ArgumentException($"The fit was run on other data: its fingerprint is {fit.Identity.Fingerprint}, theirs {fingerprint}.", nameof(data));
        if (spec.Points < 2) throw new ArgumentOutOfRangeException(nameof(spec), "The sampled route needs at least two points.");

        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        if (inputs.Outputs is not null)
            throw new ArgumentException("A functional fit's C takes a declared output metric, E8's slice 3; this readout is a scalar response's.", nameof(fit));
        var measure = ProductMeasure.Resolve(fit.Transform, inputs, spec.Measures);
        int p = inputs.Cases.ContinuousCount;
        if (p == 0) throw new ArgumentException("C needs a continuous predictor: a categorical one has no derivative.", nameof(fit));
        if (spec.Dimension < 1 || spec.Dimension > p)
            throw new ArgumentOutOfRangeException(nameof(spec), $"The dimension k is from 1 to the {p} continuous predictors.");

        int[] variables = [.. Enumerable.Range(0, p).Select(measure.OfContinuous)];
        string[] names = [.. variables.Select(v => measure.Names[v])];
        string[] excluded = [.. Enumerable.Range(0, inputs.Cases.CategoricalCount).Select(c => measure.Names[measure.OfCategorical(c)])];
        double[] widths = [.. variables.Select(v => { var s = (RangeScale)fit.Transform.Predictors[v]; return s.Max - s.Min; })];
        double[] scale = Scale(spec.Metric, measure, variables, widths);
        var map = new double[p, p];
        for (int i = 0; i < p; i++) map[i, i] = scale[i];

        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        (int[] selected, double[] weights) = Draws3.Select(spec.Draws, pooled.Length, fit.DrawWeights, nameof(spec));
        ActiveSubspaceRoute route = spec.Route == ActiveSubspaceRoute.Auto ? ActiveSubspaceRoute.Exact : spec.Route;

        var cBox = new double[selected.Length][,];
        var draws = new DrawSubspace[selected.Length];
        double[][,]? error = null;
        if (route == ActiveSubspaceRoute.Exact)
        {
            var moments = BassMoments.On(inputs, fit.Config.Structure.Degree, measure);
            for (int d = 0; d < selected.Length; d++)
            {
                BassDraw draw = pooled[selected[d]];
                BassTerm[] terms = [.. draw.Structure.TermIds.Select(id => fit.Terms[id])];
                cBox[d] = moments.GradientOuterProduct(terms, draw.Coefficients.AsSpan(1));
                draws[d] = DrawSubspace.FromC(Congruence(cBox[d], scale), spec.Dimension);
            }
        }
        else
        {
            error = new double[selected.Length][,];
            (PredictionCases cases, int n) = Points(fit, measure, variables, spec);
            for (int start = 0; start < selected.Length; start += Batch)
            {
                int[] batch = selected[start..Math.Min(selected.Length, start + Batch)];
                BassGradient g = fit.Gradient(data, cases, batch);
                for (int b = 0; b < batch.Length; b++)
                {
                    int d = start + b;
                    // The draw's gradients in box units (∂f/∂x = width·∂f/∂z), then in the metric, over √N.
                    var factor = new double[n, p];
                    var box = new double[n, p];
                    for (int i = 0; i < p; i++)
                        for (int r = 0; r < n; r++)
                        {
                            box[r, i] = widths[i] * g.Partials[i][b, r, 0];
                            factor[r, i] = box[r, i] * scale[i] / Math.Sqrt(n);
                        }
                    draws[d] = DrawSubspace.FromFactor(factor, spec.Dimension);
                    (cBox[d], double[,] se) = MeanAndError(box);
                    error[d] = Congruence(se, scale);
                }
            }
        }
        var posterior = SubspacePosterior.Of(draws, weights);
        PoincareCheck poincare = Poincare(fit, data, spec, measure, variables, selected, cBox, error is null ? null : error, scale);
        return new ActiveSubspaceAnalysis(Array.AsReadOnly(names), Array.AsReadOnly(excluded), route, spec.Metric, map, measure,
                                          route == ActiveSubspaceRoute.Sampled ? spec.Points : 0, posterior, error, poincare, Array.AsReadOnly(selected));
    }

    // L's diagonal: box coordinates per metric unit.
    private static double[] Scale(InputMetric metric, ProductMeasure measure, int[] variables, double[] widths) => metric switch
    {
        InputMetric.Box => [.. variables.Select(_ => 1.0)],
        InputMetric.Raw => [.. widths.Select(w => 1.0 / w)],
        InputMetric.Whitened => [.. variables.Select(v => StandardDeviation(measure.Continuous(v)!))],
        _ => throw new ArgumentOutOfRangeException(nameof(metric), $"No metric {metric}."),
    };

    // The law's standard deviation in box units, from its interval moments about the support's centre.
    private static double StandardDeviation(ContinuousMeasure m)
    {
        Span<double> moments = stackalloc double[3];
        double centre = 0.5 * (m.Lower + m.Upper);
        m.Measure.IntervalMoments(Math.BitDecrement(m.Lower), m.Upper, centre, 0.0, moments);
        double mean = moments[1] / moments[0], variance = moments[2] / moments[0] - mean * mean;
        if (!(variance > 0.0)) throw new ArgumentException("Whitening needs each input's law to have a positive variance.");
        return Math.Sqrt(variance);
    }

    // LᵀCL for a diagonal L.
    private static double[,] Congruence(double[,] c, double[] scale)
    {
        int p = scale.Length;
        var z = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) z[i, j] = scale[i] * c[i, j] * scale[j];
        return z;
    }

    // The mean over the rows of g gᵀ, and each entry's standard error.
    private static (double[,] Mean, double[,] Error) MeanAndError(double[,] g)
    {
        int n = g.GetLength(0), p = g.GetLength(1);
        var mean = new double[p, p];
        var error = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = i; j < p; j++)
            {
                double sum = 0.0, square = 0.0;
                for (int r = 0; r < n; r++)
                {
                    double v = g[r, i] * g[r, j];
                    sum += v;
                    square += v * v;
                }
                double m = sum / n, variance = Math.Max(0.0, (square - n * m * m) / (n - 1));
                mean[i, j] = mean[j, i] = m;
                error[i, j] = error[j, i] = Math.Sqrt(variance / n);
            }
        return (mean, error);
    }

    // N points from the product law, as cases in the inputs' own units, drawn once from the spec's seed: each continuous input from
    // its law in box units (an empirical marginal by its atoms), mapped back through the box; each categorical one by its level
    // probabilities.
    private static (PredictionCases Cases, int N) Points(BassFit fit, ProductMeasure measure, int[] variables, ActiveSubspaceSpec spec)
    {
        int n = spec.Points;
        var rng = new Xoshiro256PlusPlus(spec.Seed);
        var columns = new List<PredictionColumn>();
        for (int v = 0; v < measure.Count; v++)
        {
            if (measure.OutputOf(v) >= 0) continue;
            if (measure.Continuous(v) is { } m)
            {
                var s = (RangeScale)fit.Transform.Predictors[v];
                var values = new double[n];
                for (int r = 0; r < n; r++)
                {
                    double x = m.Measure is Empirical e ? e.Values[rng.NextInt(e.Count)] : ((IContinuousDistribution)m.Measure).Sample(rng);
                    values[r] = s.Min + x * (s.Max - s.Min);
                }
                columns.Add(new ContinuousColumn(measure.Names[v], values));
            }
            else
            {
                IReadOnlyList<double> probabilities = measure.Probabilities(v)!;
                var levels = ((TrainingLevels)fit.Transform.Predictors[v]).Levels;
                var picks = new string[n];
                for (int r = 0; r < n; r++)
                {
                    double u = rng.NextDouble(), cumulative = 0.0;
                    int l = 0;
                    for (; l < probabilities.Count - 1; l++)
                    {
                        cumulative += probabilities[l];
                        if (u < cumulative) break;
                    }
                    while (probabilities[l] == 0.0 && l > 0) l--;
                    picks[r] = levels[l];
                }
                columns.Add(new CategoricalColumn(measure.Names[v], picks));
            }
        }
        return (new PredictionCases(columns), n);
    }

    // Decision 6: Poincaré per draw and continuous input with a uniform law, and SG-4(b) on a Legendre fit under the product uniform
    // law on the boxes; the totals are the moment route's, under the same law and draws.
    private static PoincareCheck Poincare(BassFit fit, BassData data, ActiveSubspaceSpec spec, ProductMeasure measure, int[] variables,
                                          int[] selected, double[][,] cBox, double[][,]? error, double[] scale)
    {
        int draws = selected.Length, p = variables.Length;
        var total = new double[draws, p];
        var bound = new double[draws, p];
        var ratio = new double[draws, p];
        SobolDecomposition sobol = BassSensitivity.Of(fit, data, new SensitivitySpec { Measures = spec.Measures, Draws = selected, MaxOrder = 1 }).Decompositions[0];
        double[,]? jacobi = null;
        if (LegendreOnBoxes(fit, measure, variables, selected))
        {
            double[,] all = BassSensitivity.FromCoefficients(fit, data).Jacobi;
            jacobi = new double[draws, p];
            for (int d = 0; d < draws; d++)
                for (int i = 0; i < p; i++) jacobi[d, i] = all[selected[d], variables[i]];
        }
        double[,]? jacobiRatio = jacobi is null ? null : new double[draws, p];
        int violations = 0;
        for (int i = 0; i < p; i++)
        {
            int v = variables[i];
            bool uniform = measure.Declared[v] is UniformMeasure;
            ContinuousMeasure m = measure.Continuous(v)!;
            double constant = uniform ? (m.Upper - m.Lower) * (m.Upper - m.Lower) / (Math.PI * Math.PI) : double.NaN;
            for (int d = 0; d < draws; d++)
            {
                double t = sobol.TotalVariance[d, v];
                total[d, i] = t;
                bound[d, i] = constant * cBox[d][i, i];
                ratio[d, i] = uniform ? (bound[d, i] > 0.0 ? t / bound[d, i] : (t > 0.0 ? double.PositiveInfinity : 0.0)) : double.NaN;
                if (uniform)
                {
                    // Rounding on the exact route; three standard errors of the bound on the sampled route (the error is in the metric).
                    double slack = error is null ? 1e-9 * Math.Max(Math.Abs(t), Math.Abs(bound[d, i])) : 3.0 * constant * error[d][i, i] / (scale[i] * scale[i]);
                    if (t > bound[d, i] + slack) violations++;
                }
                if (jacobi is not null)
                {
                    double half = 0.5 * jacobi[d, i];
                    jacobiRatio![d, i] = half > 0.0 ? t / half : (t > 0.0 ? double.PositiveInfinity : 0.0);
                    if (t > half + 1e-9 * Math.Max(Math.Abs(t), Math.Abs(half))) violations++;
                }
            }
        }
        return new PoincareCheck(total, bound, ratio, jacobiRatio, violations);
    }

    // Whether SG-4(b)'s closed form applies: every term of the draws a Legendre product, and every input's law uniform on its box.
    private static bool LegendreOnBoxes(BassFit fit, ProductMeasure measure, int[] variables, int[] selected)
    {
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        foreach (int s in selected)
            foreach (int id in pooled[s].Structure.TermIds)
            {
                BassTerm term = fit.Terms[id];
                if (term.Legendre.Count == 0 || term.Continuous.Count > 0 || term.Categorical.Count > 0 || term.Functional.Count > 0) return false;
            }
        if (measure.Count != variables.Length) return false;
        foreach (int v in variables)
        {
            ContinuousMeasure m = measure.Continuous(v)!;
            if (measure.Declared[v] is not UniformMeasure || m.Lower != 0.0 || m.Upper != 1.0) return false;   // the box, as FromCoefficients reads it
        }
        return true;
    }
}
