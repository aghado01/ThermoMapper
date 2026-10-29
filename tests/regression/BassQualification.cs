using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;

namespace Maths.Regression.Tests;

/// <summary>
/// What S10's qualification runs share: the switch that runs them (<c>BASS_QUALIFICATION</c> = 1; each is heavy, and the
/// fitting checklist records what it read), each draw's fitted values, pointwise summaries and a Gaussian draw.
/// </summary>
internal static class BassQualification
{
    /// <summary>Whether the qualification runs are asked for.</summary>
    internal static bool Enabled => Environment.GetEnvironmentVariable("BASS_QUALIFICATION") is "1";

    /// <summary>Each retained draw's f at every training observation, the draws pooled in chain order, from the terms'
    /// explicit uncentred columns.</summary>
    internal static double[][] Fitted(BassModel model, BassFit fit)
    {
        int n = model.Data.Response.Observations;
        var columns = new Dictionary<int, double[]>();
        var fitted = new List<double[]>();
        foreach (BassDraw draw in fit.Draws.SelectMany(c => c))
        {
            var f = new double[n];
            Array.Fill(f, draw.Coefficients[0]);
            int[] ids = draw.Structure.TermIds;
            for (int q = 0; q < ids.Length; q++)
            {
                if (!columns.TryGetValue(ids[q], out double[]? column)) columns[ids[q]] = column = BassPointwiseTests.Column(model, fit.Terms[ids[q]]);
                double c = draw.Coefficients[q + 1];
                for (int i = 0; i < n; i++) f[i] += c * column[i];
            }
            fitted.Add(f);
        }
        return [.. fitted];
    }

    /// <summary>Observation <paramref name="k"/>'s posterior mean and sd over <paramref name="fitted"/>, and its central
    /// interval at <paramref name="level"/> (type-7 quantiles).</summary>
    internal static (double Mean, double Sd, double Lower, double Upper) At(double[][] fitted, int k, double level = 0.95)
    {
        double[] values = [.. fitted.Select(f => f[k]).Order()];
        double mean = values.Average();
        double sd = Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1));
        return (mean, sd, PosteriorInterval.Quantile(values, (1 - level) / 2), PosteriorInterval.Quantile(values, (1 + level) / 2));
    }

    /// <summary>A standard Gaussian draw by Box–Muller, as the BASS fixtures draw their noise.</summary>
    internal static double Gaussian(Xoshiro256PlusPlus rng)
        => Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
}
