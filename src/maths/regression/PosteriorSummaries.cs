using System;
using System.Collections.Generic;

namespace Maths.Regression;

/// <summary>A posterior median with its 95% interval: type-7 quantiles of the pooled draws.</summary>
public sealed record PosteriorInterval(double Median, double Lo, double Hi)
{
    /// <summary>The interval of a pooled sample, which is sorted in place.</summary>
    internal static PosteriorInterval Of(double[] sample)
    {
        Array.Sort(sample);
        return new PosteriorInterval(Quantile(sample, 0.5), Quantile(sample, 0.025), Quantile(sample, 0.975));
    }

    /// <summary>The type-7 quantile of a sorted sample.</summary>
    internal static double Quantile(double[] sorted, double p)
    {
        double pos = p * (sorted.Length - 1);
        int lo = (int)Math.Floor(pos);
        int hi = Math.Min(lo + 1, sorted.Length - 1);
        return sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]);
    }
}

/// <summary>A pointwise band on a fitted curve at the rows of a design: the posterior mean and sd, and the type-7 2.5%
/// and 97.5% quantiles of the draws' fits.</summary>
public sealed record PointwiseBand(double[] Mean, double[] Sd, double[] Lo, double[] Hi)
{
    /// <summary>
    /// The band at the rows of <paramref name="design"/>: the mean is the design times <paramref name="meanCoefficients"/>
    /// (the reported, affine readout), and the sd and quantiles are those of the fits design × β over the pooled
    /// coefficient draws.
    /// </summary>
    internal static PointwiseBand Of(double[,] design, double[] meanCoefficients, IReadOnlyList<double[]> draws)
    {
        ArgumentNullException.ThrowIfNull(design);
        int m = design.GetLength(0), nu = meanCoefficients.Length, count = draws.Count;
        if (design.GetLength(1) != nu)
            throw new ArgumentException($"The design needs {nu} columns, one per coefficient.", nameof(design));
        var fits = new double[m][];
        for (int i = 0; i < m; i++) fits[i] = new double[count];
        for (int s = 0; s < count; s++)
        {
            double[] beta = draws[s];
            for (int i = 0; i < m; i++)
            {
                double f = 0.0;
                for (int j = 0; j < nu; j++) f += design[i, j] * beta[j];
                fits[i][s] = f;
            }
        }
        var mean = new double[m];
        var sd = new double[m];
        var lo = new double[m];
        var hi = new double[m];
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < nu; j++) mean[i] += design[i, j] * meanCoefficients[j];
            double drawMean = 0.0;
            foreach (double f in fits[i]) drawMean += f;
            drawMean /= count;
            double ss = 0.0;
            foreach (double f in fits[i]) ss += (f - drawMean) * (f - drawMean);
            sd[i] = Math.Sqrt(ss / count);
            Array.Sort(fits[i]);
            lo[i] = PosteriorInterval.Quantile(fits[i], 0.025);
            hi[i] = PosteriorInterval.Quantile(fits[i], 0.975);
        }
        return new PointwiseBand(mean, sd, lo, hi);
    }
}
