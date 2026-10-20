using System;
using System.Collections.Generic;

namespace Maths.Regression.Bgp;

/// <summary>How <see cref="GpRegression"/> maps its inputs into the frame that the kernel and the bandwidth prior see.</summary>
public enum GpInputScaling
{
    /// <summary>
    /// Isotropic similarity <c>x′ = (x − lo)/s</c> fitted on the training inputs: lo is the per-coordinate minimum and s
    /// the largest coordinate span, so the data fill the unit cube with their shape preserved. The default.
    /// </summary>
    Similarity,

    /// <summary>Inputs used as given: the caller supplies coordinates whose squared distances are of order 1.</summary>
    None,
}

/// <summary>
/// The fitted input map shared by training and prediction: <c>x′ = (x − InputOffset)/InputScale</c>. Bandwidths live in
/// the mapped frame; a bandwidth t there is <c>InputScale²·t</c> in the original units. The
/// map is fitted on the sample, so its constants depend on n (the qualification recorded under BGP-2 of the
/// TWCD2025 lemma notes).
/// </summary>
public sealed class GpTransform
{
    private readonly double[] _offset;

    private GpTransform(GpInputScaling scaling, double[] offset, double scale)
    {
        Scaling = scaling;
        _offset = offset;
        InputScale = scale;
    }

    /// <summary>Which input map was fitted.</summary>
    public GpInputScaling Scaling { get; }

    /// <summary>Per-coordinate offset subtracted from the inputs (zero for <see cref="GpInputScaling.None"/>).</summary>
    public IReadOnlyList<double> InputOffset => _offset;

    /// <summary>The single scale s dividing every coordinate (1 for <see cref="GpInputScaling.None"/>).</summary>
    public double InputScale { get; }

    /// <summary>Fit the map on the training inputs <paramref name="x"/> (n × D).</summary>
    public static GpTransform Fit(double[,] x, GpInputScaling scaling)
    {
        ArgumentNullException.ThrowIfNull(x);
        int n = x.GetLength(0), d = x.GetLength(1);
        if (n < 1) throw new ArgumentException("Need at least one input.", nameof(x));
        foreach (double v in x)
            if (!double.IsFinite(v)) throw new ArgumentException("Inputs must be finite.", nameof(x));

        var offset = new double[d];
        double scale = 1.0;
        if (scaling == GpInputScaling.Similarity)
        {
            scale = 0.0;
            for (int k = 0; k < d; k++)
            {
                double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
                for (int i = 0; i < n; i++) { lo = Math.Min(lo, x[i, k]); hi = Math.Max(hi, x[i, k]); }
                offset[k] = lo;
                scale = Math.Max(scale, hi - lo);
            }
            if (!(scale > 0.0) || !double.IsFinite(scale))
                throw new ArgumentException("Similarity scaling needs inputs that are not all identical and a finite span.", nameof(x));
        }
        else if (scaling != GpInputScaling.None)
        {
            throw new ArgumentOutOfRangeException(nameof(scaling));
        }
        return new GpTransform(scaling, offset, scale);
    }

    /// <summary>Map inputs (m × D, the training dimension) into the fitted frame; returns a new array.</summary>
    public double[,] MapInputs(double[,] x)
    {
        ArgumentNullException.ThrowIfNull(x);
        int m = x.GetLength(0), d = x.GetLength(1);
        if (d != _offset.Length) throw new ArgumentException("Inputs must have the training dimension.", nameof(x));
        var mapped = new double[m, d];
        for (int i = 0; i < m; i++)
            for (int k = 0; k < d; k++) mapped[i, k] = (x[i, k] - _offset[k]) / InputScale;
        return mapped;
    }
}
