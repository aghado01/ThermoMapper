namespace Maths.Regression.Spline;

/// <summary>
/// The active-span geometry of a spline design, computed once and reused. Each row of a B-spline design has only
/// degree+1 contiguous non-zero basis functions, so every observation model's normal equations <c>A = ZᵀW Z</c>
/// are banded; this records each row's first/last non-zero column and the global half-bandwidth, then accumulates
/// <c>A</c> straight into LAPACK lower-band storage (and the right-hand side <c>ZᵀW r</c>) — the dense ν×ν matrix
/// is never materialized, and the bandwidth feeds <see cref="Maths.LinAlg.BandCholesky"/> directly. The spans
/// depend only on the (fixed) design, so an IRLS loop computes them once and re-bands cheaply each reweighting.
/// </summary>
public sealed class BandedDesign
{
    private readonly double[,] _design;
    private readonly int[] _lo;   // first non-zero column per row, or −1 for an all-zero row
    private readonly int[] _hi;   // last non-zero column per row

    /// <summary>Number of design rows (observations).</summary>
    public int Rows { get; }

    /// <summary>Number of basis functions ν (design columns).</summary>
    public int Dimension { get; }

    /// <summary>Lower half-bandwidth of <c>A = ZᵀW Z</c> — the spline degree, or 0 for a step (indicator) basis.</summary>
    public int Bandwidth { get; }

    public BandedDesign(double[,] design)
    {
        _design = design;
        int m = design.GetLength(0);
        int nu = design.GetLength(1);
        Rows = m;
        Dimension = nu;
        _lo = new int[m];
        _hi = new int[m];

        int bw = 0;
        for (int i = 0; i < m; i++)
        {
            int l = 0;
            while (l < nu && design[i, l] == 0.0) l++;
            if (l == nu) { _lo[i] = -1; continue; }   // degenerate all-zero row
            int h = nu - 1;
            while (design[i, h] == 0.0) h--;
            _lo[i] = l;
            _hi[i] = h;
            if (h - l > bw) bw = h - l;
        }
        Bandwidth = bw;
    }

    /// <summary>
    /// Accumulates the Gram matrix <c>A = ZᵀW Z</c> into LAPACK lower-band storage
    /// <paramref name="band"/><c>[d, q] = A(q+d, q)</c> (sized <c>[Bandwidth+1, Dimension]</c>, pre-zeroed) and,
    /// when <paramref name="b"/>/<paramref name="r"/> are supplied, the right-hand side <c>b = ZᵀW r</c> in one
    /// pass. <paramref name="w"/> null means unit weights. Pass <paramref name="r"/>=null to build the Gram only
    /// (the constant-weight Whittle case factors A once and iterates only the right-hand side).
    /// </summary>
    public void Accumulate(double[]? w, double[,] band, double[]? r, double[]? b)
    {
        bool withRhs = r is not null && b is not null;
        for (int i = 0; i < Rows; i++)
        {
            int lo = _lo[i];
            if (lo < 0) continue;
            int hi = _hi[i];
            double wi = w is null ? 1.0 : w[i];
            double ri = withRhs ? r![i] : 0.0;
            for (int p = lo; p <= hi; p++)
            {
                double wzip = wi * _design[i, p];
                if (withRhs) b![p] += wzip * ri;
                for (int q = lo; q <= p; q++)
                    band[p - q, q] += wzip * _design[i, q];
            }
        }
    }

    /// <summary>
    /// Accumulates only the right-hand side <c>b = ZᵀW r</c> (pre-zeroed), for IRLS iterations that reuse a
    /// once-factored constant Gram (Whittle, W = 1). <paramref name="w"/> null means unit weights.
    /// </summary>
    public void AccumulateRhs(double[]? w, double[] r, double[] b)
    {
        for (int i = 0; i < Rows; i++)
        {
            int lo = _lo[i];
            if (lo < 0) continue;
            int hi = _hi[i];
            double wr = (w is null ? 1.0 : w[i]) * r[i];
            for (int p = lo; p <= hi; p++)
                b[p] += _design[i, p] * wr;
        }
    }
}
