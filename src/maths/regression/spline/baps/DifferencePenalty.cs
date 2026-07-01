using System;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// r-th order difference penalty for penalized (P-spline) regression — the measure-side smoother dual to
/// free-knot's metric-side knot adaptation (Eilers &amp; Marx 1996; the penalized carrier of MRA2015's BAPS). On a
/// rich fixed B-spline basis the roughness penalty is <c>βᵀ DᵀΛD β</c> with D the r-th finite-difference operator;
/// its null space (polynomials of degree &lt; r) is left unpenalized, so λ→∞ drives the fit to that polynomial and
/// λ→0 to the unpenalized least-squares fit. <c>DᵀΛD</c> is banded at half-bandwidth r, so it adds straight into
/// the already-banded normal equations and the whole penalized system stays a <see cref="Maths.LinAlg.BandCholesky"/>
/// solve (the Reinsch system). Λ is a per-difference weight: uniform Λ = λI is global smoothing, a varying Λ is the
/// locally-adaptive penalty.
/// </summary>
public sealed class DifferencePenalty : IBandPenalty
{
    // Difference-operator row coefficients c_l = (−1)^l · C(r, l), l = 0..r (row i couples columns i..i+r).
    private readonly double[] _c;

    /// <summary>Difference order r (2 = the usual P-spline curvature penalty).</summary>
    public int Order { get; }

    /// <param name="order">Difference order r ≥ 1; defaults to 2.</param>
    public DifferencePenalty(int order = 2)
    {
        if (order < 1) throw new ArgumentOutOfRangeException(nameof(order), "Penalty order must be ≥ 1.");
        Order = order;
        _c = new double[order + 1];
        double binom = 1.0;                       // C(r, 0)
        for (int l = 0; l <= order; l++)
        {
            _c[l] = (l % 2 == 0 ? 1.0 : -1.0) * binom;
            binom = binom * (order - l) / (l + 1); // C(r, l+1)
        }
    }

    /// <summary>Penalty half-bandwidth (= <see cref="Order"/>).</summary>
    public int Bandwidth => Order;

    /// <summary>Null-space dimension (= <see cref="Order"/>): the unpenalized degree-&lt;r polynomial coefficients.</summary>
    public int Nullity => Order;

    /// <summary>
    /// Adds the global penalty <c>λ·DᵀD</c> into LAPACK lower-band storage <paramref name="band"/><c>[d, j] =
    /// A(j+d, j)</c> over <paramref name="dim"/> coefficients. <paramref name="band"/>'s first dimension must be
    /// ≥ Order+1; typically the band already holds the Gram <c>ZᵀWZ</c> and this adds the roughness term in place.
    /// </summary>
    public void AccumulateInto(double[,] band, int dim, double lambda)
        => AccumulateStrided(band, dim, stride: 1, offset: 0, lambda);

    /// <summary>
    /// Adds <c>λ·DᵀD</c> over a sequence of <paramref name="count"/> coefficients placed at flattened positions
    /// <c>offset + i·stride</c>. The tensor-product penalty <c>P_x⊗I + I⊗P_y</c> is two such strided 1-D penalties
    /// (stride = inner-dimension size for the outer factor, stride 1 for the inner) — see <see cref="TensorPenalty"/>.
    /// </summary>
    public void AccumulateStrided(double[,] band, int count, int stride, int offset, double lambda)
    {
        for (int i = 0; i + Order < count; i++)       // one row per r-th difference
            for (int a = 0; a <= Order; a++)
            {
                double la = lambda * _c[a];
                for (int b = 0; b <= a; b++)          // lower triangle: row ≥ col
                    band[(a - b) * stride, offset + (i + b) * stride] += la * _c[b];
            }
    }

    /// <summary>
    /// Adds the locally-adaptive penalty <c>Dᵀ diag(weights) D</c> into band storage. <paramref name="weights"/>
    /// has one entry per difference row (length <paramref name="dim"/>−Order) — the local smoothing strength λ_i.
    /// </summary>
    public void AccumulateInto(double[,] band, int dim, double[] weights)
    {
        ArgumentNullException.ThrowIfNull(weights);
        if (weights.Length != dim - Order)
            throw new ArgumentException($"Expected {dim - Order} penalty weights, got {weights.Length}.", nameof(weights));
        for (int i = 0; i + Order < dim; i++)
        {
            double wi = weights[i];
            for (int a = 0; a <= Order; a++)
            {
                double wa = wi * _c[a];
                for (int b = 0; b <= a; b++)
                    band[a - b, i + b] += wa * _c[b];
            }
        }
    }

    /// <summary>Number of r-th differences over <paramref name="dim"/> coefficients (= the local-weight count, ν−r).</summary>
    public int DifferenceCount(int dim) => Math.Max(0, dim - Order);

    /// <summary>
    /// The roughness quadratic <c>βᵀDᵀD β = ‖Dβ‖² = Σ_i (Δ^r β)_i²</c> — the penalized sum of squares that drives
    /// the smoothing-variance update in a Bayesian P-spline (the τ² full conditional's scale).
    /// </summary>
    public double Roughness(double[] beta)
    {
        ArgumentNullException.ThrowIfNull(beta);
        double sum = 0.0;
        for (int i = 0; i + Order < beta.Length; i++)
        {
            double d = 0.0;
            for (int l = 0; l <= Order; l++) d += _c[l] * beta[i + l];   // (Δ^r β)_i
            sum += d * d;
        }
        return sum;
    }

    /// <summary>
    /// Writes the per-difference squared values <c>into[i] = (Δ^r β)_i²</c> (length <c>DifferenceCount</c>) — the
    /// local roughnesses that drive each adaptive multiplier's conjugate update in locally-adaptive BAPS.
    /// </summary>
    public void SquaredDifferencesInto(double[] beta, double[] into)
    {
        ArgumentNullException.ThrowIfNull(beta);
        ArgumentNullException.ThrowIfNull(into);
        for (int i = 0; i + Order < beta.Length; i++)
        {
            double d = 0.0;
            for (int l = 0; l <= Order; l++) d += _c[l] * beta[i + l];
            into[i] = d * d;
        }
    }

    /// <summary>
    /// The penalty matrix <c>DᵀD</c> as a dense symmetric <paramref name="dim"/>×<paramref name="dim"/> matrix —
    /// for the 1-D eigendecomposition behind anisotropic tensor REML (the combined penalty's eigenvalues are the
    /// pairwise sums of the two factors' spectra). Small ν only; the banded form is the workhorse elsewhere.
    /// </summary>
    public double[,] ToDense(int dim)
    {
        var band = new double[Order + 1, dim];
        AccumulateInto(band, dim, 1.0);
        var m = new double[dim, dim];
        for (int j = 0; j < dim; j++)
            for (int d = 0; d <= Order && j + d < dim; d++)
            {
                m[j + d, j] = band[d, j];
                m[j, j + d] = band[d, j];
            }
        return m;
    }
}
