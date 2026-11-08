using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>How a random-intercept block keeps clear of the level-centred intercept (E21).</summary>
public enum RandomInterceptConvention
{
    /// <summary>
    /// The group effects' observation-weighted mean is zero, Σ_g n_g u_g = 0, so their contribution Zu has mean zero over the cases and
    /// the intercept alone carries the response's level. The unweighted sum-to-zero constraint is not offered: on unbalanced groups Zu
    /// keeps a nonzero mean, which is the overlap with the level-centred intercept the convention exists to remove.
    /// </summary>
    ObservationCentred,
}

/// <summary>
/// Random-effect intercepts on the data's grouping (E21; <see cref="BassData.Grouping"/>): y = α + Σⱼ βⱼBⱼ(x) + u_g + ε, with
/// u | σ² ~ N(0, σ²s²Π) on the observed levels, Π the projector onto {u : Σ_g n_g u_g = 0} (<see cref="Convention"/>), so
/// <paramref name="Scale"/> s is the group effects' standard deviation in units of σ, held. No term reads the grouping, and the coefficient
/// prior is the config's, of any form. The block is integrated out of every structure's marginal (<see cref="RandomInterceptBlock"/>), and
/// each retained draw carries its group effects (<see cref="BassDraw.GroupEffects"/>).
/// </summary>
/// <param name="Scale">s: the group effects' standard deviation over σ, positive and finite.</param>
/// <param name="Convention">The identifiability convention; the observation-centred one, the default and the only one built.</param>
public sealed record BassRandomIntercept(double Scale, RandomInterceptConvention Convention = RandomInterceptConvention.ObservationCentred);

/// <summary>
/// A random-intercept block on one data set (E21), integrated out. With Z the cases' level indicators over the L observed levels, D = ZᵀZ
/// their counts and H an orthonormal basis of {u : nᵀu = 0}, u = Hγ with γ ~ N(0, σ²s²I), so y has covariance σ²V,
/// V = I + s²(ZH)(ZH)ᵀ. Since nᵀH = 0, V1 = 1: the constant column is V's eigenvector with eigenvalue 1, so V^(−½) leaves the intercept's
/// column, every column's mean and ȳ as they are.
/// <list type="bullet">
/// <item><b>The marginal.</b> Whitening by V^(−½) turns the model into BASS's own, y_w = α + X_wβ + ε with ε ~ N(0, σ²I), at the price
/// −½ln|V| on the log density, which the working response carries as weights carry ½Σ ln wᵢ (<see cref="WorkingResponse"/>). With
/// HᵀDH = EΛEᵀ, V^(−½)x = x + Z·M·Zᵀx, M = HE·diag(((1 + s²λ)^(−½) − 1)/λ)·EᵀHᵀ: O(n + L²) a column, and ln|V| = Σ ln(1 + s²λ).
/// Every coefficient form then reads the whitened design as it reads a design.</item>
/// <item><b>The group effects.</b> Given α, β and σ², u's conditional is exact: with r = y − α − Xβ on the unwhitened columns,
/// γ ~ N((HᵀDH + I/s²)⁻¹HᵀZᵀr, σ²(HᵀDH + I/s²)⁻¹). A composed draw of (α, β, σ²) from the whitened model, u integrated, then u from
/// this conditional, is a joint draw.</item>
/// </list>
/// </summary>
public sealed class RandomInterceptBlock
{
    private readonly int[] _group;          // each case's observed-level index
    private readonly int[] _observedOf;     // each declared level's observed index, −1 for a level no case takes
    private readonly double[] _counts;      // n_g over the observed levels
    private readonly double[,] _whiten;     // M, L × L
    private readonly double[,] _mean;       // HE·diag(1/(λ + 1/s²))·EᵀHᵀ, L × L
    private readonly double[,] _spread;     // HE·diag((λ + 1/s²)^(−½)), L × (L − 1)

    private RandomInterceptBlock(BassRandomIntercept spec, int[] group, int[] observedOf, double[] counts, double[,] whiten, double[,] mean,
                                 double[,] spread, double logDeterminant)
    {
        Spec = spec;
        _group = group;
        _observedOf = observedOf;
        _counts = counts;
        _whiten = whiten;
        _mean = mean;
        _spread = spread;
        LogDeterminant = logDeterminant;
    }

    /// <summary>The block <paramref name="spec"/> declares on <paramref name="grouping"/>'s cases.</summary>
    public static RandomInterceptBlock Of(CategoricalPredictor grouping, BassRandomIntercept spec)
    {
        ArgumentNullException.ThrowIfNull(grouping);
        ArgumentNullException.ThrowIfNull(spec);
        if (!(spec.Scale > 0.0) || !double.IsFinite(spec.Scale)) throw new ArgumentOutOfRangeException(nameof(spec), "The scale is positive and finite.");
        if (spec.Convention != RandomInterceptConvention.ObservationCentred) throw new ArgumentOutOfRangeException(nameof(spec), "The convention is not built.");
        ReadOnlySpan<int> codes = grouping.Codes.Span;
        var observedOf = new int[grouping.Levels.Count];
        Array.Fill(observedOf, -1);
        int levels = 0;
        for (int l = 0; l < observedOf.Length; l++)
            if (grouping.Observed[l]) observedOf[l] = levels++;
        var group = new int[codes.Length];
        var counts = new double[levels];
        for (int i = 0; i < codes.Length; i++)
        {
            group[i] = observedOf[codes[i]];
            counts[group[i]]++;
        }

        // H: an orthonormal basis of n's complement, by Gram–Schmidt over the unit vectors with n̂ removed, twice for its digits.
        double norm = 0.0;
        foreach (double c in counts) norm += c * c;
        norm = Math.Sqrt(norm);
        var basis = new List<double[]>(levels - 1);
        for (int k = 0; k < levels && basis.Count < levels - 1; k++)
        {
            var v = new double[levels];
            v[k] = 1.0;
            for (int pass = 0; pass < 2; pass++)
            {
                Remove(v, counts, norm);
                foreach (double[] b in basis) Remove(v, b, 1.0);
            }
            double length = Math.Sqrt(Dot(v, v));
            if (length < 1e-8) continue;
            for (int g = 0; g < levels; g++) v[g] /= length;
            basis.Add(v);
        }
        int q = basis.Count;   // L − 1

        // HᵀDH = EΛEᵀ.
        var gram = new double[q, q];
        for (int a = 0; a < q; a++)
            for (int b = 0; b <= a; b++)
            {
                double sum = 0.0;
                for (int g = 0; g < levels; g++) sum += basis[a][g] * counts[g] * basis[b][g];
                gram[a, b] = gram[b, a] = sum;
            }
        EigenResult eigen = DenseEigen.DecomposeSymmetric(gram, tol: 1e-15);
        double s2 = spec.Scale * spec.Scale, logDeterminant = 0.0;
        // HE: each eigenvector carried back to the levels.
        var he = new double[levels, q];
        for (int k = 0; k < q; k++)
            for (int g = 0; g < levels; g++)
            {
                double sum = 0.0;
                for (int a = 0; a < q; a++) sum += basis[a][g] * eigen.Eigenvectors[k][a];
                he[g, k] = sum;
            }
        var whitenScale = new double[q];
        var meanScale = new double[q];
        var spread = new double[levels, q];
        for (int k = 0; k < q; k++)
        {
            double lambda = eigen.Eigenvalues[k];
            if (!(lambda > 0.0)) throw new InvalidOperationException("The groups' Gram on n's complement is not positive definite.");
            logDeterminant += Math.Log(1.0 + s2 * lambda);
            whitenScale[k] = (1.0 / Math.Sqrt(1.0 + s2 * lambda) - 1.0) / lambda;
            double precision = lambda + 1.0 / s2;
            meanScale[k] = 1.0 / precision;
            for (int g = 0; g < levels; g++) spread[g, k] = he[g, k] / Math.Sqrt(precision);
        }
        return new RandomInterceptBlock(spec, group, observedOf, counts, Outer(he, whitenScale), Outer(he, meanScale), spread, logDeterminant);
    }

    private static void Remove(double[] v, double[] direction, double length)
    {
        double along = Dot(v, direction) / (length * length);
        for (int g = 0; g < v.Length; g++) v[g] -= along * direction[g];
    }

    private static double Dot(double[] a, double[] b)
    {
        double sum = 0.0;
        for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
        return sum;
    }

    // A·diag(d)·Aᵀ for A L × q.
    private static double[,] Outer(double[,] a, double[] d)
    {
        int l = a.GetLength(0), q = a.GetLength(1);
        var m = new double[l, l];
        for (int g = 0; g < l; g++)
            for (int h = 0; h <= g; h++)
            {
                double sum = 0.0;
                for (int k = 0; k < q; k++) sum += a[g, k] * d[k] * a[h, k];
                m[g, h] = m[h, g] = sum;
            }
        return m;
    }

    /// <summary>The declaration.</summary>
    public BassRandomIntercept Spec { get; }

    /// <summary>L, the observed levels.</summary>
    public int Levels => _counts.Length;

    /// <summary>Each observed level's case count, n_g.</summary>
    public ReadOnlySpan<double> Counts => _counts;

    /// <summary>ln|V|, the block's price on the log density: −½ln|V| at b = 1.</summary>
    public double LogDeterminant { get; }

    /// <summary>Case <paramref name="i"/>'s observed level.</summary>
    public int GroupOf(int i) => _group[i];

    /// <summary>The observed level of declared level <paramref name="code"/>; −1 for a level no training case takes.</summary>
    public int ObservedOf(int code) => (uint)code < (uint)_observedOf.Length ? _observedOf[code] : -1;

    /// <summary>Zᵀx: <paramref name="x"/>'s sums over each observed level's cases.</summary>
    public double[] GroupSums(ReadOnlySpan<double> x)
    {
        if (x.Length != _group.Length) throw new ArgumentException($"The column holds {x.Length} values for {_group.Length} cases.", nameof(x));
        var sums = new double[_counts.Length];
        for (int i = 0; i < x.Length; i++) sums[_group[i]] += x[i];
        return sums;
    }

    /// <summary>V^(−½)x, in place: x + Z·M·Zᵀx.</summary>
    public void Whiten(Span<double> x)
    {
        double[] sums = GroupSums(x);
        int l = _counts.Length;
        var shift = new double[l];
        for (int g = 0; g < l; g++)
        {
            double sum = 0.0;
            for (int h = 0; h < l; h++) sum += _whiten[g, h] * sums[h];
            shift[g] = sum;
        }
        for (int i = 0; i < x.Length; i++) x[i] += shift[_group[i]];
    }

    /// <summary>
    /// The working response of <paramref name="working"/> whitened: V^(−½)y, as an offset y − V^(−½)y, its weighting carrying −ln|V|. Unit
    /// weights and no latents only: a weighted or latent working response would need V^(−½) under W, which is not built.
    /// </summary>
    public WorkingResponse Whitened(WorkingResponse working)
    {
        ArgumentNullException.ThrowIfNull(working);
        if (!working.IsUnweighted || working.Latents is not null || !working.Offset.IsEmpty || working.Response is FunctionalResponse)
            throw new ArgumentException("Random intercepts whiten a scalar response under unit weights, with no offset or latents.", nameof(working));
        ReadOnlySpan<double> y = working.Response.Values;
        double[] white = y.ToArray();
        Whiten(white);
        var offset = new double[y.Length];
        for (int i = 0; i < y.Length; i++) offset[i] = y[i] - white[i];
        return new WorkingResponse(working.Response, offset, reference: working.Reference, logJacobian: -LogDeterminant);
    }

    /// <summary>
    /// The group effects u over the observed levels from their conditional given α, β and σ² (<paramref name="sigma2"/>), with
    /// <paramref name="residualSums"/> = Zᵀ(y − α − Xβ) on the unwhitened columns: u = H·γ, γ's mean and spread as the block holds them.
    /// Σ_g n_g u_g = 0.
    /// </summary>
    public double[] DrawEffects(ReadOnlySpan<double> residualSums, double sigma2, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int l = _counts.Length, q = _spread.GetLength(1);
        if (residualSums.Length != l) throw new ArgumentException($"The sums hold {residualSums.Length} values for {l} levels.", nameof(residualSums));
        if (!(sigma2 > 0.0) || !double.IsFinite(sigma2)) throw new ArgumentOutOfRangeException(nameof(sigma2));
        var z = new double[q];
        for (int k = 0; k < q; k++) z[k] = Normal.SampleStandard(rng);
        double sigma = Math.Sqrt(sigma2);
        var u = new double[l];
        for (int g = 0; g < l; g++)
        {
            double mean = 0.0, noise = 0.0;
            for (int h = 0; h < l; h++) mean += _mean[g, h] * residualSums[h];
            for (int k = 0; k < q; k++) noise += _spread[g, k] * z[k];
            u[g] = mean + sigma * noise;
        }
        return u;
    }

    /// <summary>
    /// Zᵀr, r = y − α − Σⱼ βⱼBⱼ, on <paramref name="model"/>'s observations and the unwhitened case columns of <paramref name="terms"/>, with
    /// <paramref name="coefficients"/> = [α, β₁, …] in their order: the sums <see cref="DrawEffects"/> and <see cref="EffectMean"/> read.
    /// </summary>
    internal double[] ResidualSums(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients)
    {
        if (coefficients.Length != terms.Count + 1)
            throw new ArgumentException($"Give the intercept and one coefficient per term ({terms.Count + 1}).", nameof(coefficients));
        double[] sums = GroupSums(model.Data.Response.Values);
        for (int g = 0; g < sums.Length; g++) sums[g] -= coefficients[0] * _counts[g];
        var column = new double[_group.Length];
        for (int k = 0; k < terms.Count; k++)
        {
            model.Basis.CaseColumn(terms[k], model.Inputs.Cases, column);
            double[] termSums = GroupSums(column);
            for (int g = 0; g < sums.Length; g++) sums[g] -= coefficients[k + 1] * termSums[g];
        }
        return sums;
    }

    /// <summary>The group effects' conditional mean given Zᵀr = <paramref name="residualSums"/>: <see cref="DrawEffects"/>'s centre.</summary>
    public double[] EffectMean(ReadOnlySpan<double> residualSums)
    {
        int l = _counts.Length;
        var u = new double[l];
        for (int g = 0; g < l; g++)
            for (int h = 0; h < l; h++) u[g] += _mean[g, h] * residualSums[h];
        return u;
    }
}
