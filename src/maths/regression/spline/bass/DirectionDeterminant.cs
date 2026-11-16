using System;
using System.Collections.Generic;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bass;

/// <summary>The prior on a directional fit's directions given their number (E36's decision 2).</summary>
public enum DirectionPrior
{
    /// <summary>The paper's: each direction independent, a uniform, 𝒥 uniform, θ_𝒥 uniform on its coordinate subsphere. The default.</summary>
    Independent,

    /// <summary>
    /// The determinant: the independent law reweighted by det(ΘᵀΘ)/c_M, the squared volume the directions span over its mean
    /// (<see cref="DirectionDeterminant"/>), a labelled option for many inputs and few sparse ridges (gcbeck's repulsion, with its normalizer
    /// stated). It has no mass past M = p, nor on directions with a linear dependence: a sign pair, a repeat, r + 1 directions in an r-plane.
    /// </summary>
    Determinant,
}

/// <summary>
/// The determinant direction prior's terms (E36's slice 4): π(Θ | M) = Π_m π(θ_m) · det(ΘᵀΘ)/c_M over the structure terms' unit directions
/// θ_m ∈ ℝ^p, with c_M = E[det(ΘᵀΘ)] under the independent law. That law's directions are exchangeable over the coordinates and symmetric in
/// each sign, so E[θθᵀ] = I/p whatever a and 𝒥's laws; by Cauchy–Binet and the columns' independence E[det ΘᵀΘ] = Σ_{|S| = M} M!·det(I/p)_S,
/// so c_M = p!/((p − M)!·p^M) = Π_{j &lt; M}(p − j)/p for sparse directions as for dense ones (the Grassmann discussion derived the dense case;
/// gcbeck's code leaves c_M out, which tilts the count law by it). ln det(ΘᵀΘ) is read through the Gram's Cholesky factor; a pivot at most
/// 10⁻¹³ times its diagonal entry reads as zero, a dependence the rounding cannot resolve, which removes the mass where a direction's sin² to
/// the earlier ones' span is below 10⁻¹³.
/// </summary>
public static class DirectionDeterminant
{
    private const double PivotFloor = 1e-13;

    /// <summary>ln c_M = ln p! − ln (p − M)! − M ln p for <paramref name="count"/> = M directions in <paramref name="variables"/> = p
    /// dimensions, summed as Σ_(j &lt; M) ln(1 − j/p), exact at M ≤ 1; −∞ past M = p, where every M directions are dependent.</summary>
    public static double LogNormalizer(int variables, int count)
    {
        if (variables < 1) throw new ArgumentOutOfRangeException(nameof(variables), "p is at least 1.");
        if (count < 0) throw new ArgumentOutOfRangeException(nameof(count), "M is not negative.");
        if (count > variables) return double.NegativeInfinity;
        double sum = 0.0;
        for (int j = 1; j < count; j++) sum += SpecialFunctions.Log1p(-(double)j / variables);
        return sum;
    }

    /// <summary>ln det(ΘᵀΘ) for <paramref name="directions"/>' Gram, 0 for none; −∞ where they are linearly dependent.</summary>
    public static double LogDeterminant(IReadOnlyList<ProjectionDirection> directions)
    {
        ArgumentNullException.ThrowIfNull(directions);
        int m = directions.Count;
        if (m == 0) return 0.0;
        var l = new double[m, m];
        double logDet = 0.0;
        for (int i = 0; i < m; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = Dot(directions[i], directions[j]);
                for (int k = 0; k < j; k++) sum -= l[i, k] * l[j, k];
                if (i == j)
                {
                    if (!(sum > PivotFloor * Dot(directions[i], directions[i]))) return double.NegativeInfinity;
                    l[i, i] = Math.Sqrt(sum);
                    logDet += Math.Log(sum);
                }
                else l[i, j] = sum / l[j, j];
            }
        return logDet;
    }

    /// <summary>ln det(ΘᵀΘ) − ln c_(M_θ) over the directions of <paramref name="terms"/>' M_θ direction-bearing structure terms, in
    /// <paramref name="variables"/> = p dimensions: the determinant prior's coupling beside the terms' own priors; −∞ past M_θ = p or on a
    /// dependence. A term without a direction (a categorical part alone, E36's slice 5) takes no part in it, and c is read at M_θ, so the
    /// prior on which terms bear a direction stays the term measure's and only their directions are reweighted.</summary>
    public static double LogCoupling(IReadOnlyList<BassTerm> terms, int variables)
    {
        ArgumentNullException.ThrowIfNull(terms);
        var directions = new List<ProjectionDirection>(terms.Count);
        foreach (BassTerm term in terms)
            if (term.Directional is { } hinge) directions.Add(hinge.Direction);
        if (directions.Count > variables) return double.NegativeInfinity;
        double logDet = LogDeterminant(directions);
        return double.IsNegativeInfinity(logDet) ? logDet : logDet - LogNormalizer(variables, directions.Count);
    }

    // θᵀθ′ over the two directions' active variables, each ascending.
    private static double Dot(ProjectionDirection a, ProjectionDirection b)
    {
        IReadOnlyList<int> va = a.Variables, vb = b.Variables;
        IReadOnlyList<double> wa = a.Weights, wb = b.Weights;
        double sum = 0.0;
        for (int i = 0, j = 0; i < va.Count && j < vb.Count;)
        {
            if (va[i] == vb[j]) sum += wa[i++] * wb[j++];
            else if (va[i] < vb[j]) i++;
            else j++;
        }
        return sum;
    }
}
