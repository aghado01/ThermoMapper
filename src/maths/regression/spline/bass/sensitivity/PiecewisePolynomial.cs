using System;
using System.Collections.Generic;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// A piecewise polynomial in one variable (O-X2): breakpoints b₀ &lt; … &lt; b_{m−1} cut the line into (−∞, b₀],
/// (b₀, b₁], …, (b_{m−1}, ∞), and each piece expands its polynomial about whatever centre it is asked for. Its expectation
/// under a measure reduces piece by piece to the measure's interval moments (<see cref="IIntervalMoments"/>) about the
/// centre of the piece clipped to the support, so no power is formed far from where the measure lives. BASS's hinge
/// products are one family (<see cref="HingeProduct"/>); another spline's pieces take the same shape.
/// </summary>
public sealed class PiecewisePolynomial
{
    /// <summary>
    /// Writes piece <paramref name="piece"/>'s coefficients about <paramref name="centre"/>, cₖ of (x − centre)ᵏ, into
    /// <paramref name="coefficients"/> (zeroed, one per power up to the degree), and returns whether the piece is nonzero.
    /// </summary>
    public delegate bool PieceExpansion(int piece, double centre, Span<double> coefficients);

    private readonly double[] _breaks;
    private readonly PieceExpansion _expand;

    /// <summary>The polynomial with ascending <paramref name="breakpoints"/>, pieces of degree at most <paramref name="degree"/>.</summary>
    public PiecewisePolynomial(IReadOnlyList<double> breakpoints, int degree, PieceExpansion expansion)
    {
        ArgumentNullException.ThrowIfNull(breakpoints);
        ArgumentNullException.ThrowIfNull(expansion);
        if (degree < 0) throw new ArgumentOutOfRangeException(nameof(degree), "A degree is not negative.");
        _breaks = [.. breakpoints];
        for (int i = 0; i < _breaks.Length; i++)
        {
            if (!double.IsFinite(_breaks[i])) throw new ArgumentException("Breakpoints are finite.", nameof(breakpoints));
            if (i > 0 && !(_breaks[i] > _breaks[i - 1])) throw new ArgumentException("Breakpoints ascend strictly.", nameof(breakpoints));
        }
        Degree = degree;
        _expand = expansion;
    }

    /// <summary>The pieces' largest degree.</summary>
    public int Degree { get; }

    /// <summary>The breakpoints, ascending.</summary>
    public ReadOnlySpan<double> Breakpoints => _breaks;

    /// <summary>The pieces, one more than the breakpoints.</summary>
    public int Pieces => _breaks.Length + 1;

    /// <summary>
    /// ∫ p dμ for <paramref name="measure"/>, whose support is [<paramref name="lower"/>, <paramref name="upper"/>], with
    /// a bound on its error from the measure's own estimates (0 when its moments are exact). Each piece passes its own
    /// bounds, so the measure clips them and counts an atom on a breakpoint once.
    /// </summary>
    public (double Value, double Error) Expectation(IIntervalMoments measure, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(measure);
        if (!double.IsFinite(lower) || !double.IsFinite(upper) || lower > upper)
            throw new ArgumentException($"A support needs finite lower ≤ upper, got [{lower}, {upper}].");
        Span<double> coefficients = stackalloc double[Degree + 1];
        Span<double> moments = stackalloc double[Degree + 1];
        double value = 0.0, error = 0.0;
        for (int p = 0; p < Pieces; p++)
        {
            double lo = p == 0 ? double.NegativeInfinity : _breaks[p - 1];
            double hi = p == _breaks.Length ? double.PositiveInfinity : _breaks[p];
            double a = Math.Max(lo, lower), b = Math.Min(hi, upper);
            if (b < a) continue;   // a point piece is kept: it may hold an atom
            double centre = 0.5 * (a + b);
            coefficients.Clear();
            if (!_expand(p, centre, coefficients)) continue;
            double e = measure.IntervalMoments(lo, hi, centre, 0.0, moments);
            double sum = 0.0, size = 0.0;
            for (int k = 0; k <= Degree; k++)
            {
                sum += coefficients[k] * moments[k];
                size += Math.Abs(coefficients[k]);
            }
            value += sum;
            error += size * e;
        }
        return (value, error);
    }

    /// <summary>
    /// Πᵢ ((sᵢ(x − tᵢ))₊^q − mᵢ), the product of hinges of degree <paramref name="degree"/> with signs and knots
    /// <paramref name="hinges"/>, each less its offset in <paramref name="offsets"/> (none: all 0). Breakpoints at the
    /// distinct knots; a hinge is its polynomial on the pieces it faces and 0 on the others. With offsets the hinges' means,
    /// it is the integrand of their covariance, which then forms no difference of large moments.
    /// </summary>
    public static PiecewisePolynomial HingeProduct(int degree, IReadOnlyList<(int Sign, double Knot)> hinges, IReadOnlyList<double>? offsets = null)
    {
        ArgumentNullException.ThrowIfNull(hinges);
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), "A hinge's degree is a positive integer.");
        if (offsets is not null && offsets.Count != hinges.Count) throw new ArgumentException("One offset per hinge.", nameof(offsets));
        var signs = new int[hinges.Count];
        var knots = new double[hinges.Count];
        var shift = new double[hinges.Count];
        for (int i = 0; i < hinges.Count; i++)
        {
            if (hinges[i].Sign is not (1 or -1)) throw new ArgumentOutOfRangeException(nameof(hinges), "A hinge's sign is +1 or −1.");
            signs[i] = hinges[i].Sign;
            knots[i] = hinges[i].Knot;
            shift[i] = offsets?[i] ?? 0.0;
        }
        var breaks = new SortedSet<double>(knots);
        double[] b = [.. breaks];
        bool Expand(int piece, double centre, Span<double> coefficients)
        {
            double lo = piece == 0 ? double.NegativeInfinity : b[piece - 1];
            double hi = piece == b.Length ? double.PositiveInfinity : b[piece];
            coefficients[0] = 1.0;
            int top = 0;   // the product's degree so far
            Span<double> factor = stackalloc double[degree + 1];
            for (int i = 0; i < signs.Length; i++)
            {
                bool faces = signs[i] > 0 ? lo >= knots[i] : hi <= knots[i];
                // The hinge about the centre: (s(c − t) + s·y)^q, less its offset; or −m off its side.
                factor.Clear();
                int width;
                if (faces)
                {
                    double alpha = signs[i] * (centre - knots[i]), beta = signs[i];
                    factor[0] = 1.0;
                    for (int r = 0; r < degree; r++)
                    {
                        for (int k = r + 1; k > 0; k--) factor[k] = factor[k] * alpha + factor[k - 1] * beta;
                        factor[0] *= alpha;
                    }
                    factor[0] -= shift[i];
                    width = degree;
                }
                else
                {
                    if (shift[i] == 0.0) return false;
                    factor[0] = -shift[i];
                    width = 0;
                }
                for (int k = top + width; k >= 0; k--)
                {
                    double sum = 0.0;
                    for (int r = Math.Max(0, k - top); r <= Math.Min(k, width); r++) sum += coefficients[k - r] * factor[r];
                    coefficients[k] = sum;
                }
                top += width;
            }
            return true;
        }
        return new PiecewisePolynomial(b, degree * hinges.Count, Expand);
    }
}
