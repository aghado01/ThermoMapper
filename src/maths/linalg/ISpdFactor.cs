using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// A factorization of a symmetric positive-definite matrix A that a consumer can solve with and draw from
    /// without knowing its storage: <see cref="BandCholesky"/> (in LDLᵀ mode for draws) or <see cref="DenseCholesky"/>.
    /// The seam the collapsed g-prior block uses, so BARS's banded Gram and BASS's dense one share one set of
    /// closed forms. A consumer reads <see cref="Status"/> before trusting the factor; what to admit is its decision.
    /// </summary>
    public interface ISpdFactor
    {
        /// <summary>Order n of the factored matrix.</summary>
        int Dimension { get; }

        /// <summary>What the factorization found about its input.</summary>
        FactorStatus Status { get; }

        /// <summary>ln|A|.</summary>
        double LogDet { get; }

        /// <summary>Solves <c>A·x = b</c>.</summary>
        double[] Solve(ReadOnlySpan<double> b);

        /// <summary>In place: <paramref name="x"/> enters as b and leaves as <c>A⁻¹b</c>.</summary>
        void SolveInPlace(Span<double> x);

        /// <summary>Maps standard normals <paramref name="z"/> to a draw from <c>N(0, A⁻¹)</c>.</summary>
        void SampleInnovation(ReadOnlySpan<double> z, Span<double> result);
    }

    /// <summary>Quadratic forms through an <see cref="ISpdFactor"/>.</summary>
    public static class SpdFactorExtensions
    {
        /// <summary>
        /// <c>scale · z(t)ᵀA⁻¹z(t)</c> for each row z(t) of <paramref name="rows"/>, each by one solve: the pointwise
        /// variance shape of a readout whose coefficient covariance is <c>scale·A⁻¹</c>.
        /// </summary>
        public static double[] InverseQuadraticForms(this ISpdFactor factor, double[,] rows, double scale)
        {
            ArgumentNullException.ThrowIfNull(factor);
            ArgumentNullException.ThrowIfNull(rows);
            int g = rows.GetLength(0), n = rows.GetLength(1);
            if (n != factor.Dimension)
                throw new ArgumentException("Each row must have the factored matrix's order.", nameof(rows));
            var v = new double[g];
            var row = new double[n];
            for (int t = 0; t < g; t++)
            {
                for (int j = 0; j < n; j++) row[j] = rows[t, j];
                double[] sol = factor.Solve(row);
                double q = 0.0;
                for (int j = 0; j < n; j++) q += row[j] * sol[j];
                v[t] = scale * q;
            }
            return v;
        }
    }
}
