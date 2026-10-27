using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace Maths.LinAlg
{
    /// <summary>
    /// Contractions with the column-wise Khatri–Rao product X = A ⊙ H, whose column j is a_j ⊗ h_j, computed through
    /// its factors without forming X: BASS's direct functional design, with A the case-side factors (n × M) and H the
    /// output-side factors (m × M), rows of X indexed (case, location) with the case outer. Xᵀ is the face-splitting
    /// product of Aᵀ and Hᵀ, the name the BASS port's numerics brief uses.
    /// <list type="bullet">
    /// <item><b>Cartesian samples</b>, every case at every location, take the Hadamard shortcut
    /// G = (AᵀA)∘(HᵀH) and h_j = a_jᵀ·Y·h_j: O(M²(n + m)) and O(M·n·m) instead of O(M²·n·m).</item>
    /// <item><b>Centred</b> forms serve the Liang intercept convention, where the intercept leaves the design and the
    /// other columns are centred over all n·m rows. With α, η the centred factors and ā, h̄ their means,
    /// G_c = (A_cᵀA_c)∘(H_cᵀH_c) + n·(āāᵀ)∘(H_cᵀH_c) + m·(h̄h̄ᵀ)∘(A_cᵀA_c), a sum of three positive semidefinite terms,
    /// so the diagonal forms without the cancellation of XᵀX − nm·x̄x̄ᵀ (the BASS port's fixture F-O7).</item>
    /// <item><b>Scattered samples</b>, given as (case, location) pairs, accumulate each row's Hadamard product
    /// a_i ∘ h_l, skipping its zeros.</item>
    /// </list>
    /// Every inner product is a <see cref="PairwiseSum"/>. Factors are passed as columns: a[j] of length n, h[j] of
    /// length m.
    /// </summary>
    public static class KhatriRao
    {
        /// <summary>G = XᵀX = (AᵀA)∘(HᵀH).</summary>
        public static double[,] Gram(double[][] a, double[][] h)
        {
            int count = Validate(a, h, out _, out _);
            var g = new double[count, count];
            for (int i = 0; i < count; i++)
                for (int j = 0; j <= i; j++)
                    g[i, j] = g[j, i] = PairwiseSum.Dot(a[i], a[j]) * PairwiseSum.Dot(h[i], h[j]);
            return g;
        }

        /// <summary>A new column a_new ⊗ h_new against the existing ones, in order (a birth's border), and its squared
        /// norm in <paramref name="diagonal"/>.</summary>
        public static double[] Border(double[][] a, double[][] h, double[] aNew, double[] hNew, out double diagonal)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireLength(aNew, n, nameof(aNew));
            RequireLength(hNew, m, nameof(hNew));
            var u = new double[count];
            for (int k = 0; k < count; k++) u[k] = PairwiseSum.Dot(a[k], aNew) * PairwiseSum.Dot(h[k], hNew);
            diagonal = PairwiseSum.Dot(aNew, aNew) * PairwiseSum.Dot(hNew, hNew);
            return u;
        }

        /// <summary>Xᵀvec(Y): h_j = a_jᵀ·Y·h_j, with Y n × m.</summary>
        public static double[] CrossProducts(double[][] a, double[][] h, double[,] y)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireResponse(y, n, m);
            var result = new double[count];
            var yh = new double[n];
            for (int j = 0; j < count; j++)
            {
                for (int i = 0; i < n; i++) yh[i] = PairwiseSum.Dot(Row(y, i, m), h[j]);
                result[j] = PairwiseSum.Dot(a[j], yh);
            }
            return result;
        }

        /// <summary>X_cᵀX_c for the columns centred over all n·m rows, by the three-term form.</summary>
        public static double[,] CentredGram(double[][] a, double[][] h)
        {
            int count = Validate(a, h, out int n, out int m);
            Centre(a, out double[][] alpha, out double[] aBar);
            Centre(h, out double[][] eta, out double[] hBar);
            var g = new double[count, count];
            for (int i = 0; i < count; i++)
                for (int j = 0; j <= i; j++)
                    g[i, j] = g[j, i] = CentredEntry(alpha[i], eta[i], aBar[i], hBar[i], alpha[j], eta[j], aBar[j], hBar[j], n, m);
            return g;
        }

        /// <summary>The centred counterpart of <see cref="Border"/>.</summary>
        public static double[] CentredBorder(double[][] a, double[][] h, double[] aNew, double[] hNew, out double diagonal)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireLength(aNew, n, nameof(aNew));
            RequireLength(hNew, m, nameof(hNew));
            Centre(a, out double[][] alpha, out double[] aBar);
            Centre(h, out double[][] eta, out double[] hBar);
            double[] alphaNew = Centred(aNew, out double aBarNew);
            double[] etaNew = Centred(hNew, out double hBarNew);
            var u = new double[count];
            for (int k = 0; k < count; k++)
                u[k] = CentredEntry(alpha[k], eta[k], aBar[k], hBar[k], alphaNew, etaNew, aBarNew, hBarNew, n, m);
            diagonal = CentredEntry(alphaNew, etaNew, aBarNew, hBarNew, alphaNew, etaNew, aBarNew, hBarNew, n, m);
            return u;
        }

        /// <summary>
        /// X_cᵀvec(Y) for the centred columns. Each centred column splits as ā_j(1 ⊗ η_j) + h̄_j(α_j ⊗ 1) + α_j ⊗ η_j, so
        /// the entry is ā_j·1ᵀYη_j + h̄_j·α_jᵀY1 + α_jᵀYη_j, read from Y less its grand mean (which the centred columns
        /// do not see) to keep the large common level out of the sums.
        /// </summary>
        public static double[] CentredCrossProducts(double[][] a, double[][] h, double[,] y)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireResponse(y, n, m);
            Centre(a, out double[][] alpha, out double[] aBar);
            Centre(h, out double[][] eta, out double[] hBar);

            var rows = new double[n][];                     // Y less its grand mean, in two passes
            for (int i = 0; i < n; i++) rows[i] = Row(y, i, m).ToArray();
            var rowMeans = new double[n];
            for (int pass = 0; pass < 2; pass++)
            {
                for (int i = 0; i < n; i++) rowMeans[i] = PairwiseSum.Of(rows[i]) / m;
                double level = PairwiseSum.Of(rowMeans) / n;
                for (int i = 0; i < n; i++) for (int l = 0; l < m; l++) rows[i][l] -= level;
            }

            var rowSums = new double[n];                    // Y_c·1
            for (int i = 0; i < n; i++) rowSums[i] = PairwiseSum.Of(rows[i]);
            var columnSums = new double[m];                 // 1ᵀY_c
            var column = new double[n];
            for (int l = 0; l < m; l++)
            {
                for (int i = 0; i < n; i++) column[i] = rows[i][l];
                columnSums[l] = PairwiseSum.Of(column);
            }

            var result = new double[count];
            var yEta = new double[n];
            for (int j = 0; j < count; j++)
            {
                for (int i = 0; i < n; i++) yEta[i] = PairwiseSum.Dot(rows[i], eta[j]);
                result[j] = aBar[j] * PairwiseSum.Dot(columnSums, eta[j])
                          + hBar[j] * PairwiseSum.Dot(alpha[j], rowSums)
                          + PairwiseSum.Dot(alpha[j], yEta);
            }
            return result;
        }

        /// <summary>
        /// G = Σ_r z_r z_rᵀ over scattered samples, z_r = a_{i_r} ∘ h_{l_r} the Hadamard product of case i_r's and
        /// location l_r's factor rows, each row's zeros skipped (ragged observations; ABMRS's scattered inputs).
        /// </summary>
        public static double[,] ScatteredGram(double[][] a, double[][] h, ReadOnlySpan<int> cases, ReadOnlySpan<int> locations)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireSamples(cases, locations, n, m);
            // Blocks of PairwiseSum.Block samples are summed directly, and the blocks pairwise, by a binary cascade of
            // partial sums (one per level), so memory stays O(M² log(samples)).
            int packed = count * (count + 1) / 2;
            var levels = new List<double[]?>();
            var block = new double[packed];
            var z = new double[count];
            var nonzero = new int[count];
            int inBlock = 0;
            for (int r = 0; r < cases.Length; r++)
            {
                int k = 0;
                for (int j = 0; j < count; j++)
                {
                    z[j] = a[j][cases[r]] * h[j][locations[r]];
                    if (z[j] != 0.0) nonzero[k++] = j;
                }
                for (int p = 0; p < k; p++)
                    for (int q = 0; q <= p; q++)
                    {
                        int i = nonzero[p], j = nonzero[q];
                        block[i * (i + 1) / 2 + j] += z[i] * z[j];
                    }
                if (++inBlock == PairwiseSum.Block)
                {
                    Carry(levels, block);
                    block = new double[packed];
                    inBlock = 0;
                }
            }
            if (inBlock > 0) Carry(levels, block);
            double[]? total = null;
            foreach (double[]? level in levels)
            {
                if (level is null) continue;
                if (total is null) total = level;
                else for (int t = 0; t < packed; t++) total[t] += level[t];
            }
            var g = new double[count, count];
            if (total is null) return g;
            for (int i = 0; i < count; i++)
                for (int j = 0; j <= i; j++) g[i, j] = g[j, i] = total[i * (i + 1) / 2 + j];
            return g;
        }

        /// <summary>Xᵀy over scattered samples: h_j = Σ_r a_j[i_r]·h_j[l_r]·y_r.</summary>
        public static double[] ScatteredCrossProducts(double[][] a, double[][] h, ReadOnlySpan<int> cases, ReadOnlySpan<int> locations, ReadOnlySpan<double> y)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireSamples(cases, locations, n, m);
            if (y.Length != cases.Length) throw new ArgumentException("Expected one response per sample.", nameof(y));
            var result = new double[count];
            var terms = new double[cases.Length];
            for (int j = 0; j < count; j++)
            {
                for (int r = 0; r < cases.Length; r++) terms[r] = a[j][cases[r]] * h[j][locations[r]] * y[r];
                result[j] = PairwiseSum.Of(terms);
            }
            return result;
        }

        // Adds a block's partial sums to the cascade as a binary counter adds one: equal-sized partial sums merge.
        private static void Carry(List<double[]?> levels, double[] sums)
        {
            for (int k = 0; ; k++)
            {
                if (k == levels.Count)
                {
                    levels.Add(sums);
                    return;
                }
                double[]? held = levels[k];
                if (held is null)
                {
                    levels[k] = sums;
                    return;
                }
                for (int t = 0; t < sums.Length; t++) sums[t] = held[t] + sums[t];
                levels[k] = null;
            }
        }

        // ⟨x_c,i, x_c,j⟩ = ⟨α_i,α_j⟩⟨η_i,η_j⟩ + n·ā_iā_j⟨η_i,η_j⟩ + m·h̄_ih̄_j⟨α_i,α_j⟩; the cross terms vanish because the
        // centred factors are orthogonal to the constant.
        private static double CentredEntry(double[] alphaI, double[] etaI, double aBarI, double hBarI,
                                           double[] alphaJ, double[] etaJ, double aBarJ, double hBarJ, int n, int m)
        {
            double caseInner = PairwiseSum.Dot(alphaI, alphaJ);
            double locationInner = PairwiseSum.Dot(etaI, etaJ);
            return caseInner * locationInner + n * aBarI * aBarJ * locationInner + m * hBarI * hBarJ * caseInner;
        }

        private static void Centre(double[][] columns, out double[][] centred, out double[] means)
        {
            centred = new double[columns.Length][];
            means = new double[columns.Length];
            for (int j = 0; j < columns.Length; j++) centred[j] = Centred(columns[j], out means[j]);
        }

        // Two-pass centring: subtract the mean, then the mean of what is left.
        private static double[] Centred(double[] x, out double mean)
        {
            int n = x.Length;
            mean = PairwiseSum.Of(x) / n;
            var c = new double[n];
            for (int i = 0; i < n; i++) c[i] = x[i] - mean;
            double residual = PairwiseSum.Of(c) / n;
            for (int i = 0; i < n; i++) c[i] -= residual;
            return c;
        }

        private static ReadOnlySpan<double> Row(double[,] y, int i, int m) =>
            MemoryMarshal.CreateReadOnlySpan(ref y[i, 0], m);

        private static int Validate(double[][] a, double[][] h, out int n, out int m)
        {
            ArgumentNullException.ThrowIfNull(a);
            ArgumentNullException.ThrowIfNull(h);
            if (a.Length != h.Length) throw new ArgumentException("A and H need one column per term each.", nameof(h));
            n = a.Length > 0 ? a[0].Length : 0;
            m = h.Length > 0 ? h[0].Length : 0;
            for (int j = 0; j < a.Length; j++)
            {
                RequireLength(a[j], n, nameof(a));
                RequireLength(h[j], m, nameof(h));
            }
            return a.Length;
        }

        private static void RequireLength(double[] x, int length, string name)
        {
            ArgumentNullException.ThrowIfNull(x, name);
            if (x.Length != length) throw new ArgumentException($"Expected length {length}, got {x.Length}.", name);
        }

        private static void RequireResponse(double[,] y, int n, int m)
        {
            ArgumentNullException.ThrowIfNull(y);
            if (y.GetLength(0) != n || y.GetLength(1) != m)
                throw new ArgumentException($"Expected the response as {n} cases × {m} locations.", nameof(y));
        }

        private static void RequireSamples(ReadOnlySpan<int> cases, ReadOnlySpan<int> locations, int n, int m)
        {
            if (cases.Length != locations.Length) throw new ArgumentException("Expected one location per sampled case.", nameof(locations));
            for (int r = 0; r < cases.Length; r++)
                if ((uint)cases[r] >= (uint)n || (uint)locations[r] >= (uint)m)
                    throw new ArgumentOutOfRangeException(nameof(cases), $"Sample {r} is ({cases[r]}, {locations[r]}), outside {n} × {m}.");
        }
    }
}
