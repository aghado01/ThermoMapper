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
    /// <item><b>Tensor designs</b> (N-X1), whose rows are themselves face-splitting products: row i of X = B₁ ∗ … ∗ B_d is
    /// the row-wise Kronecker product b¹ᵢ ⊗ … ⊗ bᵈᵢ of the axes' basis rows, as in a free-knot tensor-product spline surface.
    /// Scattered samples accumulate G = Σᵢ xᵢxᵢᵀ over each row's nonzeros, O(n·(p+1)^(2d)) for B-splines of degree p, and a
    /// full grid takes the Kronecker shortcut G = ⊗ₐ(BₐᵀBₐ). The Gram's columns are lexicographic with the first axis
    /// slowest, so it is block-banded, with half-bandwidth <see cref="FaceSplittingHalfBandwidth"/>.</item>
    /// </list>
    /// Every inner product is a <see cref="PairwiseSum"/>. Factors are passed as columns: a[j] of length n, h[j] of
    /// length m.
    /// <para><b>The output-side inner product (F-X6).</b> The direct contractions are written with Q, the declared inner
    /// product on the output index (<see cref="IInnerProduct"/>): G = (AᵀA)∘(HᵀQH) and h_j = a_jᵀ·Y·Q·h_j, the weight of a
    /// residual with covariance σ²(I ⊗ R) being I ⊗ Q with Q = R⁻¹. Only the output-side factors change; the case side and
    /// the separable algebra stay. A null product is Q = I, which takes the unweighted arithmetic above exactly. The centred
    /// forms serve Q = I: under Q the intercept's orthogonal complement is the Q-centred one, and the three-term form
    /// carries over with m replaced by 1ᵀQ1, the means by Q-weighted means (1ᵀQh/1ᵀQ1) and the factors' inner products by
    /// hᵀQh, which is the correlated-residual extension's to write.</para>
    /// </summary>
    public static class KhatriRao
    {
        /// <summary>G = XᵀW X = (AᵀA)∘(HᵀQH), with W = I ⊗ Q; XᵀX = (AᵀA)∘(HᵀH) for a null <paramref name="output"/>.</summary>
        public static double[,] Gram(double[][] a, double[][] h, IInnerProduct? output = null)
        {
            int count = Validate(a, h, out _, out int m);
            IInnerProduct q = Product(output, m);
            var g = new double[count, count];
            for (int i = 0; i < count; i++)
                for (int j = 0; j <= i; j++)
                    g[i, j] = g[j, i] = PairwiseSum.Dot(a[i], a[j]) * q.Dot(h[i], h[j]);
            return g;
        }

        /// <summary>A new column a_new ⊗ h_new against the existing ones, in order (a birth's border), and its squared
        /// norm in <paramref name="diagonal"/>.</summary>
        public static double[] Border(double[][] a, double[][] h, double[] aNew, double[] hNew, out double diagonal,
                                      IInnerProduct? output = null)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireLength(aNew, n, nameof(aNew));
            RequireLength(hNew, m, nameof(hNew));
            IInnerProduct q = Product(output, m);
            var u = new double[count];
            for (int k = 0; k < count; k++) u[k] = PairwiseSum.Dot(a[k], aNew) * q.Dot(h[k], hNew);
            diagonal = PairwiseSum.Dot(aNew, aNew) * q.Dot(hNew, hNew);
            return u;
        }

        /// <summary>XᵀW vec(Y): h_j = a_jᵀ·Y·Q·h_j, with Y n × m; Xᵀvec(Y) for a null <paramref name="output"/>.</summary>
        public static double[] CrossProducts(double[][] a, double[][] h, double[,] y, IInnerProduct? output = null)
        {
            int count = Validate(a, h, out int n, out int m);
            RequireResponse(y, n, m);
            IInnerProduct q = Product(output, m);
            var result = new double[count];
            var yh = new double[n];
            for (int j = 0; j < count; j++)
            {
                for (int i = 0; i < n; i++) yh[i] = q.Dot(Row(y, i, m), h[j]);
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
            int[] caseOf = cases.ToArray(), locationOf = locations.ToArray();
            return AccumulateGram(count, caseOf.Length, (r, z, nonzero) =>
            {
                int k = 0;
                for (int j = 0; j < count; j++)
                {
                    z[j] = a[j][caseOf[r]] * h[j][locationOf[r]];
                    if (z[j] != 0.0) nonzero[k++] = j;
                }
                return k;
            });
        }

        /// <summary>
        /// G = XᵀX for the tensor design X = B₁ ∗ … ∗ B_d (N-X1), the face-splitting product of the axes' bases
        /// <paramref name="axes"/>, each n × νₐ for the same n samples: row i of X is the row-wise Kronecker product
        /// b¹ᵢ ⊗ … ⊗ bᵈᵢ, columns lexicographic with axis 0 slowest. Each row's zeros are skipped, so a B-spline basis, with
        /// p + 1 nonzeros a row on each axis, costs O(n·(p+1)^(2d)). Dense, for the designs whose Πνₐ a dense Gram holds.
        /// </summary>
        public static double[,] FaceSplittingGram(IReadOnlyList<double[,]> axes)
        {
            int rows = FaceSplittingShape(axes, out int[] widths, out int count);
            var row = new FaceSplittingRow(axes, widths, count);
            return AccumulateGram(count, rows, row.Fill);
        }

        /// <summary>Xᵀy for the tensor design of <see cref="FaceSplittingGram"/> and one response per sample.</summary>
        public static double[] FaceSplittingCrossProducts(IReadOnlyList<double[,]> axes, ReadOnlySpan<double> y)
        {
            int rows = FaceSplittingShape(axes, out int[] widths, out int count);
            if (y.Length != rows) throw new ArgumentException("Expected one response per sample.", nameof(y));
            var row = new FaceSplittingRow(axes, widths, count);
            double[] response = y.ToArray();
            var levels = new List<double[]?>();
            var block = new double[count];
            var z = new double[count];
            var nonzero = new int[count];
            int inBlock = 0;
            for (int r = 0; r < rows; r++)
            {
                int k = row.Fill(r, z, nonzero);
                for (int p = 0; p < k; p++) block[nonzero[p]] += z[nonzero[p]] * response[r];
                if (++inBlock == PairwiseSum.Block)
                {
                    Carry(levels, block);
                    block = new double[count];
                    inBlock = 0;
                }
            }
            if (inBlock > 0) Carry(levels, block);
            var total = new double[count];
            foreach (double[]? level in levels)
            {
                if (level is null) continue;
                for (int t = 0; t < count; t++) total[t] += level[t];
            }
            return total;
        }

        /// <summary>
        /// The half-bandwidth of the Gram of the tensor design of <see cref="FaceSplittingGram"/>: the widest span between a
        /// row's first and last nonzero column, which bounds |i − j| over the Gram's nonzero entries and meets it. For linear
        /// B-splines on two axes it is p·ν₂ + p at p = 1; a caller that declares it chooses the
        /// factor of the Gram by it.
        /// </summary>
        public static int FaceSplittingHalfBandwidth(IReadOnlyList<double[,]> axes)
        {
            int rows = FaceSplittingShape(axes, out int[] widths, out int count);
            var stride = new long[widths.Length];
            long s = 1;
            for (int a = widths.Length - 1; a >= 0; a--)
            {
                stride[a] = s;
                s *= widths[a];
            }
            long band = 0;
            for (int r = 0; r < rows; r++)
            {
                long first = 0, last = 0;
                bool empty = false;
                for (int a = 0; a < widths.Length && !empty; a++)
                {
                    int lo = -1, hi = -1;
                    for (int j = 0; j < widths[a]; j++)
                        if (axes[a][r, j] != 0.0)
                        {
                            if (lo < 0) lo = j;
                            hi = j;
                        }
                    if (lo < 0) empty = true;
                    else
                    {
                        first += lo * stride[a];
                        last += hi * stride[a];
                    }
                }
                if (!empty) band = Math.Max(band, last - first);
            }
            return checked((int)band);
        }

        /// <summary>
        /// G = ⊗ₐ(BₐᵀBₐ) for a tensor design on a full grid (N-X1's Cartesian shortcut): the samples are every combination of
        /// the axes' own samples, so X = B₁ ⊗ … ⊗ B_d and its Gram is the Kronecker product of the axes' Grams. The axes' row
        /// counts are their own grids' sizes. Columns lexicographic with axis 0 slowest, as <see cref="FaceSplittingGram"/>.
        /// </summary>
        public static double[,] KroneckerGram(IReadOnlyList<double[,]> axes)
        {
            ArgumentNullException.ThrowIfNull(axes);
            if (axes.Count == 0) throw new ArgumentException("A tensor design needs an axis.", nameof(axes));
            double[,]? gram = null;
            foreach (double[,] axis in axes)
            {
                ArgumentNullException.ThrowIfNull(axis);
                int rows = axis.GetLength(0), width = axis.GetLength(1);
                if (rows == 0 || width == 0) throw new ArgumentException("Every axis needs samples and basis functions.", nameof(axes));
                var columns = new double[width][];
                for (int j = 0; j < width; j++)
                {
                    columns[j] = new double[rows];
                    for (int i = 0; i < rows; i++) columns[j][i] = axis[i, j];
                }
                var g = new double[width, width];
                for (int i = 0; i < width; i++)
                    for (int j = 0; j <= i; j++) g[i, j] = g[j, i] = PairwiseSum.Dot(columns[i], columns[j]);
                if (gram is null)
                {
                    gram = g;
                    continue;
                }
                int outer = gram.GetLength(0);
                if ((long)outer * width > MaxDenseOrder) throw new ArgumentException("The Gram is too large to hold densely.", nameof(axes));
                var product = new double[outer * width, outer * width];
                for (int i = 0; i < outer; i++)
                    for (int j = 0; j < outer; j++)
                        for (int k = 0; k < width; k++)
                            for (int l = 0; l < width; l++) product[i * width + k, j * width + l] = gram[i, j] * g[k, l];
                gram = product;
            }
            return gram!;
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

        private const int MaxDenseOrder = 46340;

        // One row's nonzero columns, ascending, and their values: fills z at those columns and returns how many.
        private delegate int RowFiller(int row, double[] z, int[] nonzero);

        // G = Σ_r z_r z_rᵀ over `rows` rows of `count` columns. Blocks of PairwiseSum.Block rows are summed directly, and the
        // blocks pairwise, by a binary cascade of partial sums (one per level), so memory stays O(count² log rows).
        private static double[,] AccumulateGram(int count, int rows, RowFiller fill)
        {
            if (count > MaxDenseOrder) throw new ArgumentException($"{count} columns are too many for a dense Gram.");
            int packed = count * (count + 1) / 2;
            var levels = new List<double[]?>();
            var block = new double[packed];
            var z = new double[count];
            var nonzero = new int[count];
            int inBlock = 0;
            for (int r = 0; r < rows; r++)
            {
                int k = fill(r, z, nonzero);
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

        // The axes' shared sample count, widths and the product of the widths.
        private static int FaceSplittingShape(IReadOnlyList<double[,]> axes, out int[] widths, out int count)
        {
            ArgumentNullException.ThrowIfNull(axes);
            if (axes.Count == 0) throw new ArgumentException("A tensor design needs an axis.", nameof(axes));
            int rows = -1;
            widths = new int[axes.Count];
            long product = 1;
            for (int a = 0; a < axes.Count; a++)
            {
                double[,] axis = axes[a] ?? throw new ArgumentException($"Axis {a} is null.", nameof(axes));
                if (rows < 0) rows = axis.GetLength(0);
                if (axis.GetLength(0) != rows) throw new ArgumentException("Every axis needs one row per sample.", nameof(axes));
                widths[a] = axis.GetLength(1);
                if (widths[a] == 0) throw new ArgumentException($"Axis {a} has no basis function.", nameof(axes));
                product *= widths[a];
                if (product > MaxDenseOrder) throw new ArgumentException("The Gram is too large to hold densely.", nameof(axes));
            }
            count = (int)product;
            return rows;
        }

        // Row i of a face-splitting product, expanded axis by axis from the axes' nonzeros, which come out ascending.
        private sealed class FaceSplittingRow
        {
            private readonly IReadOnlyList<double[,]> _axes;
            private readonly int[] _widths;
            private int[] _columns, _nextColumns;
            private double[] _values, _nextValues;

            public FaceSplittingRow(IReadOnlyList<double[,]> axes, int[] widths, int count)
            {
                _axes = axes;
                _widths = widths;
                _columns = new int[count];
                _nextColumns = new int[count];
                _values = new double[count];
                _nextValues = new double[count];
            }

            public int Fill(int row, double[] z, int[] nonzero)
            {
                int total = 1;
                _columns[0] = 0;
                _values[0] = 1.0;
                for (int a = 0; a < _axes.Count; a++)
                {
                    int next = 0, width = _widths[a];
                    for (int t = 0; t < total; t++)
                        for (int j = 0; j < width; j++)
                        {
                            double v = _axes[a][row, j];
                            if (v == 0.0) continue;
                            _nextColumns[next] = _columns[t] * width + j;
                            _nextValues[next++] = _values[t] * v;
                        }
                    (_columns, _nextColumns) = (_nextColumns, _columns);
                    (_values, _nextValues) = (_nextValues, _values);
                    total = next;
                }
                for (int t = 0; t < total; t++)
                {
                    z[_columns[t]] = _values[t];
                    nonzero[t] = _columns[t];
                }
                return total;
            }
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

        // The declared product, the identity when none; one that acts on another output length is refused.
        private static IInnerProduct Product(IInnerProduct? output, int m)
        {
            if (output is null) return IdentityInnerProduct.Instance;
            if (output.Dimension != 0 && output.Dimension != m)
                throw new ArgumentException($"The inner product acts on {output.Dimension} output locations, not {m}.", nameof(output));
            return output;
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
