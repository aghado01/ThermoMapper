using System;
using Maths.Distributions;
using Maths.LinAlg;

namespace Maths.Geometry.Inference
{
    /// <summary>
    /// Maximum pairwise Bayes factor (mxPBF) tests for high-dimensional two-sample problems,
    /// after Lee, You &amp; Lin (arXiv 2112.02580). Each test scans every coordinate (or
    /// coordinate pair) and reports the log of the single largest pairwise Bayes factor — the
    /// most discriminating feature — as the global statistic, controlling dimension-wise
    /// multiplicity without a Bonferroni penalty.
    /// </summary>
    public static class MxPbf
    {
        /// <summary>
        /// Log mxPBF for the two-sample mean test: H0: mu1 = mu2 vs H1: mu1 != mu2,
        /// in high-dimensional settings.
        /// </summary>
        /// <param name="x">Data matrix 1 (n1 x p)</param>
        /// <param name="y">Data matrix 2 (n2 x p)</param>
        /// <param name="alpha">Hyperparameter for the prior variance scaling, default 2.0</param>
        /// <returns>The log of the maximum pairwise Bayes factor (mxPBF)</returns>
        public static double TwoSampleMean(double[][] x, double[][] y, double alpha = 2.0)
        {
            int n1 = x.Length;
            int n2 = y.Length;
            if (n1 == 0 || n2 == 0) throw new ArgumentException("Empty matrices");
            int p = x[0].Length;
            if (p != y[0].Length) throw new ArgumentException("Dimension mismatch between x and y");

            int n = n1 + n2;
            double gamma = Math.Pow(Math.Max(n, p), -alpha);
            double priorTerm = 0.5 * Math.Log(gamma / (1.0 + gamma));

            double maxLogBF = double.NegativeInfinity;

            for (int j = 0; j < p; j++)
            {
                double sumX = 0, sumY = 0;
                for (int i = 0; i < n1; i++) sumX += x[i][j];
                for (int i = 0; i < n2; i++) sumY += y[i][j];

                double meanX = sumX / n1;
                double meanY = sumY / n2;
                double meanZ = (sumX + sumY) / n;

                double varX = 0, varY = 0, varZ = 0;

                for (int i = 0; i < n1; i++)
                {
                    double val = x[i][j];
                    varX += (val - meanX) * (val - meanX);
                    varZ += (val - meanZ) * (val - meanZ);
                }

                for (int i = 0; i < n2; i++)
                {
                    double val = y[i][j];
                    varY += (val - meanY) * (val - meanY);
                    varZ += (val - meanZ) * (val - meanZ);
                }

                if (varZ == 0) continue;

                double denom = varX + varY;
                if (denom == 0) denom = 1e-12;

                double dataTerm = (n / 2.0) * Math.Log(varZ / denom);
                double logBF = priorTerm + dataTerm;

                if (logBF > maxLogBF)
                {
                    maxLogBF = logBF;
                }
            }

            return maxLogBF;
        }

        /// <summary>
        /// Log mxPBF for the two-sample covariance test: H0: Σ₁ = Σ₂ vs H1: Σ₁ ≠ Σ₂ (paper §3.1).
        /// For each ordered pair (i, j), i ≠ j, the covariance comparison is reparametrized as a
        /// regression of coordinate i on coordinate j (Lee et al. 2021): the residual variances τ̂
        /// carry the partial-covariance signal, under an Inverse-Gamma(a₀, b₀) prior on each.
        /// </summary>
        /// <remarks>
        /// Implementation uses Gram-matrix factorization: G_X = X̃ᵀX̃, G_Y = ỸᵀỸ are computed once
        /// via SIMD-accelerated dot products on column-major transposed data. G_Z = G_X + G_Y is free
        /// (stacked Z = [X;Y] has column inner products that sum). Each τ̂(i,j) then becomes O(1)
        /// arithmetic from three Gram entries, eliminating O(p²n) redundant inner products.
        /// </remarks>
        /// <param name="x">Data matrix 1 (n1 x p)</param>
        /// <param name="y">Data matrix 2 (n2 x p)</param>
        /// <param name="alpha">Prior variance scaling exponent: γ = (n ∨ p)^(−alpha). Default 2.0.</param>
        /// <param name="a0">Inverse-Gamma prior shape. PROVISIONAL default — pin against R oracle.</param>
        /// <param name="b0">Inverse-Gamma prior scale. PROVISIONAL default.</param>
        /// <returns>The log of the maximum pairwise Bayes factor (mxPBF)</returns>
        public static double TwoSampleCovariance(
            double[][] x, double[][] y, double alpha = 2.0, double a0 = 1.0, double b0 = 1.0)
        {
            int n1 = x.Length;
            int n2 = y.Length;
            if (n1 == 0 || n2 == 0) throw new ArgumentException("Empty matrices");
            int p = x[0].Length;
            if (p != y[0].Length) throw new ArgumentException("Dimension mismatch between x and y");
            int n = n1 + n2;

            // 1. Column means for centering.
            double[] meanX = ColumnMeans(x, p);
            double[] meanY = ColumnMeans(y, p);

            // 2. Center + transpose to column-major: O(np) pass each, cache-friendly for Gram.
            double[][] colsX = TransposeCentered(x, n1, p, meanX);
            double[][] colsY = TransposeCentered(y, n2, p, meanY);

            // 3. Gram matrices via SIMD dot products: G[i,j] = Dot(col_i, col_j).
            double[,] gramX = MatrixOps.ColumnGramMatrix(colsX, n1, p);
            double[,] gramY = MatrixOps.ColumnGramMatrix(colsY, n2, p);

            // 4. Pooled Gram: Z = [X;Y] stacked → column inner products sum: G_Z = G_X + G_Y.
            //    (Z is centered by [meanX;meanY] weighted, but τ̂ formula uses centered inner products.)
            //    Actually, we need G_Z on centered Z. Since colsX/colsY are centered, we sum.
            double[,] gramZ = new double[p, p];
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                    gramZ[i, j] = gramX[i, j] + gramY[i, j];

            // 5. Center gramZ properly: Z centered mean = (n1*meanX + n2*meanY)/n, but since
            //    colsX and colsY are already centered at their own means, we need to adjust.
            //    For now, compute pooled tau from gramZ as-is (assumes proper centering handled).
            //    NOTE: The pooled centering adjustment: E[Z_i Z_j] with global center is:
            //    GramZ[i,j] - n * meanZ[i] * meanZ[j] where meanZ = weighted mean.
            //    However, tau uses centered-by-meanZ columns. We need centered Gram.
            //    Simpler: directly compute pooled centered Gram from stacked centered data.
            //    Optimization: adjust gramZ by centering correction.
            double[] meanZ = new double[p];
            for (int j = 0; j < p; j++) meanZ[j] = (meanX[j] * n1 + meanY[j] * n2) / n;

            // Adjust gramZ for global centering: G_Z_centered[i,j] = G_Z[i,j] - n * meanZ[i] * meanZ[j]
            // Only diagonal and below/above needed since symmetric.
            for (int i = 0; i < p; i++)
            {
                for (int j = i; j < p; j++)
                {
                    double centered = gramZ[i, j] - n * meanZ[i] * meanZ[j];
                    gramZ[i, j] = gramZ[j, i] = centered;
                }
            }

            double gamma = Math.Pow(Math.Max(n, p), -alpha);

            // Terms constant across pairs.
            double constTerm =
                -0.5 * Math.Log(gamma / (1.0 + gamma))
                + SpecialFunctions.LogGamma(n1 / 2.0 + a0)
                + SpecialFunctions.LogGamma(n2 / 2.0 + a0)
                - SpecialFunctions.LogGamma(n / 2.0 + a0)
                + a0 * Math.Log(b0) - SpecialFunctions.LogGamma(a0);

            double maxLogBF = double.NegativeInfinity;

            // 6. O(p²) scalar loop: tau from Gram entries.
            for (int i = 0; i < p; i++)
            {
                for (int j = 0; j < p; j++)
                {
                    if (i == j) continue;

                    double tau1 = ResidualVarFromGram(gramX, i, j, n1);
                    double tau2 = ResidualVarFromGram(gramY, i, j, n2);
                    double tauP = ResidualVarFromGram(gramZ, i, j, n);

                    double logBF = constTerm
                        - (n1 / 2.0 + a0) * Math.Log(b0 + (n1 / 2.0) * tau1)
                        - (n2 / 2.0 + a0) * Math.Log(b0 + (n2 / 2.0) * tau2)
                        + (n / 2.0 + a0) * Math.Log(b0 + (n / 2.0) * tauP);

                    if (logBF > maxLogBF) maxLogBF = logBF;
                }
            }

            return maxLogBF;
        }

        /// <summary>
        /// Computes residual variance τ̂(i,j) = (1/n) · (‖x̃ᵢ‖² − (x̃ᵢᵀx̃ⱼ)²/‖x̃ⱼ‖²) from Gram matrix entries.
        /// </summary>
        private static double ResidualVarFromGram(double[,] gram, int i, int j, int n)
        {
            double gii = gram[i, i];
            double gjj = gram[j, j];
            double gij = gram[i, j];
            double rss = gjj > 1e-12 ? gii - gij * gij / gjj : gii;
            return rss / n;
        }

        private static double[] ColumnMeans(double[][] data, int p)
        {
            int m = data.Length;
            var mean = new double[p];
            for (int k = 0; k < m; k++)
                for (int j = 0; j < p; j++) mean[j] += data[k][j];
            for (int j = 0; j < p; j++) mean[j] /= m;
            return mean;
        }

        /// <summary>
        /// Transposes row-major data to column-major, centering each column by subtracting its mean.
        /// Result: p columns, each a contiguous double[n] suitable for SIMD dot products.
        /// </summary>
        private static double[][] TransposeCentered(double[][] data, int n, int p, double[] mean)
        {
            var cols = new double[p][];
            for (int j = 0; j < p; j++)
            {
                var col = new double[n];
                double mu = mean[j];
                for (int i = 0; i < n; i++)
                    col[i] = data[i][j] - mu;
                cols[j] = col;
            }
            return cols;
        }
    }
}
