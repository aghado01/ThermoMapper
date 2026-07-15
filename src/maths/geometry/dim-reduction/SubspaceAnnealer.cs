using System;
using System.Threading;
using Maths.Rng;
using Maths.LinAlg;

namespace Maths.Geometry.DimReduction
{
    /// <summary>
    /// A black-box objective on a candidate subspace: takes an orthonormal targetDim×d projection
    /// (each row a basis vector of the subspace) and returns a scalar to be minimized.
    /// </summary>
    public delegate double SubspaceObjectiveFunction(double[][] projectionMatrix);

    /// <summary>
    /// Simulated-annealing search over the Grassmann manifold Gr(k, d) for a k-dimensional subspace
    /// (returned as an orthonormal k×d projection) that minimizes an arbitrary caller-supplied
    /// <see cref="SubspaceObjectiveFunction"/>. The objective is a black box, so the engine is agnostic
    /// to what "good" means; proposals move along Grassmann geodesics via
    /// <see cref="GrassmannManifold.ExpMap"/>.
    ///
    /// <para>This is the engine extracted from SPRED (Shape-Preserving Dimensionality Reduction,
    /// Kisung You, arXiv:2106.02096). SPRED is recovered by supplying a persistent-homology objective —
    /// the Wasserstein distance between the barcodes of the projected and the ambient cloud. That
    /// objective, and the SPRED driver that wires it in, are a consumer sitting above this engine and
    /// belong with the barcodes in TDA.Ph; the annealer itself carries no homology dependency.</para>
    /// </summary>
    public static class SubspaceAnnealer
    {
        /// <summary>Anneal a projection that minimizes <paramref name="objective"/>.</summary>
        /// <param name="data">Row-major samples; every row has the same ambient dimension d.</param>
        /// <param name="targetDim">Subspace dimension k (1 ≤ k ≤ d), e.g. 2 or 3 for visualization.</param>
        /// <param name="objective">Scalar objective value of a candidate k×d orthonormal projection.</param>
        /// <param name="maxIters">Number of simulated-annealing steps.</param>
        /// <param name="seed">RNG seed for a reproducible annealing stream; null draws OS entropy.</param>
        /// <param name="cancellationToken">Cancellation observed before setup and between annealing steps.</param>
        /// <returns>The best k×d orthonormal projection found (rows orthonormal).</returns>
        public static double[][] Compute(
            double[][] data, int targetDim, SubspaceObjectiveFunction objective,
            int maxIters = 1000, int? seed = null, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            int nSamples = data.Length;
            if (nSamples == 0) throw new ArgumentException("Empty data", nameof(data));
            int d = data[0].Length;
            int k = targetDim;
            if (k < 1 || k > d)
                throw new ArgumentOutOfRangeException(nameof(targetDim),
                    "targetDim must satisfy 1 ≤ targetDim ≤ d.");

            var manifold = new GrassmannManifold(ambientN: d, subspaceR: k);

            // Warm start on the PCA subspace, packed as a d×k column-major Grassmann representative
            // (each PCA component becomes one orthonormal column / one projection row).
            var pca = Pca.Compute(data, k, center: true, whiten: false);
            cancellationToken.ThrowIfCancellationRequested();
            double[] current = PackColumnMajor(pca.Components, k, d);
            MatrixOps.Orthonormalize(current, d, k);   // clean Stiefel representative for a valid point

            double[][] currentProj = ToProjection(current, k, d);
            double currentValue = objective(currentProj);
            cancellationToken.ThrowIfCancellationRequested();

            double[][] bestProj = currentProj;
            double bestValue = currentValue;

            var rng = new Xoshiro256PlusPlus(seed);
            double[] tangent = new double[d * k];
            double[] proposal = new double[d * k];
            const double initialTemp = 1.0;

            for (int iter = 0; iter < maxIters; iter++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                double temp = initialTemp * Math.Pow(0.99, iter);   // geometric cooling
                double step = temp * 0.1;                            // geodesic step length

                HorizontalTangent(current, d, k, rng, step, tangent);
                manifold.ExpMap(current, tangent, proposal);         // retract along the geodesic

                double[][] proposalProj = ToProjection(proposal, k, d);
                double proposalValue = objective(proposalProj);
                cancellationToken.ThrowIfCancellationRequested();

                if (proposalValue < currentValue ||
                    rng.NextDouble() < Math.Exp((currentValue - proposalValue) / temp))
                {
                    (current, proposal) = (proposal, current);       // adopt; recycle the old buffer
                    currentProj = proposalProj;
                    currentValue = proposalValue;

                    if (currentValue < bestValue)
                    {
                        bestProj = currentProj;
                        bestValue = currentValue;
                    }
                }
            }

            return bestProj;
        }

        /// <summary>
        /// Random horizontal tangent at the subspace <paramref name="y"/> (d×k column-major), scaled
        /// to Frobenius norm <paramref name="step"/> so <see cref="GrassmannManifold.ExpMap"/> travels
        /// that geodesic distance. Isotropic ambient Gaussian, projected off Y's span (Δ ← Δ − Y(YᵀΔ)).
        /// </summary>
        private static void HorizontalTangent(
            double[] y, int d, int k, Xoshiro256PlusPlus rng, double step, double[] dst)
        {
            for (int i = 0; i < d * k; i++) dst[i] = NextGaussian(rng);

            // YᵀΔ (k×k), then Δ -= Y·(YᵀΔ) to strip the vertical (in-span) component.
            double[] ytd = new double[k * k];
            for (int b = 0; b < k; b++)
                for (int a = 0; a < k; a++)
                {
                    double s = 0.0;
                    int ya = a * d, db = b * d;
                    for (int row = 0; row < d; row++) s += y[ya + row] * dst[db + row];
                    ytd[a + b * k] = s;
                }
            for (int b = 0; b < k; b++)
                for (int a = 0; a < k; a++)
                {
                    double coeff = ytd[a + b * k];
                    int ya = a * d, db = b * d;
                    for (int row = 0; row < d; row++) dst[db + row] -= coeff * y[ya + row];
                }

            double norm = 0.0;
            for (int i = 0; i < d * k; i++) norm += dst[i] * dst[i];
            norm = Math.Sqrt(norm);
            double scale = norm > 1e-12 ? step / norm : 0.0;
            for (int i = 0; i < d * k; i++) dst[i] *= scale;
        }

        /// <summary>Standard normal via Box–Muller (one draw per call; deterministic on the seeded stream).</summary>
        private static double NextGaussian(Xoshiro256PlusPlus rng)
        {
            double u1 = 1.0 - rng.NextDouble();   // (0, 1], so Log is finite
            double u2 = rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }

        // k×d row-major components → d×k column-major point (row j → column j).
        private static double[] PackColumnMajor(double[][] components, int k, int d)
        {
            double[] point = new double[d * k];
            for (int j = 0; j < k; j++)
                Array.Copy(components[j], 0, point, j * d, d);
            return point;
        }

        // d×k column-major point → k×d projection (column j → row j) for the objective function.
        private static double[][] ToProjection(double[] point, int k, int d)
        {
            double[][] proj = new double[k][];
            for (int j = 0; j < k; j++)
            {
                proj[j] = new double[d];
                Array.Copy(point, j * d, proj[j], 0, d);
            }
            return proj;
        }
    }
}
