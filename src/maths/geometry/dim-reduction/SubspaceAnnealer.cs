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

    /// <summary>The outcome of a subspace anneal.</summary>
    /// <param name="Projection">The best k×d orthonormal projection found (rows orthonormal).</param>
    /// <param name="Objective">The objective value at <paramref name="Projection"/>, as tracked during
    /// the anneal; equals a fresh evaluation because the objective is a deterministic function of the
    /// projection.</param>
    public sealed record SubspaceAnnealerResult(double[][] Projection, double Objective);

    /// <summary>
    /// Simulated-annealing search over the Grassmann manifold Gr(k, d) for a k-dimensional subspace
    /// (returned as an orthonormal k×d projection) that minimizes an arbitrary caller-supplied
    /// <see cref="SubspaceObjectiveFunction"/>. The objective is a black box, so the engine is agnostic
    /// to what "good" means; proposals move along Grassmann geodesics — primarily two-plane Givens
    /// rotations (the rank-1 closed form of <see cref="GrassmannManifold.ExpMap"/>), optionally mixed
    /// with isotropic horizontal tangents per <see cref="SubspaceAnnealerOptions"/>. Geodesic step
    /// scales adapt per move coordinate toward a target acceptance rate; cooling governs the
    /// Metropolis temperature.
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
        /// <param name="seed">RNG seed for a reproducible annealing stream; null draws OS entropy.
        /// One Xoshiro stream drives proposals, Metropolis, and (indirectly) step adaptation, so
        /// same-seed runs are bit-identical.</param>
        /// <param name="options">Proposal mixture, step adaptation, and cooling; null takes the
        /// <see cref="SubspaceAnnealerOptions"/> defaults.</param>
        /// <param name="cancellationToken">Cancellation observed before setup and between annealing steps.</param>
        /// <returns>The best k×d orthonormal projection found and its objective value.</returns>
        public static SubspaceAnnealerResult Compute(
            double[][] data, int targetDim, SubspaceObjectiveFunction objective,
            int maxIters = 1000, int? seed = null,
            SubspaceAnnealerOptions? options = null,
            CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            options ??= new SubspaceAnnealerOptions();
            options.Validate();
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
            double[] direction = new double[d];

            // One adaptive step scale per move coordinate — each retained column's Givens angle,
            // plus one for the isotropic kind (index k). Column-local acceptance diverges sharply
            // once some columns converge (their moves reject at any scale while mobile columns
            // keep improving), so a pooled scale would be strangled below the target by the
            // converged majority and starve the mobile columns.
            double[] stepByMove = new double[k + 1];
            Array.Fill(stepByMove, Math.Clamp(options.InitialStep, options.StepFloor, options.StepCeiling));
            // Multiplicative controller with zero expected log-step drift exactly at the target:
            // grow^p · shrink^(1−p) = 1 ⇔ p = TargetAcceptance.
            double growOnAccept = Math.Exp(AdaptationGain * (1.0 - options.TargetAcceptance));
            double shrinkOnReject = Math.Exp(-AdaptationGain * options.TargetAcceptance);

            for (int iter = 0; iter < maxIters; iter++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                double temp = options.InitialTemperature * Math.Pow(options.CoolingRate, iter);

                int move;   // controller index: the rotated column for Givens, k for isotropic
                if (rng.NextDouble() < options.IsotropicFraction)
                {
                    move = k;
                    HorizontalTangent(current, d, k, rng, stepByMove[move], tangent);
                    manifold.ExpMap(current, tangent, proposal);     // retract along the geodesic
                }
                else
                {
                    move = rng.NextInt(k);
                    GivensProposal(current, d, k, move, rng, stepByMove[move], direction, proposal);
                }

                double[][] proposalProj = ToProjection(proposal, k, d);
                double proposalValue = objective(proposalProj);
                cancellationToken.ThrowIfCancellationRequested();

                bool accepted = proposalValue < currentValue ||
                    rng.NextDouble() < Math.Exp((currentValue - proposalValue) / temp);
                if (accepted)
                {
                    (current, proposal) = (proposal, current);       // adopt; recycle the old buffer
                    currentValue = proposalValue;

                    if (currentValue < bestValue)
                    {
                        bestProj = proposalProj;
                        bestValue = currentValue;
                    }
                }

                stepByMove[move] = accepted
                    ? Math.Min(options.StepCeiling, stepByMove[move] * growOnAccept)
                    : Math.Max(options.StepFloor, stepByMove[move] * shrinkOnReject);
            }

            return new SubspaceAnnealerResult(bestProj, bestValue);
        }

        // Step-controller gain: each decision nudges log-step by ±gain-scaled amounts, so the
        // scale equilibrates within tens of iterations without chattering.
        private const double AdaptationGain = 0.1;

        /// <summary>
        /// Two-plane Givens proposal: rotate the retained column <paramref name="col"/> of
        /// <paramref name="y"/> by angle θ = step × standard normal toward a uniformly random unit
        /// direction v orthogonal to the current span (y_i ← y_i·cosθ + v·sinθ). This is the closed
        /// form of <see cref="GrassmannManifold.ExpMap"/> for the rank-1 horizontal tangent θ·v·e_iᵀ
        /// (single singular value θ) — the classic subspace-search move. Its improving fraction
        /// survives high codimension, where a fixed-length isotropic tangent's directional
        /// derivative thins out like 1/√dim against an O(step²) curvature penalty.
        /// </summary>
        private static void GivensProposal(
            double[] y, int d, int k, int col, Xoshiro256PlusPlus rng, double step,
            double[] direction, double[] dst)
        {
            double theta = step * NextGaussian(rng);

            for (int row = 0; row < d; row++) direction[row] = NextGaussian(rng);

            // v ← v − Y(Yᵀv): strip the in-span component so the rotated column stays orthogonal
            // to every other column.
            for (int a = 0; a < k; a++)
            {
                double dot = 0.0;
                int ya = a * d;
                for (int row = 0; row < d; row++) dot += y[ya + row] * direction[row];
                for (int row = 0; row < d; row++) direction[row] -= dot * y[ya + row];
            }

            double norm = 0.0;
            for (int row = 0; row < d; row++) norm += direction[row] * direction[row];
            norm = Math.Sqrt(norm);

            Array.Copy(y, dst, d * k);
            if (norm < 1e-12) return;   // k = d: the complement is empty and Gr(d, d) is a point

            double cos = Math.Cos(theta);
            double sin = Math.Sin(theta) / norm;   // folds v's normalization into the rotation
            int yc = col * d;
            double renorm = 0.0;
            for (int row = 0; row < d; row++)
            {
                double value = y[yc + row] * cos + direction[row] * sin;
                dst[yc + row] = value;
                renorm += value * value;
            }
            // Exact unit length keeps roundoff from compounding over long anneals; the residual
            // cross-column error stays second-order.
            renorm = 1.0 / Math.Sqrt(renorm);
            for (int row = 0; row < d; row++) dst[yc + row] *= renorm;
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
