// ============================================================================
// Optimization/Irls.cs
// ============================================================================
// Unified IRLS solver with three-axis dispatch:
//   Axis 1 — TLoss.IsClosedForm   → closed-form Karcher (L²)
//   Axis 2 — TManifold.IsFlat     → ambient weighted average vs tangent round-trip
//   Axis 3 — HybridMode / distance threshold → Weiszfeld vs projected subgradient
//
// Changes from previous version:
//   - Distance() calls now use rented array buffers (not Span aliases) so the
//     'in double[]' signature is satisfied.
//   - Subgradient sign fixed: accumulate +w (descent toward data), then step
//     along -subgradient (standard gradient descent).
//   - Added missing using directives.
// ============================================================================
using System;
using System.Buffers;
using System.Numerics.Tensors;
using Maths.Geometry;

namespace Maths.Geometry.Solver
{
    public static class Irls
    {
        /// <summary>
        /// Hot-path IRLS solver. <paramref name="destination"/> is mutated in place;
        /// supply a warm-start initial value before calling.
        /// <para>
        /// Pass a <paramref name="finalIrlsWeights"/> span of length data.Length to
        /// capture the converged per-point IRLS weights for downstream scatter
        /// computation. Pass <c>default</c> to skip capture.
        /// </para>
        /// </summary>
        public static void Solve<TManifold, TLoss>(
            TManifold manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> weights,
            double[] destination,
            IrlsOptions opts,
            Span<double> finalIrlsWeights = default)
            where TManifold : struct, IRiemannianManifold
            where TLoss : struct, IRobustLoss
        {
            // ── Sentinel: default(IrlsOptions) zeroes MaxIterations → apply defaults ──
            if (opts.MaxIterations == 0) opts = IrlsOptions.Default;

            // ── Closed-form short-circuit (L²) ────────────────────────────────
            if (TLoss.IsClosedForm)
            {
                KarcherStep<TManifold>(manifold, data, weights, destination, opts);

                if (!finalIrlsWeights.IsEmpty)
                {
                    // L² weights are all 1 — populate for API symmetry.
                    for (int i = 0; i < data.Length; i++)
                        finalIrlsWeights[i] = weights[i];
                }
                return;
            }

            // ── Iterative IRLS path ───────────────────────────────────────────
            int n = data.Length;
            int dim = manifold.Dimension;

            double[] irlsWArr = ArrayPool<double>.Shared.Rent(n);
            double[] tangArr = ArrayPool<double>.Shared.Rent(dim);
            double[] logArr = ArrayPool<double>.Shared.Rent(dim);
            double[] nextArr = ArrayPool<double>.Shared.Rent(dim);

            try
            {
                Span<double> wirls = irlsWArr.AsSpan(0, n);
                Span<double> tangSum = tangArr.AsSpan(0, dim);
                Span<double> logV = logArr.AsSpan(0, dim);
                Span<double> next = nextArr.AsSpan(0, dim);

                ReadOnlySpan<double> current = destination;

                for (int iter = 0; iter < opts.MaxIterations; iter++)
                {
                    // ── Step 1: IRLS weights ──────────────────────────────────
                    double minDist = double.MaxValue;
                    int coincidentIdx = -1;

                    for (int i = 0; i < n; i++)
                    {
                        double r = manifold.Distance(current, data[i]);
                        if (r < minDist) { minDist = r; }

                        if (r <= opts.Epsilon)
                        {
                            wirls[i] = 0.0;
                            coincidentIdx = i;
                            continue;
                        }

                        double rUsed = (TLoss.IsSingularAtZero &&
                                        opts.SingularityPolicy == SingularityPolicy.Regularise)
                                     ? Math.Max(r, opts.Epsilon)
                                     : r;

                        wirls[i] = weights[i] * TLoss.Weight(rUsed);
                    }

                    // ── Step 1.5: optimality check at exact coincidence ───────
                    if (TLoss.IsSingularAtZero
                        && coincidentIdx >= 0
                        && opts.SingularityPolicy == SingularityPolicy.OptimalityCheck
                        && opts.HybridMode != HybridMode.SubgradientOnly)
                    {
                        if (CheckOptimality<TManifold>(
                                manifold, data, wirls, current, coincidentIdx,
                                weights[coincidentIdx], tangSum, logV))
                        {
                            // destination already holds the median; capture weights and exit.
                            if (!finalIrlsWeights.IsEmpty) wirls.CopyTo(finalIrlsWeights);
                            return;
                        }
                        // Not optimal: coincident weight stays 0, iteration drifts away.
                    }

                    // ── Step 2: update step ───────────────────────────────────
                    bool useSubgradient =
                        opts.HybridMode == HybridMode.SubgradientOnly ||
                        (opts.HybridMode == HybridMode.Hybrid
                         && TLoss.IsSingularAtZero
                         && minDist < opts.SubgradientThreshold);

                    if (!useSubgradient)
                    {
                        if (TManifold.IsFlat)
                            AmbientWeightedAverage(data, wirls, current, next);
                        else
                            TangentRoundTrip<TManifold>(manifold, data, wirls, current, tangSum, logV, next);
                    }
                    else
                    {
                        // Decaying step size: η_k = η₀ / √(k+1)
                        double eta = opts.Eta0 / Math.Sqrt(iter + 1);
                        SubgradientStep<TManifold, TLoss>(
                            manifold, data, weights, current, eta, tangSum, logV, next);
                    }

                    // ── Step 3: convergence ───────────────────────────────────
                    double shift = manifold.Distance(current, next);
                    next.CopyTo(destination);
                    current = destination;  // re-alias after copy

                    if (HasConverged(shift, destination, opts)) break;
                }

                // Recompute weights at the converged location so scatter reflects
                // the final position, not the pre-update position of the last iteration.
                if (!finalIrlsWeights.IsEmpty)
                {
                    ReadOnlySpan<double> converged = destination;
                    for (int i = 0; i < n; i++)
                    {
                        double r = manifold.Distance(converged, data[i]);
                        if (r <= opts.Epsilon) { wirls[i] = 0.0; continue; }

                        double rUsed = (TLoss.IsSingularAtZero &&
                                        opts.SingularityPolicy == SingularityPolicy.Regularise)
                                     ? Math.Max(r, opts.Epsilon) : r;
                        wirls[i] = weights[i] * TLoss.Weight(rUsed);
                    }
                    wirls.CopyTo(finalIrlsWeights);
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(irlsWArr);
                ArrayPool<double>.Shared.Return(tangArr);
                ArrayPool<double>.Shared.Return(logArr);
                ArrayPool<double>.Shared.Return(nextArr);
            }
        }

        // ── Path A: flat ambient weighted average ─────────────────────────────
        private static void AmbientWeightedAverage(
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> w,
            ReadOnlySpan<double> fallback,
            Span<double> dst)
        {
            dst.Clear();
            double sumW = 0;
            for (int i = 0; i < data.Length; i++)
            {
                double wi = w[i];
                if (wi == 0) continue;
                sumW += wi;
                var xi = data[i].AsSpan();
                for (int d = 0; d < dst.Length; d++) dst[d] += wi * xi[d];
            }
            if (sumW > 0)
            {
                double inv = 1.0 / sumW;
                for (int d = 0; d < dst.Length; d++) dst[d] *= inv;
            }
            else
            {
                // All weights zero (degenerate configuration) — hold current position.
                fallback.CopyTo(dst);
            }
        }

        // ── Path B: tangent-space round-trip ──────────────────────────────────
        private static void TangentRoundTrip<TManifold>(
            TManifold manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> w,
            ReadOnlySpan<double> pCurrent,
            Span<double> tangSum,
            Span<double> logBuf,
            Span<double> dst)
            where TManifold : struct, IRiemannianManifold
        {
            tangSum.Clear();
            double sumW = TensorPrimitives.Sum(w);
            if (sumW == 0) { pCurrent.CopyTo(dst); return; }
            double inv = 1.0 / sumW;

            for (int i = 0; i < data.Length; i++)
            {
                if (w[i] == 0) continue;
                manifold.LogMap(pCurrent, data[i], logBuf);
                manifold.AddScaled(tangSum, logBuf, w[i] * inv);
            }
            manifold.ExpMap(pCurrent, tangSum, dst);
        }

        // ── Subgradient step ──────────────────────────────────────────────────
        // FIX: accumulate +w (descent direction toward data), then step along
        //      -subgradient (standard gradient descent on the L¹ objective).
        private static void SubgradientStep<TManifold, TLoss>(
            TManifold manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> outerWeights,
            ReadOnlySpan<double> pCurrent,
            double eta,
            Span<double> subgrad,
            Span<double> logBuf,
            Span<double> dst)
            where TManifold : struct, IRiemannianManifold
            where TLoss : struct, IRobustLoss
        {
            subgrad.Clear();
            for (int i = 0; i < data.Length; i++)
            {
                manifold.LogMap(pCurrent, data[i], logBuf);
                double r = manifold.Norm(pCurrent, logBuf);
                if (r < 1e-12) continue;

                // Accumulate the descent direction: +w(r) * log_p(x_i)
                double w = outerWeights[i] * TLoss.Weight(r);
                manifold.AddScaled(subgrad, logBuf, w);
            }

            // Step: p_{k+1} = exp_p( η * subgrad )  (subgrad already points toward data)
            logBuf.Clear();
            for (int d = 0; d < subgrad.Length; d++) logBuf[d] = eta * subgrad[d];
            manifold.ExpMap(pCurrent, logBuf, dst);
        }

        // ── Closed-form Karcher step (L²) ─────────────────────────────────────
        private static void KarcherStep<TManifold>(
            TManifold manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> weights,
            double[] destination,
            IrlsOptions opts)
            where TManifold : struct, IRiemannianManifold
        {
            int dim = manifold.Dimension;
            double[] tangArr = ArrayPool<double>.Shared.Rent(dim);
            double[] logArr = ArrayPool<double>.Shared.Rent(dim);
            double[] nextArr = ArrayPool<double>.Shared.Rent(dim);
            try
            {
                Span<double> tang = tangArr.AsSpan(0, dim);
                Span<double> logV = logArr.AsSpan(0, dim);
                Span<double> next = nextArr.AsSpan(0, dim);

                ReadOnlySpan<double> current = destination;
                int maxIter = TManifold.IsFlat ? 1 : opts.MaxIterations;

                for (int iter = 0; iter < maxIter; iter++)
                {
                    if (TManifold.IsFlat)
                        AmbientWeightedAverage(data, weights, current, next);
                    else
                        TangentRoundTrip<TManifold>(manifold, data, weights, current, tang, logV, next);

                    double shift = manifold.Distance(current, next);
                    next.CopyTo(destination);
                    current = destination;
                    if (TManifold.IsFlat || HasConverged(shift, destination, opts)) break;
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(tangArr);
                ArrayPool<double>.Shared.Return(logArr);
                ArrayPool<double>.Shared.Return(nextArr);
            }
        }

        // ── Optimality check at exact coincidence ─────────────────────────────
        // Tests: ||Σ_{j≠i} wirls[j] · log_p(x_j)|| ≤ w_coincident.
        // wirls[j] already encodes weights[j]/r_j so log_p(x_j) is unnormalized here.
        private static bool CheckOptimality<TManifold>(
            TManifold manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> wirls,
            ReadOnlySpan<double> pCurrent,
            int coincidentIdx,
            double wCoincident,
            Span<double> gradBuf,
            Span<double> logBuf)
            where TManifold : struct, IRiemannianManifold
        {
            gradBuf.Clear();
            for (int j = 0; j < data.Length; j++)
            {
                if (j == coincidentIdx || wirls[j] == 0) continue;
                manifold.LogMap(pCurrent, data[j], logBuf);
                manifold.AddScaled(gradBuf, logBuf, wirls[j]);
            }
            double gradNorm = manifold.Norm(pCurrent, gradBuf);
            return gradNorm <= wCoincident;
        }

        // ── Convergence test ──────────────────────────────────────────────────
        private static bool HasConverged(
            double shift,
            ReadOnlySpan<double> p,
            IrlsOptions opts)
        {
            switch (opts.ConvergenceCriterion)
            {
                case ConvergenceCriterion.Absolute:
                    return shift < opts.Tolerance;
                case ConvergenceCriterion.RelativeToNorm:
                    double pNorm = 0;
                    for (int i = 0; i < p.Length; i++) pNorm += p[i] * p[i];
                    return shift < opts.Tolerance * (1.0 + Math.Sqrt(pNorm));
                default:
                    return shift < opts.Tolerance;
            }
        }
    }
}
