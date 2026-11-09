// ============================================================================
// Optimization/Irls.cs
// ============================================================================
// Unified IRLS solver with three-axis dispatch:
//   Axis 1 — TLoss.IsClosedForm   → closed-form Karcher (L²)
//   Axis 2 — TManifold.IsFlat     → ambient weighted average vs tangent round-trip
//   Axis 3 — HybridMode / distance threshold → Weiszfeld vs projected subgradient
// The distance-valued constants are read in the unit λ that IrlsOptions.Scale
// picks: the spread σ (default) or the data's own units.
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
        /// Stall radius under <see cref="IrlsScale.Spread"/>, per unit of the iterate's ambient norm: a few ulps. A step
        /// that moves the iterate no further has fallen below the manifold's Exp resolution or the rounding of the iterate.
        /// </summary>
        internal const double RoundingFloor = 4 * 2.220446049250313e-16;

        // A Log map more than this fraction short of a distance above CutLocusFloor was cut at the cut locus.
        private const double CutLocusGap = 1e-3;

        // Below this a Log short of the distance is rounding, not the cut locus: every curved manifold reads its
        // distances through atan2 or asinh, accurate to about √(n·r)·ε ≈ 2e-15, so at 1e-10 the 0.1% gap is clear of it.
        private const double CutLocusFloor = 1e-10;

        /// <summary>
        /// Hot-path IRLS solver. <paramref name="destination"/> is mutated in place;
        /// supply a warm-start initial value before calling.
        /// <para>
        /// The distance-valued constants of <paramref name="opts"/> are read in the unit that
        /// <see cref="IrlsOptions.Scale"/> picks: by default the spread, the weighted median distance
        /// from the iterate to the data, so the median and mean are equivariant under rescaling the data.
        /// </para>
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
            bool spread = opts.Scale == IrlsScale.Spread;

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
            double[] radiiArr = ArrayPool<double>.Shared.Rent(n);
            double[] keyArr = ArrayPool<double>.Shared.Rent(n);
            double[] keyWArr = ArrayPool<double>.Shared.Rent(n);
            double[] tangArr = ArrayPool<double>.Shared.Rent(dim);
            double[] logArr = ArrayPool<double>.Shared.Rent(dim);
            double[] nextArr = ArrayPool<double>.Shared.Rent(dim);
            // Spread keeps the iteration's n log vectors: they give the radii and the step.
            double[]? logsArr = spread ? ArrayPool<double>.Shared.Rent(n * dim) : null;
            bool[] testedArr = ArrayPool<bool>.Shared.Rent(n);

            try
            {
                Span<double> wirls = irlsWArr.AsSpan(0, n);
                Span<double> radii = radiiArr.AsSpan(0, n);
                Span<double> keys = keyArr.AsSpan(0, n);
                Span<double> keyW = keyWArr.AsSpan(0, n);
                Span<double> tangSum = tangArr.AsSpan(0, dim);
                Span<double> logV = logArr.AsSpan(0, dim);
                Span<double> next = nextArr.AsSpan(0, dim);
                Span<double> logs = spread ? logsArr.AsSpan(0, n * dim) : default;
                Span<bool> tested = testedArr.AsSpan(0, n);
                tested.Clear();

                ReadOnlySpan<double> current = destination;
                double totalWeight = spread ? TensorPrimitives.Sum(weights.Slice(0, n)) : 1.0;

                for (int iter = 0; iter < opts.MaxIterations; iter++)
                {
                    // ── Step 1: radii, the unit λ, IRLS weights ───────────────
                    // Spread reads each radius as ‖Log_p xᵢ‖, which resolves small spreads that an
                    // arccos-based Distance cannot, except at the cut locus (see Radius); Data keeps Distance.
                    if (spread)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            Span<double> li = logs.Slice(i * dim, dim);
                            manifold.LogMap(current, data[i], li);
                            radii[i] = Radius(manifold, current, data[i], li, out _);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < n; i++)
                            radii[i] = manifold.Distance(current, data[i]);
                    }

                    double unit = 1.0;
                    if (spread)
                    {
                        unit = SpreadOf(radii, weights, keys, keyW);
                        // Every weighted point sits on the iterate, which is then the median.
                        if (unit == 0) break;
                    }
                    double coincidence = opts.Epsilon * unit;

                    double minDist = double.MaxValue;
                    int nearest = -1;
                    double irlsSum = 0, nearestWeight = 0;   // nearestWeight: IRLS weight at radius minDist

                    for (int i = 0; i < n; i++)
                    {
                        double r = radii[i];
                        if (r < minDist) { minDist = r; nearest = i; nearestWeight = 0; }

                        if (r <= coincidence)
                        {
                            wirls[i] = 0.0;
                            continue;
                        }

                        double rUsed = (TLoss.IsSingularAtZero &&
                                        opts.SingularityPolicy == SingularityPolicy.Regularise)
                                     ? Math.Max(r, coincidence)
                                     : r;

                        wirls[i] = weights[i] * TLoss.Weight(rUsed);
                        irlsSum += wirls[i];
                        if (r == minDist) nearestWeight += wirls[i];
                    }

                    // ── Step 1.5: optimality check at the nearest data point ──
                    // Weiszfeld reaches a median that is a data point only in the limit, and the Hybrid switch
                    // throws it off again. Once the nearest point x_m coincides with the iterate or holds, with
                    // whatever ties its radius (its duplicates), at least half the IRLS weight, test it directly;
                    // each point is tested once.
                    if (TLoss.IsSingularAtZero
                        && nearest >= 0
                        && !tested[nearest]
                        && opts.SingularityPolicy == SingularityPolicy.OptimalityCheck
                        && opts.HybridMode != HybridMode.SubgradientOnly
                        && (minDist <= coincidence || 2.0 * nearestWeight >= irlsSum))
                    {
                        if (VertexIsOptimal<TManifold>(
                                manifold, data, weights, nearest, coincidence, tested, tangSum, logV))
                        {
                            // The median is that data point: return it, with weights recomputed there.
                            data[nearest].AsSpan(0, dim).CopyTo(destination);
                            break;
                        }
                        // Not optimal: coincident weights stay 0, and the iteration drifts away.
                    }

                    // ── Step 2: update step (left in tangSum, except the Data-unit flat average) ──
                    bool useSubgradient =
                        opts.HybridMode == HybridMode.SubgradientOnly ||
                        (opts.HybridMode == HybridMode.Hybrid
                         && TLoss.IsSingularAtZero
                         && minDist < opts.SubgradientThreshold * unit);

                    if (spread)
                    {
                        if (!useSubgradient)
                        {
                            // Weiszfeld: v = Σ ωᵢ Log_p xᵢ / Σ ω.
                            double sumW = TensorPrimitives.Sum(wirls);
                            StoredStep(manifold, logs, wirls, sumW > 0 ? 1.0 / sumW : 0.0, tangSum);
                        }
                        else
                        {
                            // Subgradient: v = λη₀/√(k+1) · Σ ωᵢ Log_p xᵢ / Σw, where ωᵢ = wᵢ w(rᵢ) already
                            // excludes the coincident points.
                            double eta = opts.Eta0 * unit / Math.Sqrt(iter + 1) / totalWeight;
                            StoredStep(manifold, logs, wirls, eta, tangSum);
                        }
                        manifold.ExpMap(current, tangSum, next);
                    }
                    else if (!useSubgradient)
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
                    // Spread measures the shift as the step's norm, and also stops when the step no
                    // longer moves the iterate (below the manifold's Exp resolution or the rounding of p).
                    bool stalled = spread && Stalled(current, next);
                    double shift = spread
                        ? manifold.Norm(current, tangSum)
                        : manifold.Distance(current, next);
                    next.CopyTo(destination);
                    current = destination;  // re-alias after copy

                    if (stalled || HasConverged(shift, destination, opts, unit)) break;
                }

                // Recompute weights at the converged location so scatter reflects
                // the final position, not the pre-update position of the last iteration.
                if (!finalIrlsWeights.IsEmpty)
                {
                    ReadOnlySpan<double> converged = destination;
                    for (int i = 0; i < n; i++)
                    {
                        if (spread)
                        {
                            manifold.LogMap(converged, data[i], logV);
                            radii[i] = Radius(manifold, converged, data[i], logV, out _);
                        }
                        else
                        {
                            radii[i] = manifold.Distance(converged, data[i]);
                        }
                    }
                    double coincidence = opts.Epsilon * (spread ? SpreadOf(radii, weights, keys, keyW) : 1.0);

                    for (int i = 0; i < n; i++)
                    {
                        double r = radii[i];
                        if (r <= coincidence) { wirls[i] = 0.0; continue; }

                        double rUsed = (TLoss.IsSingularAtZero &&
                                        opts.SingularityPolicy == SingularityPolicy.Regularise)
                                     ? Math.Max(r, coincidence) : r;
                        wirls[i] = weights[i] * TLoss.Weight(rUsed);
                    }
                    wirls.CopyTo(finalIrlsWeights);
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(irlsWArr);
                ArrayPool<double>.Shared.Return(radiiArr);
                ArrayPool<double>.Shared.Return(keyArr);
                ArrayPool<double>.Shared.Return(keyWArr);
                ArrayPool<double>.Shared.Return(tangArr);
                ArrayPool<double>.Shared.Return(logArr);
                ArrayPool<double>.Shared.Return(nextArr);
                if (logsArr is not null) ArrayPool<double>.Shared.Return(logsArr);
                ArrayPool<bool>.Shared.Return(testedArr);
            }
        }

        // The radius of x from p under Spread: ‖Log_p x‖, which holds at small distances where an arccos-based
        // Distance loses digits (relative error about ε/2θ²). The sphere's and Grassmann's Log maps drop the
        // direction of an angle at the cut locus (an antipode, a principal angle of π/2), leaving a vector short of
        // the distance, even zero. There (cut: a Log more than CutLocusGap short of a distance above CutLocusFloor)
        // the radius is the Distance, and the shortened Log stays the direction. Flat manifolds have no cut locus.
        private static double Radius<TManifold>(
            TManifold manifold,
            ReadOnlySpan<double> p,
            ReadOnlySpan<double> x,
            ReadOnlySpan<double> log,
            out bool cut)
            where TManifold : struct, IRiemannianManifold
        {
            double r = manifold.Norm(p, log);
            cut = false;
            if (TManifold.IsFlat) return r;
            double d = manifold.Distance(p, x);
            cut = d > CutLocusFloor && r < (1.0 - CutLocusGap) * d;
            return cut ? d : r;
        }

        /// <summary>
        /// The spread σ of <see cref="IrlsScale.Spread"/>: the weighted lower median of <paramref name="radii"/> over the
        /// positively weighted points; their weighted mean when at least half the weight sits at radius 0; 0 when all of
        /// it does, or when no weight is positive. <paramref name="keys"/> and <paramref name="keyWeights"/> are scratch
        /// of at least radii.Length.
        /// </summary>
        internal static double SpreadOf(
            ReadOnlySpan<double> radii,
            ReadOnlySpan<double> weights,
            Span<double> keys,
            Span<double> keyWeights)
        {
            int m = 0;
            double total = 0, weightedSum = 0;
            for (int i = 0; i < radii.Length; i++)
            {
                double w = weights[i];
                if (!(w > 0)) continue;
                keys[m] = radii[i];
                keyWeights[m] = w;
                m++;
                total += w;
                weightedSum += w * radii[i];
            }
            if (m == 0) return 0.0;

            keys.Slice(0, m).Sort(keyWeights.Slice(0, m));
            double half = 0.5 * total, cumulative = 0, median = keys[m - 1];
            for (int j = 0; j < m; j++)
            {
                cumulative += keyWeights[j];
                if (cumulative >= half) { median = keys[j]; break; }
            }
            return median > 0 ? median : weightedSum / total;
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
        // Leaves the step Σ wᵢ Log_p xᵢ / Σ w in tangSum. With radii non-empty, also writes
        // ‖Log_p xᵢ‖ for each point it visits (those of nonzero weight).
        private static void TangentRoundTrip<TManifold>(
            TManifold manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> w,
            ReadOnlySpan<double> pCurrent,
            Span<double> tangSum,
            Span<double> logBuf,
            Span<double> dst,
            Span<double> radii = default)
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
                if (!radii.IsEmpty) radii[i] = manifold.Norm(pCurrent, logBuf);
                manifold.AddScaled(tangSum, logBuf, w[i] * inv);
            }
            manifold.ExpMap(pCurrent, tangSum, dst);
        }

        // ── Spread step from the stored logs ──────────────────────────────────
        // step = scale · Σ wᵢ Log_p xᵢ over the stored logs (n blocks of Dimension).
        private static void StoredStep<TManifold>(
            TManifold manifold,
            ReadOnlySpan<double> logs,
            ReadOnlySpan<double> w,
            double scale,
            Span<double> step)
            where TManifold : struct, IRiemannianManifold
        {
            int dim = step.Length;
            step.Clear();
            for (int i = 0; i < w.Length; i++)
            {
                if (w[i] == 0) continue;
                manifold.AddScaled(step, logs.Slice(i * dim, dim), w[i] * scale);
            }
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
            bool spread = opts.Scale == IrlsScale.Spread;
            int n = data.Length;
            int dim = manifold.Dimension;
            double[] tangArr = ArrayPool<double>.Shared.Rent(dim);
            double[] logArr = ArrayPool<double>.Shared.Rent(dim);
            double[] nextArr = ArrayPool<double>.Shared.Rent(dim);
            double[] radiiArr = ArrayPool<double>.Shared.Rent(n);
            double[] keyArr = ArrayPool<double>.Shared.Rent(n);
            double[] keyWArr = ArrayPool<double>.Shared.Rent(n);
            try
            {
                Span<double> tang = tangArr.AsSpan(0, dim);
                Span<double> logV = logArr.AsSpan(0, dim);
                Span<double> next = nextArr.AsSpan(0, dim);
                // Curved Karcher flow under Spread reads σ from the log norms it already computes.
                Span<double> radii = spread && !TManifold.IsFlat ? radiiArr.AsSpan(0, n) : default;

                ReadOnlySpan<double> current = destination;
                int maxIter = TManifold.IsFlat ? 1 : opts.MaxIterations;

                for (int iter = 0; iter < maxIter; iter++)
                {
                    if (TManifold.IsFlat && !spread)
                        AmbientWeightedAverage(data, weights, current, next);
                    else
                        TangentRoundTrip<TManifold>(manifold, data, weights, current, tang, logV, next, radii);

                    double unit = 1.0;
                    double shift;
                    bool stalled = false;
                    if (spread)
                    {
                        shift = manifold.Norm(current, tang);
                        stalled = Stalled(current, next);
                        if (!radii.IsEmpty)
                            unit = SpreadOf(radii, weights, keyArr.AsSpan(0, n), keyWArr.AsSpan(0, n));
                    }
                    else
                    {
                        shift = manifold.Distance(current, next);
                    }
                    next.CopyTo(destination);
                    current = destination;
                    if (TManifold.IsFlat || stalled || HasConverged(shift, destination, opts, unit)) break;
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(tangArr);
                ArrayPool<double>.Shared.Return(logArr);
                ArrayPool<double>.Shared.Return(nextArr);
                ArrayPool<double>.Shared.Return(radiiArr);
                ArrayPool<double>.Shared.Return(keyArr);
                ArrayPool<double>.Shared.Return(keyWArr);
            }
        }

        // ── Optimality check at a data point ─────────────────────────────────
        // The L¹ condition at x_m (Kuhn 1973; Vardi and Zhang 2000): x_m is a median iff
        // ‖Σ_{d(x_m,x_j) > ε} w_j Log_{x_m} x_j / ‖Log_{x_m} x_j‖‖ ≤ Σ_{d(x_m,x_j) ≤ ε} w_j, the weight sitting on x_m,
        // duplicates included. A point on x_m's cut locus has no single direction; each distance is 1-Lipschitz, so
        // such a point enters at its full weight, ‖G‖ + Σ_cut w_j ≤ the weight on x_m, which still certifies a local
        // minimum. Marks x_m and its duplicates as tested.
        private static bool VertexIsOptimal<TManifold>(
            TManifold manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double> weights,
            int m,
            double coincidence,
            Span<bool> tested,
            Span<double> gradBuf,
            Span<double> logBuf)
            where TManifold : struct, IRiemannianManifold
        {
            ReadOnlySpan<double> xm = data[m];
            gradBuf.Clear();
            double weightOn = 0, cutWeight = 0;
            for (int j = 0; j < data.Length; j++)
            {
                double w = weights[j];
                manifold.LogMap(xm, data[j], logBuf);
                double r = Radius(manifold, xm, data[j], logBuf, out bool cut);
                if (j == m || r <= coincidence)
                {
                    weightOn += w;
                    tested[j] = true;
                    continue;
                }
                if (w == 0) continue;
                if (cut) { cutWeight += w; continue; }
                manifold.AddScaled(gradBuf, logBuf, w / r);
            }
            return manifold.Norm(xm, gradBuf) + cutWeight <= weightOn;
        }

        // ── Convergence test ──────────────────────────────────────────────────
        // Absolute: shift < τλ. RelativeToNorm: shift < τ(λ + ‖p‖).
        private static bool HasConverged(
            double shift,
            ReadOnlySpan<double> p,
            IrlsOptions opts,
            double unit)
        {
            switch (opts.ConvergenceCriterion)
            {
                case ConvergenceCriterion.RelativeToNorm:
                    double pNorm = 0;
                    for (int i = 0; i < p.Length; i++) pNorm += p[i] * p[i];
                    return shift < opts.Tolerance * (unit + Math.Sqrt(pNorm));
                default:
                    return shift < opts.Tolerance * unit;
            }
        }

        // The iterate moved by no more than RoundingFloor·‖p‖ in ambient coordinates.
        private static bool Stalled(ReadOnlySpan<double> p, ReadOnlySpan<double> next)
        {
            double move = 0, norm = 0;
            for (int i = 0; i < next.Length; i++)
            {
                double d = next[i] - p[i];
                move += d * d;
                norm += p[i] * p[i];
            }
            return move <= RoundingFloor * RoundingFloor * norm;
        }
    }
}
