// ============================================================================
// Coupling/BandwidthEstimation.cs
// ============================================================================
// Bandwidth (Î´) estimation for coupling kernels from a span of scalar edge
// distances (e.g. k-NN edge lengths per node).
//
// Three layers:
//   Per-kernel routes    â€” ForGaussian / ForLaplacian / ForCauchy / ForLinear
//                          + ForMixture. Each accepts an optional
//                          BandwidthStrategy that selects how the scale
//                          is extracted from the NN distance sample before
//                          the kernel's consistency factor is applied:
//                            MadConsistencyFactor â€” raw MAD (Euclidean default)
//                            QuantileNormalized   â€” MAD on q95-clamped sample
//                                                   (bounded metrics: JSD,
//                                                   FisherRao/Simplex, Cosine)
//                            LogScaleHyperbolic   â€” back-transformed log-MAD
//                                                   (PoincarÃ©, hyperbolic)
//   Agnostic primitives  â€” Mean, Median, Mad, Max
//                          exposed for callers who need custom routing.
// ============================================================================
using System;
using Graphs.Coupling;
using Graphs.Distance;
using Maths.Geometry.Bandwidth;

namespace Graphs
{
    public static class BandwidthEstimation
    {
        // â”€â”€ Consistency factors â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        // Scale raw MAD into an unbiased scale estimate under each distribution.
        //   Gaussian:  1 / Î¦â»Â¹(0.75)  â‰ˆ 1.4826
        //   Laplacian: 1 / ln(2)       â‰ˆ 1.4427
        //   Cauchy:    MAD(Cauchy(0,Î³)) = Î³ exactly â†’ 1.0
        public const double GaussianFactor = 1.4826;
        public const double LaplacianFactor = 1.4427;
        public const double CauchyFactor = 1.0;
        // Hyperbolic (GeodesicLinear) MADâ†’scale factors, per kernel. Derived in
        // opus-hyperbolic-bandwidth-factors.md against the radial reference
        // p(r) âˆ K(r)Â·sinh^{d-1}(r) (You 2604.24895). Dimension is absorbed by the
        // empirical log-median/log-MAD, so these are dimension-free.
        public const double GaussianHyperbolicFactor = 1.4826;
        // MAD consistency factor for Laplacian kernel on hyperbolic space H^d.
        // Hyperbolic-volume-derived and regime-conditional; do not treat it
        // as a universal invariant or transport it unchanged to other curvatures.
        // Derived via moment-matching on the intrinsic Gaussian reference measure
        // p(r) âˆ K(r)Â·sinh^{d-1}(r) [You 2604.24895], with computation in
        // Hyperbolic.MatchBeta().
        // Unlike 1.4826 for Gaussian, it is tied to the hyperbolic volume-growth
        // regime (Î² > d - 1) and should not be reused as a cross-curvature constant.
        // The live moment-matching path ForIntrinsicGaussianHyperbolic() provides
        // adaptive refinement for the exact data sample.
        public const double LaplacianHyperbolicFactor = 1.67;
        // Intentional quantile fallback: the hyperbolic Cauchy radial reference is
        // non-normalizable, so no moment-matched scale factor exists. factor = 0
        // yields Î´ = exp(logMedian) = median(R).
        public const double CauchyHyperbolicFactor = 0.0;
        // Spherical MADâ†’scale factors, per kernel.
        public const double GaussianSphericalFactor = 1.4826;
        public const double LaplacianSphericalFactor = 1.34;
        public const double CauchySphericalFactor = 1.0;

        /// <summary>
        /// Default quantile for <see cref="BandwidthStrategy.QuantileNormalized"/>'s
        /// outlier clamping â€” clamps NN distances above this quantile to the
        /// quantile value before computing MAD. Chosen as a fixed statistical
        /// convention rather than a free parameter (per Pass 3 design notes).
        /// </summary>
        public const double QuantileClampLevel = 0.95;

        /// <summary>
        /// Floor applied before <c>log()</c> so coincident points (distance 0)
        /// do not produce <c>log(0)</c> = -âˆž. Pure log space: no additive offset.
        /// </summary>
        private const double LogScaleFloor = 1e-12;

        // â”€â”€ Per-kernel routes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

        public static double ForGaussian(
            ReadOnlySpan<double> distances, Span<double> scratch,
            BandwidthStrategy strategy = BandwidthStrategy.MadConsistencyFactor,
            SpaceGeometry geometry = SpaceGeometry.Euclidean,
            double fallback = 1.0)
            => distances.Length == 0 ? fallback : BandwidthFor(distances, scratch, strategy, GaussianFactor, GaussianHyperbolicFactor, GaussianSphericalFactor, geometry);

        public static double ForIntrinsicGaussianHyperbolic(
            ReadOnlySpan<double> distances,
            int ambientDimension,
            double fallback = 1.0)
        {
            if (distances.Length == 0)
                return fallback;
            if (ambientDimension < 1)
                throw new ArgumentOutOfRangeException(nameof(ambientDimension), ambientDimension, "Ambient dimension must be positive.");

            double observedSecondMoment = MeanSquare(distances, fallback * fallback);
            double beta = Hyperbolic.MatchBeta(observedSecondMoment, ambientDimension);
            return HeatKernel.BandwidthFromBeta(beta);
        }

        public static double ForIntrinsicGaussianSpherical(
            ReadOnlySpan<double> distances,
            int intrinsicDimension,
            double fallback = 1.0)
        {
            if (distances.Length == 0)
                return fallback;
            if (intrinsicDimension < 1)
                throw new ArgumentOutOfRangeException(nameof(intrinsicDimension), intrinsicDimension, "Intrinsic dimension must be positive.");

            double observedSecondMoment = MeanSquare(distances, fallback * fallback);
            double beta = Spherical.MatchBeta(observedSecondMoment, intrinsicDimension);
            return HeatKernel.BandwidthFromBeta(beta);
        }

        public static double IntrinsicHeatTimeFromBandwidth(double bandwidth)
        {
            if (!double.IsFinite(bandwidth) || bandwidth <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(bandwidth), bandwidth, "Bandwidth must be finite and positive.");

            return 0.5 * bandwidth * bandwidth;
        }

        public static double ForLaplacian(
            ReadOnlySpan<double> distances, Span<double> scratch,
            BandwidthStrategy strategy = BandwidthStrategy.MadConsistencyFactor,
            SpaceGeometry geometry = SpaceGeometry.Euclidean,
            double fallback = 1.0)
            => distances.Length == 0 ? fallback : BandwidthFor(distances, scratch, strategy, LaplacianFactor, LaplacianHyperbolicFactor, LaplacianSphericalFactor, geometry);

        public static double ForCauchy(
            ReadOnlySpan<double> distances, Span<double> scratch,
            BandwidthStrategy strategy = BandwidthStrategy.MadConsistencyFactor,
            SpaceGeometry geometry = SpaceGeometry.Euclidean,
            double fallback = 1.0)
            => distances.Length == 0 ? fallback : BandwidthFor(distances, scratch, strategy, CauchyFactor, CauchyHyperbolicFactor, CauchySphericalFactor, geometry);

        public static MixtureBandwidth ForMixture(
            ReadOnlySpan<double> distances, Span<double> scratch,
            BandwidthStrategy strategy = BandwidthStrategy.MadConsistencyFactor,
            SpaceGeometry geometry = SpaceGeometry.Euclidean,
            double fallback = 1.0)
        {
            if (distances.Length == 0)
                return new MixtureBandwidth(fallback, fallback, fallback);

            // Compute the strategy-dependent scale once, then apply each kernel's
            // consistency factor. For LogScaleHyperbolic the back-transform
            // depends on the kernel factor, so we go through BandwidthFor three
            // times — still cheap (the dominant cost is the MAD sort).
            return new MixtureBandwidth(
                Gaussian:  BandwidthFor(distances, scratch, strategy, GaussianFactor, GaussianHyperbolicFactor, GaussianSphericalFactor, geometry),
                Cauchy:    BandwidthFor(distances, scratch, strategy, CauchyFactor, CauchyHyperbolicFactor, CauchySphericalFactor, geometry),
                Laplacian: BandwidthFor(distances, scratch, strategy, LaplacianFactor, LaplacianHyperbolicFactor, LaplacianSphericalFactor, geometry));
        }

        // Linear kernel: compact support — δ should span the neighbourhood
        // exactly. Bandwidth strategy doesn't apply (no MAD involved); always
        // uses Max.
        public static double ForLinear(
            ReadOnlySpan<double> distances, double fallback = 1.0)
            => distances.Length == 0 ? fallback : Max(distances);

        // ── Strategy-aware dispatch ──────────────────────────────────────────
        public static double BandwidthFor(
            ReadOnlySpan<double> distances, Span<double> scratch,
            BandwidthStrategy strategy, double euclideanFactor, double hyperbolicFactor, double sphericalFactor, SpaceGeometry geometry)
        {
            double factor = ResolveConsistencyFactor(geometry, euclideanFactor, hyperbolicFactor, sphericalFactor);
            return strategy switch
            {
                BandwidthStrategy.MadConsistencyFactor => Mad(distances, scratch) * factor,
                BandwidthStrategy.QuantileNormalized   => MadOfQuantileClamped(distances, scratch, QuantileClampLevel) * factor,
                BandwidthStrategy.LogScaleHyperbolic   => LogScaleBandwidth(distances, scratch, factor),
                BandwidthStrategy.MeanEdgeDistance     => Mean(distances),
                _ => Mad(distances, scratch) * factor,
            };
        }

        private static double ResolveConsistencyFactor(
            SpaceGeometry geometry, double euclideanFactor, double hyperbolicFactor, double sphericalFactor)
            => geometry switch
            {
                SpaceGeometry.Hyperbolic => hyperbolicFactor,
                SpaceGeometry.Spherical  => sphericalFactor,
                // Euclidean: linear MAD→σ factor applies directly.
                _ => euclideanFactor,
            };

        /// <summary>
        /// MAD of distances after clamping values above the chosen quantile
        /// to the quantile value. Truncates the long tail that destabilises
        /// MAD for bounded metrics with degenerate near-boundary distances.
        /// </summary>
        public static double MadOfQuantileClamped(
            ReadOnlySpan<double> distances, Span<double> scratch, double quantile)
        {
            int n = distances.Length;
            if (n == 0) return 0.0;
            if (quantile < 0.0 || quantile > 1.0)
                throw new ArgumentOutOfRangeException(nameof(quantile), "Must be in [0, 1].");

            EnsureScratch(scratch, n);
            Span<double> work = scratch.Slice(0, n);
            for (int i = 0; i < n; i++) { Validate(distances[i]); work[i] = distances[i]; }
            work.Sort();

            // Quantile lookup (nearest-rank). n=1 maps to index 0 unconditionally.
            int qIdx = Math.Min(n - 1, (int)Math.Ceiling(quantile * n) - 1);
            if (qIdx < 0) qIdx = 0;
            double qValue = work[qIdx];

            // Clamp + MAD in one pass: write clamped values back into scratch,
            // compute median from the (still-sorted) prefix below qValue plus
            // qValue itself for clamped tail.
            for (int i = qIdx + 1; i < n; i++) work[i] = qValue;
            // work is non-decreasing after the clamp; re-derive sample median.
            double location = MedianOfSorted(work);

            for (int i = 0; i < n; i++)
            {
                double clamped = Math.Min(distances[i], qValue);
                work[i] = Math.Abs(clamped - location);
            }
            work.Sort();
            return MedianOfSorted(work);
        }

        /// <summary>
        /// Bandwidth derived from MAD on log-distances, back-transformed to
        /// distance units. For hyperbolic metrics whose NN distance
        /// distribution is roughly log-normal (volume of a hyperbolic ball
        /// grows exponentially with radius, so distances compress near the
        /// origin and stretch near the boundary).
        /// </summary>
        /// <remarks>
        /// Formula: <c>logD[i] = log(max(d[i], LogScaleFloor))</c>; compute
        /// <c>med = median(logD)</c> and <c>mad = mad(logD)</c>; bandwidth =
        /// <c>exp(med + madÂ·kernelFactor)</c>. Operates in pure log space with
        /// only a tiny floor for coincident points so <c>log(0)</c> does not
        /// produce -âˆž. An additive <c>d + 1</c> shift would distort both the
        /// location and spread once the median and MAD are computed in log
        /// space, especially near the origin where this hyperbolic route is
        /// intended to help.
        /// </remarks>
        public static double LogScaleBandwidth(
            ReadOnlySpan<double> distances, Span<double> scratch, double kernelFactor)
        {
            int n = distances.Length;
            if (n == 0) return 0.0;

            EnsureScratch(scratch, n);
            Span<double> work = scratch.Slice(0, n);
            for (int i = 0; i < n; i++)
            {
                Validate(distances[i]);
                work[i] = Math.Log(Math.Max(distances[i], LogScaleFloor));
            }
            work.Sort();
            double logMedian = MedianOfSorted(work);

            for (int i = 0; i < n; i++) work[i] = Math.Abs(Math.Log(Math.Max(distances[i], LogScaleFloor)) - logMedian);
            work.Sort();
            double logMad = MedianOfSorted(work);

            return Math.Exp(logMedian + logMad * kernelFactor);
        }

        // â”€â”€ Kernel-agnostic primitives â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

        public static double Mean(ReadOnlySpan<double> distances, double fallback = 1.0)
        {
            if (distances.Length == 0) return fallback;
            double sum = 0;
            int count = 0;
            for (int i = 0; i < distances.Length; i++)
            {
                if (!double.IsFinite(distances[i])) continue;
                Validate(distances[i]);
                sum += distances[i];
                count++;
            }
            return count == 0 ? fallback : sum / count;
        }

        public static double MeanSquare(ReadOnlySpan<double> distances, double fallback = 1.0)
        {
            if (distances.Length == 0) return fallback;
            double sum = 0.0;
            for (int i = 0; i < distances.Length; i++)
            {
                Validate(distances[i]);
                sum += distances[i] * distances[i];
            }

            return sum / distances.Length;
        }

        public static double Median(
            ReadOnlySpan<double> distances, Span<double> scratch, double fallback = 1.0)
        {
            if (distances.Length == 0) return fallback;
            int n = distances.Length;
            EnsureScratch(scratch, n);
            Span<double> work = scratch.Slice(0, n);
            for (int i = 0; i < n; i++) { Validate(distances[i]); work[i] = distances[i]; }
            work.Sort();
            return MedianOfSorted(work);
        }

        /// <summary>Raw MAD with caller-supplied location.</summary>
        public static double Mad(
            ReadOnlySpan<double> distances, double location, Span<double> scratch, double fallback = 1.0)
        {
            if (distances.Length == 0) return fallback;
            int n = distances.Length;
            EnsureScratch(scratch, n);
            Span<double> work = scratch.Slice(0, n);
            for (int i = 0; i < n; i++) { Validate(distances[i]); work[i] = Math.Abs(distances[i] - location); }
            work.Sort();
            return MedianOfSorted(work);
        }

        /// <summary>Raw MAD with internally computed sample median as location.</summary>
        public static double Mad(
            ReadOnlySpan<double> distances, Span<double> scratch, double fallback = 1.0)
        {
            if (distances.Length == 0) return fallback;
            int n = distances.Length;
            EnsureScratch(scratch, n);
            Span<double> work = scratch.Slice(0, n);
            for (int i = 0; i < n; i++) { Validate(distances[i]); work[i] = distances[i]; }
            work.Sort();
            double location = MedianOfSorted(work);
            for (int i = 0; i < n; i++) work[i] = Math.Abs(distances[i] - location);
            work.Sort();
            return MedianOfSorted(work);
        }

        public static double Max(ReadOnlySpan<double> distances, double fallback = 1.0)
        {
            if (distances.Length == 0) return fallback;
            double max = 0;
            for (int i = 0; i < distances.Length; i++)
            {
                Validate(distances[i]);
                if (distances[i] > max) max = distances[i];
            }
            return max;
        }

        // â”€â”€ Helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

        private static double MedianOfSorted(ReadOnlySpan<double> sorted)
            => Graphs.Primitives.Statistics.MedianOfSorted(sorted);

        private static void Validate(double d)
        {
            if (double.IsNaN(d) || double.IsInfinity(d) || d < 0)
                throw new ArgumentOutOfRangeException("distances",
                    "Distances must be finite and non-negative.");
        }

        private static void EnsureScratch(Span<double> scratch, int n)
        {
            if (scratch.Length < n)
                throw new ArgumentException(
                    "Scratch span is shorter than the distance sample.", "scratch");
        }
    }
}
