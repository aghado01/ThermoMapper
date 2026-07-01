using System;
using System.Collections.Generic;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Primitives;
using Maths.Geometry;

namespace Graphs.Pipeline.Scalers;

/// <summary>
/// Stage 5 — global-bandwidth kernel scaler. Estimates a single
/// MAD-based bandwidth over the NN-distance sample and applies the
/// configured kernel (Gaussian, Cauchy, Laplacian, Linear, or a kernel
/// mixture) to every edge to produce coupling weights.
/// </summary>
/// <remarks>
/// <para>The MAD sample is drawn from the selection's k-th neighbor
/// distance vector (<see cref="NeighborSelection.KthNeighborDistances"/>) —
/// the kNN search radius the kernel is meant to scale to.</para>
///
/// <para>An explicit bandwidth (constructor argument) bypasses
/// estimation entirely. A non-positive bandwidth (the default)
/// triggers MAD estimation via
/// <see cref="BandwidthEstimation"/>.</para>
/// </remarks>
public sealed class GlobalBandwidthScaler : IEdgeScaler
{
    private readonly KernelType        _kernel;
    private readonly double            _explicitBandwidth;
    private readonly BandwidthStrategy _strategy;
    private readonly SpaceGeometry     _geometry;
    private readonly CouplingFidelity  _fidelity;
    private readonly int?              _ambientDimension;
    private readonly MixtureWeights?   _mixtureWeights;
    private readonly MixtureBandwidth? _explicitMixtureBandwidth;
    private readonly double[][]?       _features;
    private readonly IRiemannianManifold? _manifold;
    private readonly SphericalIntrinsicMode _sphericalMode;

    /// <summary>Single-kernel scaler.</summary>
    public GlobalBandwidthScaler(
        KernelType        kernel,
        double            bandwidth = 0.0,
        BandwidthStrategy strategy  = BandwidthStrategy.MadConsistencyFactor,
        SpaceGeometry     geometry  = SpaceGeometry.Euclidean,
        CouplingFidelity  fidelity  = CouplingFidelity.GeodesicLinear,
        int?              ambientDimension = null,
        double[][]?       features = null,
        IRiemannianManifold? manifold = null,
        SphericalIntrinsicMode sphericalMode = SphericalIntrinsicMode.Auto)
    {
        _kernel = kernel;
        _explicitBandwidth = bandwidth;
        _strategy = strategy;
        _geometry = geometry;
        _fidelity = fidelity;
        _ambientDimension = ambientDimension;
        _mixtureWeights = null;
        _explicitMixtureBandwidth = null;
        _features = features;
        _manifold = manifold;
        _sphericalMode = sphericalMode;
    }

    /// <summary>Mixture-kernel scaler.</summary>
    public GlobalBandwidthScaler(
        MixtureWeights    mixtureWeights,
        MixtureBandwidth? mixtureBandwidth = null,
        BandwidthStrategy strategy         = BandwidthStrategy.MadConsistencyFactor,
        SpaceGeometry     geometry         = SpaceGeometry.Euclidean,
        CouplingFidelity  fidelity         = CouplingFidelity.GeodesicLinear,
        int?              ambientDimension = null,
        double[][]?       features = null,
        IRiemannianManifold? manifold = null,
        SphericalIntrinsicMode sphericalMode = SphericalIntrinsicMode.Auto)
    {
        _kernel = default;                       // unused on the mixture path
        _explicitBandwidth = 0.0;
        _strategy = strategy;
        _geometry = geometry;
        _fidelity = fidelity;
        _ambientDimension = ambientDimension;
        _mixtureWeights = mixtureWeights;
        _explicitMixtureBandwidth = mixtureBandwidth;
        _features = features;
        _manifold = manifold;
        _sphericalMode = sphericalMode;
    }

    public ScalerResult Scale(NeighborSelection refined, int n)
    {
        // Bandwidth is estimated from the selection's k-th neighbor sample
        // (the kNN search radius the kernel scales to). The selection carries
        // this field; the scaler stays agnostic to the pipeline stages that
        // produced it.
        return _mixtureWeights is { } weights
            ? ScaleMixture(refined, n, weights)
            : ScaleSingle(refined, n);
    }

    private ScalerResult ScaleSingle(NeighborSelection refined, int n)
    {
        ValidateFidelitySupport();

        // MeanEdgeDistance (BWD replication): a = mean over ALL retained
        // neighbor-pair distances (the papers' "average nearest-neighbor
        // distance" over the mutual-K graph), not the K-th-only sample.
        double delta = _explicitBandwidth > 0.0
            ? _explicitBandwidth
            : _strategy == BandwidthStrategy.MeanEdgeDistance
                ? EstimateMeanEdgeDistance(refined)
                : EstimateSingle(refined, refined.KthNeighborDistances);

        Func<double, double, double> kernelFn = _kernel switch
        {
            KernelType.Gaussian  => GaussianKernel.Evaluate,
            KernelType.Cauchy    => CauchyKernel.Evaluate,
            KernelType.Laplacian => LaplacianKernel.Evaluate,
            KernelType.Linear    => LinearKernel.Evaluate,
            _ => throw new NotSupportedException($"KernelType {_kernel} not implemented."),
        };

        bool isLocalParametrix = _geometry == SpaceGeometry.Spherical
            && _fidelity == CouplingFidelity.Intrinsic
            && _sphericalMode == SphericalIntrinsicMode.LocalParametrix;

        bool isGlobalSchoenberg = _geometry == SpaceGeometry.Spherical
            && _fidelity == CouplingFidelity.Intrinsic
            && _sphericalMode == SphericalIntrinsicMode.GlobalSchoenberg;

        int intrinsicDimension = 0;
        if (isLocalParametrix || isGlobalSchoenberg)
        {
            intrinsicDimension = ResolveSphericalIntrinsicDimension();
        }

        Neighbor[][] allNeighbors = refined.AllNeighbors;
        var edges = new List<Edge>(n * (allNeighbors.Length > 0 ? allNeighbors[0].Length : 0));

        // BWD replication kernel: J = (1/K̂)·exp(−d²/2a²) with K̂ = the actual
        // average number of neighbors per site (WBD1998 Table I). A uniform J
        // rescale only moves the T-axis (structure-neutral at the chosen state
        // point) but makes T comparable across data sets and enables the
        // q-only T_ps bracket. Gated to the replication bandwidth strategy so
        // production T-axes don't shift.
        double couplingScale = 1.0;
        if (_strategy == BandwidthStrategy.MeanEdgeDistance)
        {
            long totalNeighbors = 0;
            for (int i = 0; i < n; i++) totalNeighbors += allNeighbors[i].Length;
            if (totalNeighbors > 0) couplingScale = n / (double)totalNeighbors;
        }

        for (int i = 0; i < n; i++)
        {
            foreach (var nb in allNeighbors[i])
            {
                if (i > nb.Index) continue;          // store each undirected edge once

                double coupling;
                if (isGlobalSchoenberg)
                {
                    coupling = EvaluateGlobalSchoenbergNormalized(nb.Distance, delta, intrinsicDimension);
                }
                else
                {
                    if (isLocalParametrix && nb.Distance >= 0.9 * Math.PI) continue;
                    coupling = kernelFn(nb.Distance, delta) * VolumeCorrection(nb.Distance);
                }

                coupling *= couplingScale;

                if (coupling > 1e-10)
                    edges.Add(new Edge(i, nb.Index, coupling));
            }
        }

        return new ScalerResult(
            Graph:            CsrGraph.FromEdges(edges.ToArray(), n),
            SingleBandwidth:  delta,
            MixtureBandwidth: null);
    }

    private double EvaluateGlobalSchoenbergNormalized(double r, double delta, int m)
    {
        double t = 0.5 * delta * delta;
        double K_r = EvaluateGlobalSchoenberg(r, t, m);
        double K_0 = EvaluateGlobalSchoenbergZero(t, m);
        return K_0 > 1e-12 ? K_r / K_0 : 0.0;
    }

    private static double EvaluateGlobalSchoenberg(double r, double t, int m)
    {
        if (m == 1)
        {
            double sum = 1.0;
            for (int n = 1; n <= 50; n++)
            {
                sum += 2.0 * Math.Cos(n * r) * Math.Exp(-n * n * t);
            }
            return sum;
        }
        else
        {
            double lambda = (m - 1) / 2.0;
            double cosR = Math.Cos(r);
            
            double sum = 1.0; // n = 0 term
            
            double cPrev2 = 1.0;
            double cPrev1 = 2.0 * lambda * cosR;
            
            // n = 1 term:
            sum += ((1.0 + lambda) / lambda) * cPrev1 * Math.Exp(-(1.0 + 2.0 * lambda) * t);
            
            for (int n = 2; n <= 50; n++)
            {
                double cn = (2.0 * (n - 1.0 + lambda) * cosR * cPrev1 - (n - 2.0 + 2.0 * lambda) * cPrev2) / n;
                double coeff = (n + lambda) / lambda;
                sum += coeff * cn * Math.Exp(-n * (n + 2.0 * lambda) * t);
                
                cPrev2 = cPrev1;
                cPrev1 = cn;
            }
            return sum;
        }
    }

    private static double EvaluateGlobalSchoenbergZero(double t, int m)
    {
        if (m == 1)
        {
            double sum = 1.0;
            for (int n = 1; n <= 50; n++)
            {
                sum += 2.0 * Math.Exp(-n * n * t);
            }
            return sum;
        }
        else
        {
            double lambda = (m - 1) / 2.0;
            double sum = 1.0; // n = 0 term
            
            double gPrev2 = 1.0;
            double gPrev1 = 2.0 * lambda;
            
            // n = 1 term:
            sum += ((1.0 + lambda) / lambda) * gPrev1 * Math.Exp(-(1.0 + 2.0 * lambda) * t);
            
            for (int n = 2; n <= 50; n++)
            {
                double gn = (2.0 * (n - 1.0 + lambda) * gPrev1 - (n - 2.0 + 2.0 * lambda) * gPrev2) / n;
                double coeff = (n + lambda) / lambda;
                sum += coeff * gn * Math.Exp(-n * (n + 2.0 * lambda) * t);
                
                gPrev2 = gPrev1;
                gPrev1 = gn;
            }
            return sum;
        }
    }

    private ScalerResult ScaleMixture(
        NeighborSelection refined, int n, MixtureWeights weights)
    {
        ValidateFidelitySupport();

        MixtureBandwidth resolved = _explicitMixtureBandwidth
            ?? EstimateMixture(refined, refined.KthNeighborDistances);

        Neighbor[][] allNeighbors = refined.AllNeighbors;
        var edges = new List<Edge>(n * (allNeighbors.Length > 0 ? allNeighbors[0].Length : 0));

        for (int i = 0; i < n; i++)
        {
            foreach (var nb in allNeighbors[i])
            {
                if (i > nb.Index) continue;
                double coupling = MixtureKernel.Evaluate(nb.Distance, resolved, weights);
                if (coupling > 1e-10)
                    edges.Add(new Edge(i, nb.Index, coupling));
            }
        }

        return new ScalerResult(
            Graph:            CsrGraph.FromEdges(edges.ToArray(), n),
            SingleBandwidth:  null,
            MixtureBandwidth: resolved);
    }

    private double EstimateMeanEdgeDistance(NeighborSelection refined)
    {
        double sum = 0.0;
        int count = 0;
        foreach (var neighbors in refined.AllNeighbors)
        {
            foreach (var nb in neighbors)
            {
                sum += nb.Distance;
                count++;
            }
        }
        return count == 0 ? 1.0 : sum / count;
    }

    private double EstimateSingle(NeighborSelection refined, double[] sample)
    {
        if (_kernel == KernelType.Gaussian
            && _fidelity == CouplingFidelity.Intrinsic
            && _geometry == SpaceGeometry.Hyperbolic)
        {
            int ambientDimension = _ambientDimension
                ?? throw new InvalidOperationException(
                    "Intrinsic coupling on hyperbolic geometry requires the ambient dimension; build via GraphMetric.FromFeatures.");
            return BandwidthEstimation.ForIntrinsicGaussianHyperbolic(sample, ambientDimension);
        }

        if (_kernel == KernelType.Gaussian
            && _fidelity == CouplingFidelity.Intrinsic
            && _geometry == SpaceGeometry.Spherical)
        {
            int intrinsicDimension = ResolveSphericalIntrinsicDimension();
            return BandwidthEstimation.ForIntrinsicGaussianSpherical(sample, intrinsicDimension);
        }

        ReadOnlySpan<double> span = sample;
        Span<double> scratch = stackalloc double[sample.Length <= 512 ? sample.Length : 0];
        double[]? scratchArr = scratch.Length == 0 ? new double[sample.Length] : null;
        Span<double> effective = scratch.Length > 0 ? scratch : scratchArr!.AsSpan();
        return _kernel switch
        {
            KernelType.Gaussian  => BandwidthEstimation.ForGaussian(span, effective, _strategy, _geometry),
            KernelType.Laplacian => BandwidthEstimation.ForLaplacian(span, effective, _strategy, _geometry),
            KernelType.Cauchy    => BandwidthEstimation.ForCauchy(span, effective, _strategy, _geometry),
            KernelType.Linear    => BandwidthEstimation.ForLinear(span),
            _ => throw new NotSupportedException($"KernelType {_kernel} not implemented."),
        };
    }

    private MixtureBandwidth EstimateMixture(NeighborSelection refined, double[] sample)
    {
        ReadOnlySpan<double> span = sample;
        Span<double> scratch = stackalloc double[sample.Length <= 512 ? sample.Length : 0];
        double[]? scratchArr = scratch.Length == 0 ? new double[sample.Length] : null;
        Span<double> effective = scratch.Length > 0 ? scratch : scratchArr!.AsSpan();
        return BandwidthEstimation.ForMixture(span, effective, _strategy, _geometry);
    }

    private void ValidateFidelitySupport()
    {
        if (_fidelity != CouplingFidelity.Intrinsic)
            return;

        if (_mixtureWeights is not null || _kernel != KernelType.Gaussian)
        {
            throw new NotSupportedException(
                "Intrinsic coupling is implemented for the Gaussian kernel only (hyperbolic/spherical heat kernel).");
        }
    }

    // H^d Van Vleck volume correction: w_heat = w_Gaussian · (r/sinh(r))^((d-1)/2).
    // The d=3 path stays exact to the prior H^3 POC by returning r/sinh(r) directly.
    private double VolumeCorrection(double r)
    {
        if (_fidelity != CouplingFidelity.Intrinsic)
            return 1.0;

        if (r < 1e-12)
            return 1.0;

        if (_geometry == SpaceGeometry.Hyperbolic)
        {
            int d = _ambientDimension
                ?? throw new InvalidOperationException(
                    "Intrinsic coupling on hyperbolic geometry requires the ambient dimension; build via GraphMetric.FromFeatures.");

            double baseFactor = r / Math.Sinh(r);
            if (d == 3)
                return baseFactor;

            double exponent = (d - 1) / 2.0;
            return exponent <= 0.0 ? 1.0 : Math.Pow(baseFactor, exponent);
        }

        if (_geometry != SpaceGeometry.Spherical)
            return 1.0;

        int intrinsicDimension = ResolveSphericalIntrinsicDimension();
        double sphericalBaseFactor = r / Math.Sin(r);
        // Shortcut only when the Van Vleck exponent (m-1)/2 == 1, i.e. m == 3 (mirrors
        // the hyperbolic d == 3 case above). NOT m == 1, where the exponent is 0 and the
        // correction must be 1.0 — the general path below handles that.
        if (intrinsicDimension == 3)
            return sphericalBaseFactor;

        double sphericalExponent = (intrinsicDimension - 1) / 2.0;
        return sphericalExponent <= 0.0 ? 1.0 : Math.Pow(sphericalBaseFactor, sphericalExponent);
    }

    private int ResolveSphericalIntrinsicDimension()
    {
        if (_manifold is SphericalManifold sphericalManifold)
            return sphericalManifold.IntrinsicDimension;

        int ambientDimension = _ambientDimension
            ?? throw new InvalidOperationException(
                "Intrinsic coupling on spherical geometry requires the ambient dimension or a spherical manifold; build via GraphMetric.FromFeatures.");

        return ambientDimension - 1;
    }
}
