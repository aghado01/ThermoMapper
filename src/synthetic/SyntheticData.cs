using System;
using System.Collections.Generic;

namespace Synthetic;

// ── Public DTOs ──────────────────────────────────────────────────────────────
// Previously nested inside the SyntheticData partial class. Now standalone so
// callers need only `using Synthetic;` — no `using static` required.

/// <summary>
/// Generators for labeled synthetic datasets used to characterize how
/// different metric / proximity / delta-estimator combinations behave on
/// controlled geometries. No generator prescribes a configuration; each
/// dataset exposes a specific structural feature (density voids, sparse
/// supports, simplex-valued points, anisotropy, parametric manifolds,
/// hierarchical scales, non-convex shapes) so experiment harnesses can
/// sweep configurations and measure outcomes directly.
/// </summary>
public sealed class SyntheticDataset
{
    public double[][] Features { get; set; } = Array.Empty<double[]>();
    public int[] Labels { get; set; } = Array.Empty<int>();
    public int ClusterCount { get; set; }
    public Dictionary<string, object> Parameters { get; set; } = new();
    public int[][]? LabelsByLevel { get; set; }
    public SyntheticDatasetMeta? Metadata { get; set; }
    public ClusterGeometry[]? ClusterGeometries { get; set; }
}

public sealed record SyntheticDatasetMeta(
    string GeneratorName,
    string GeometryClass,
    string TopologyTag,
    string HierarchyTag,
    int GTNumClusters,
    int AmbientDimensionality,
    string? LiteratureReference,
    string? SuggestedMetric = null,
    string? FutureMetric = null);

public abstract class ClusterGeometry { }

public sealed class EllipsoidGeometry : ClusterGeometry
{
    public required double[] Center { get; init; }
    public required double[,] Covariance { get; init; }
}

public sealed class ArcGeometry : ClusterGeometry
{
    public required double[][] SpineSamples { get; init; }
    public double Radius { get; init; }
    public double NoiseScale { get; init; }
}

public sealed class ManifoldGeometry : ClusterGeometry
{
    public required double[][] SpineSamples { get; init; }
    public required double[][][] TangentBases { get; init; }
}

// ── Internal shared utilities ────────────────────────────────────────────────
// Used by generators in Synthetic.Euclidean and Synthetic.Manifolds.
// Internal to the assembly; not part of the public API.

internal static class SyntheticData
{
    internal static double SampleStandardNormal(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    internal static double SampleGamma(double alpha, Random rng)
    {
        if (alpha < 1.0)
        {
            double u = rng.NextDouble();
            return SampleGamma(alpha + 1.0, rng) * Math.Pow(u, 1.0 / alpha);
        }

        double d = alpha - 1.0 / 3.0;
        double c = 1.0 / Math.Sqrt(9.0 * d);
        while (true)
        {
            double x = SampleStandardNormal(rng);
            double v = 1.0 + c * x;
            if (v <= 0) continue;
            v = v * v * v;

            double uniform = rng.NextDouble();
            double xSq = x * x;
            if (uniform < 1.0 - 0.0331 * xSq * xSq) return d * v;
            if (Math.Log(uniform) < 0.5 * xSq + d * (1.0 - v + Math.Log(v)))
                return d * v;
        }
    }

    internal static double[] SampleDirichlet(double[] alpha, Random rng)
    {
        var sample = new double[alpha.Length];
        double sum = 0;
        for (int i = 0; i < alpha.Length; i++)
        {
            sample[i] = SampleGamma(alpha[i], rng);
            sum += sample[i];
        }

        if (sum <= 0) sum = 1.0;
        for (int i = 0; i < alpha.Length; i++) sample[i] /= sum;
        return sample;
    }

    internal static void Normalize(double[] v)
    {
        double sum = 0;
        for (int i = 0; i < v.Length; i++) sum += v[i];
        if (sum <= 0) return;
        for (int i = 0; i < v.Length; i++) v[i] /= sum;
    }

    internal static double[][] PlaceCentroidsOnSphere(
        int count, int dim, double radius, Random rng)
    {
        var centroids = new double[count][];
        for (int c = 0; c < count; c++)
        {
            var centroid = new double[dim];
            double normSq = 0;
            for (int d = 0; d < dim; d++)
            {
                centroid[d] = SampleStandardNormal(rng);
                normSq += centroid[d] * centroid[d];
            }

            double norm = Math.Sqrt(normSq);
            for (int d = 0; d < dim; d++)
                centroid[d] = radius * centroid[d] / norm;
            centroids[c] = centroid;
        }

        return centroids;
    }

    internal static double[,] RandomRotationMatrix(int dim, Random rng)
    {
        if (dim == 2)
        {
            double angle = 2.0 * Math.PI * rng.NextDouble();
            double c = Math.Cos(angle), s = Math.Sin(angle);
            return new double[,] { { c, -s }, { s, c } };
        }

        var a = new double[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                a[i, j] = SampleStandardNormal(rng);
        return GramSchmidtOrthonormalize(a, dim);
    }

    internal static double[,] GramSchmidtOrthonormalize(double[,] a, int dim)
    {
        var q = new double[dim, dim];
        for (int k = 0; k < dim; k++)
        {
            var v = new double[dim];
            for (int i = 0; i < dim; i++) v[i] = a[i, k];
            for (int j = 0; j < k; j++)
            {
                double dot = 0;
                for (int i = 0; i < dim; i++) dot += v[i] * q[i, j];
                for (int i = 0; i < dim; i++) v[i] -= dot * q[i, j];
            }

            double normSq = 0;
            for (int i = 0; i < dim; i++) normSq += v[i] * v[i];
            double norm = Math.Sqrt(normSq);
            if (norm < 1e-12) norm = 1.0;
            for (int i = 0; i < dim; i++) q[i, k] = v[i] / norm;
        }

        return q;
    }

    internal static double[] MultiplyMatrixVector(double[,] m, double[] v)
    {
        int rows = m.GetLength(0), cols = m.GetLength(1);
        var result = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            double sum = 0;
            for (int j = 0; j < cols; j++) sum += m[i, j] * v[j];
            result[i] = sum;
        }
        return result;
    }

    internal static double[,] BuildCovariance(double[,] r, double[] scale)
    {
        int dim = scale.Length;
        var cov = new double[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                double sum = 0;
                for (int k = 0; k < dim; k++)
                    sum += r[i, k] * scale[k] * scale[k] * r[j, k];
                cov[i, j] = sum;
            }
        return cov;
    }

    /// <summary>Euler rotation Rz(rz) * Ry(ry) * Rx(rx), intrinsic XYZ convention.
    /// Shared between CrescentAndEllipsoid and MobiusAndEllipsoid generators.</summary>
    internal static double[,] EulerToRotation(double rx, double ry, double rz)
    {
        double cx = Math.Cos(rx), sx = Math.Sin(rx);
        double cy = Math.Cos(ry), sy = Math.Sin(ry);
        double cz = Math.Cos(rz), sz = Math.Sin(rz);

        return new double[,] {
            {  cy*cz,  cz*sx*sy - cx*sz,  cx*cz*sy + sx*sz },
            {  cy*sz,  cx*cz + sx*sy*sz,  cx*sy*sz - cz*sx },
            { -sy,     cy*sx,             cx*cy            }
        };
    }

    /// <summary>Lower-triangular Cholesky factor of a dim×dim positive-definite matrix.
    /// L s.t. L * Lᵀ = A. Shared between CrescentAndEllipsoid and MobiusAndEllipsoid.</summary>
    internal static double[,] CholeskyLower(double[,] A, int dim)
    {
        var L = new double[dim, dim];
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double sum = A[i, j];
                for (int k = 0; k < j; k++)
                    sum -= L[i, k] * L[j, k];
                L[i, j] = i == j
                    ? Math.Sqrt(Math.Max(sum, 1e-14))
                    : sum / (L[j, j] + 1e-14);
            }
        }
        return L;
    }
}
