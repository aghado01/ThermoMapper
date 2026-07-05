using System;
using System.Collections.Generic;
using Graphs;
using Maths.Geometry.DimReduction;
using TDA.Ph;

namespace TDA.DimReduction;

/// <summary>
/// The SPRED persistent-homology objective (Yu &amp; You, arXiv:2106.02096): for a candidate k×d
/// projection P, the weighted multi-order Wasserstein distance between the barcode of the projected
/// cloud P·X and the reference barcode of the ambient cloud X. Both barcodes are built by the injected
/// <see cref="GraphCompiler"/> recipe → <see cref="RipsFiltration"/> → <see cref="PersistentHomology"/>
/// (recompile-per-proposal, faithful to the paper's D(d_X) vs D(d_Y) construction).
///
/// <para><see cref="Evaluate"/> is a <see cref="SubspaceObjectiveFunction"/> by signature, so it drops
/// into <see cref="SubspaceAnnealer"/> as a method group. The reference barcode is computed once and
/// exposed via <see cref="ReferenceBarcode"/> (the §4 μ_quasi-iso measure will want it).</para>
/// </summary>
public sealed class PersistenceObjective
{
    private readonly double[][] _data;
    private readonly PersistenceObjectiveConfig _config;
    private readonly int _ambientDim;
    private readonly DiagramMetrics.EssentialPolicy _essential;
    private readonly double[][]? _covariance;   // Σ_X — materialized only when the regularizer is active

    /// <summary>The ambient reference barcode D(d_X), built once at construction.</summary>
    public Barcode ReferenceBarcode { get; }

    public PersistenceObjective(double[][] data, PersistenceObjectiveConfig config)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(config);
        if (data.Length == 0) throw new ArgumentException("Empty data", nameof(data));

        _data = data;
        _config = config;
        _ambientDim = data[0].Length;

        // Reference: same recipe (or the ReferenceGraph override) applied to the ambient cloud.
        ReferenceBarcode = BarcodeFor(data, config.ReferenceGraph ?? config.Graph);

        _essential = config.Essential
            ?? DiagramMetrics.EssentialPolicy.FinitePenalty(Diameter(data) / 2.0);

        _covariance = config.VarianceRegularizer != 0.0 ? Covariance(data, _ambientDim) : null;
    }

    /// <summary>Objective value at a candidate k×d orthonormal projection (each row a basis vector).</summary>
    public double Evaluate(double[][] projection)
    {
        Barcode projected;
        try
        {
            projected = BarcodeFor(Project(_data, projection), _config.Graph);
        }
        catch (GraphPathologyException)
        {
            return _config.PathologyPenalty;   // degenerate projected graph — reject this proposal
        }

        double total = 0.0;
        foreach (var (dim, weight) in _config.Dimensions)
            total += weight * DiagramMetrics.Wasserstein(
                projected, ReferenceBarcode, dim, _config.WassersteinOrder, _essential);

        if (_covariance is not null)
            total += _config.VarianceRegularizer * TraceProjectedCovariance(projection, _covariance);

        return total;
    }

    // ── pipeline ─────────────────────────────────────────────────────────────

    private Barcode BarcodeFor(double[][] features, GraphCompilerConfig recipe)
    {
        GraphMetric metric = GraphMetric.FromFeatures(features, _config.ProjectedMetric);
        GraphBuildResult built = GraphCompiler.Build(recipe, features.Length, metric);
        var filtration = RipsFiltration.RipsFromGraph(built.Graph, _config.Filtration, _config.MaxDimension);
        Barcode bc = PersistentHomology.Compute(filtration, _config.MaxDimension);
        return _config.MinPersistence > 0.0 ? Prune(bc, _config.MinPersistence) : bc;
    }

    // Drop finite bars below the persistence threshold; essential (infinite) bars are always kept.
    private static Barcode Prune(Barcode bc, double minPersistence)
    {
        var kept = new List<Bar>(bc.Bars.Count);
        foreach (Bar bar in bc.Bars)
            if (bar.Persistence >= minPersistence) kept.Add(bar);
        return new Barcode(kept, bc.AxisLabel);
    }

    // Y = P·X : project each ambient row x_i (length d) through the k×d frame → length-k row.
    private static double[][] Project(double[][] data, double[][] p)
    {
        int n = data.Length, k = p.Length, d = p[0].Length;
        var y = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double[] xi = data[i];
            var yi = new double[k];
            for (int r = 0; r < k; r++)
            {
                double s = 0.0;
                double[] pr = p[r];
                for (int j = 0; j < d; j++) s += pr[j] * xi[j];
                yi[r] = s;
            }
            y[i] = yi;
        }
        return y;
    }

    // ── helpers ──────────────────────────────────────────────────────────────

    // Largest pairwise ambient distance — the a_n = diam(X)/2 barcode extent the paper uses (§4).
    private static double Diameter(double[][] data)
    {
        double maxSq = 0.0;
        for (int i = 0; i < data.Length; i++)
        {
            double[] xi = data[i];
            for (int j = i + 1; j < data.Length; j++)
            {
                double[] xj = data[j];
                double s = 0.0;
                for (int c = 0; c < xi.Length; c++) { double dd = xi[c] - xj[c]; s += dd * dd; }
                if (s > maxSq) maxSq = s;
            }
        }
        return Math.Sqrt(maxSq);
    }

    private static double[][] Covariance(double[][] data, int d)
    {
        int n = data.Length;
        var mean = new double[d];
        foreach (double[] x in data) for (int j = 0; j < d; j++) mean[j] += x[j];
        for (int j = 0; j < d; j++) mean[j] /= n;

        var cov = new double[d][];
        for (int a = 0; a < d; a++) cov[a] = new double[d];
        foreach (double[] x in data)
            for (int a = 0; a < d; a++)
            {
                double xa = x[a] - mean[a];
                for (int b = 0; b < d; b++) cov[a][b] += xa * (x[b] - mean[b]);
            }
        double inv = 1.0 / n;
        for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) cov[a][b] *= inv;
        return cov;
    }

    // tr(P Σ Pᵀ) = Σ_r Σ_{a,b} P[r][a] Σ[a][b] P[r][b].
    private static double TraceProjectedCovariance(double[][] p, double[][] cov)
    {
        double trace = 0.0;
        for (int r = 0; r < p.Length; r++)
        {
            double[] pr = p[r];
            for (int a = 0; a < pr.Length; a++)
            {
                double ca = 0.0;
                double[] covA = cov[a];
                for (int b = 0; b < pr.Length; b++) ca += covA[b] * pr[b];
                trace += pr[a] * ca;
            }
        }
        return trace;
    }
}
