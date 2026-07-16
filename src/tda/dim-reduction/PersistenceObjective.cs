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
        for (int i = 1; i < data.Length; i++)
            if (data[i].Length != data[0].Length)
                throw new ArgumentException("All data rows must have the same dimension.", nameof(data));
        ValidateConfig(config);

        _data = data;
        _config = config;
        _ambientDim = data[0].Length;

        // Reference: same recipe (or the ReferenceGraph override) applied to the ambient cloud.
        ReferenceBarcode = BarcodeFor(data, config.ReferenceGraph ?? config.Graph);

        _essential = config.Essential
            ?? DiagramMetrics.EssentialPolicy.FinitePenalty(Diameter(data) / 2.0);

        _covariance = config.VarianceRegularizer != 0.0 ? Covariance(data, _ambientDim) : null;
    }

    // Reject silently-degenerate recipes at construction, so every entry point — Spred,
    // DistributedSpred, direct use — inherits the same guarantees.
    private static void ValidateConfig(PersistenceObjectiveConfig config)
    {
        if (config.MaxDimension < 1)
            throw new ArgumentException(
                "MaxDimension must be >= 1 — below that the Rips complex carries no edges and there is no homology to match.",
                nameof(config));

        if (config.Dimensions is not { Length: > 0 })
            throw new ArgumentException(
                "Dimensions must weight at least one homological dimension — an empty combination makes the objective constant and the anneal optimizes nothing.",
                nameof(config));

        foreach (var (dim, weight) in config.Dimensions)
        {
            if (dim < 0 || dim >= config.MaxDimension)
                throw new ArgumentException(
                    $"Dimensions entry H{dim} needs 0 <= Dim < MaxDimension ({config.MaxDimension}): H_k requires (k+1)-simplices as fillers and the Rips complex only builds up to MaxDimension, so this term compares degenerate diagrams and goes constant.",
                    nameof(config));
            if (!double.IsFinite(weight) || weight <= 0.0)
                throw new ArgumentException(
                    $"Weight {weight} on H{dim} must be finite and > 0: a zero weight still pays the Wasserstein matching and turns an infinite essential penalty into NaN (0·∞) — drop the entry instead; a negative weight rewards barcode mismatch.",
                    nameof(config));
        }

        if (!double.IsFinite(config.WassersteinOrder) || config.WassersteinOrder < 1.0)
            throw new ArgumentException(
                "WassersteinOrder must be finite and >= 1 — DiagramMetrics requires p >= 1 (rejecting at evaluation time otherwise), and p = ∞ is Bottleneck, not Wasserstein.",
                nameof(config));

        if (config.SlicedDirections < 1)
            throw new ArgumentException(
                "SlicedDirections must be >= 1 — the sliced distance averages 1-D transports over slice directions and is undefined over zero slices.",
                nameof(config));

        if (!double.IsFinite(config.SinkhornEpsilon) || config.SinkhornEpsilon <= 0.0)
            throw new ArgumentException(
                "SinkhornEpsilon must be finite and > 0 — it is the entropic smoothing scale; 0, NaN, or ∞ degenerates the Sinkhorn kernel.",
                nameof(config));

        if (config.SinkhornMaxIters < 1)
            throw new ArgumentException(
                "SinkhornMaxIters must be >= 1 — zero iterations returns the unscaled kernel, not a transport cost.",
                nameof(config));

        if (!(config.MinPersistence >= 0.0))
            throw new ArgumentException(
                "MinPersistence must be >= 0 — it thresholds finite-bar persistence before matching (0 = no pruning).",
                nameof(config));

        if (!double.IsFinite(config.PathologyPenalty) || config.PathologyPenalty <= 0.0)
            throw new ArgumentException(
                "PathologyPenalty must be finite and > 0 — the annealer compares it as an ordinary objective value, so NaN, ∞, or a reward poisons the anneal.",
                nameof(config));
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
            total += weight * DiagramDistance(projected, ReferenceBarcode, dim);

        if (_covariance is not null)
            total += _config.VarianceRegularizer * TraceProjectedCovariance(projection, _covariance);

        return total;
    }

    // ── pipeline ─────────────────────────────────────────────────────────────

    // One dispatch point for the per-dimension diagram distance, selected by the config's
    // DiagramDistance backend (exact Hungarian vs the two screening-scale alternatives).
    private double DiagramDistance(Barcode projected, Barcode reference, int dimension) =>
        _config.DiagramDistance switch
        {
            DiagramDistanceKind.SlicedWasserstein => DiagramMetrics.SlicedWasserstein(
                projected, reference, dimension, _config.WassersteinOrder, _essential,
                _config.SlicedDirections),
            DiagramDistanceKind.SinkhornWasserstein => DiagramMetrics.SinkhornWasserstein(
                projected, reference, dimension, _config.WassersteinOrder, _essential,
                _config.SinkhornEpsilon, _config.SinkhornMaxIters),
            _ => DiagramMetrics.Wasserstein(
                projected, reference, dimension, _config.WassersteinOrder, _essential),
        };

    private Barcode BarcodeFor(double[][] features, GraphCompilerConfig recipe)
    {
        GraphMetric metric = GraphMetric.FromFeatures(features, _config.ProjectedMetric);
        GraphBuildResult built = GraphCompiler.Build(recipe, features.Length, metric);
        var filtration = RipsFiltration.GraphRips(built.Graph, _config.Filtration, _config.MaxDimension);
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
