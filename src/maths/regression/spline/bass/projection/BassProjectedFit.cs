using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Rotate-and-refit's refit (E8, slice 7b): a BASS fit on the active variables of a <see cref="BassProjection"/>, read in the
/// original inputs' terms. 2307.11241 projects and refits so (§4.1.3, §4.3); it calls the idea a projected Bayesian MARS.
/// <list type="bullet">
/// <item><b>Prediction</b> takes cases in the original units and projects them.</item>
/// <item><b>The gradient</b> in the original inputs is the chain rule through u = Aᵀx_box: ∂f/∂x_i = Σⱼ (a_ij/width_i)·∂f/∂u_j, with
/// width_i the input's box.</item>
/// <item><b>E8's readout</b> runs on the refit in u, on the sampled route at a pushforward joint sample, since the pushforward of the
/// analysis's law is in general not a product; each draw's C then maps to the analysis's metric by C = W₁D·C_u·DW₁ᵀ, with C_u in u's
/// box units and D = diag(1/width of u's box), since LᵀA = W₁. So the mapped C lies in span(W₁), of rank at most d.</item>
/// </list>
/// The frame is a plug-in, as the PCA basis is: the refit does not carry the frame's posterior spread, which the analysis reports,
/// and the frame was read from the same data, so in-sample criteria favour the refit (slice 7c's report).
/// </summary>
public sealed class BassProjectedFit
{
    private BassProjectedFit(BassProjection projection, BassData data, BassFit fit)
    {
        Projection = projection;
        Data = data;
        Fit = fit;
    }

    /// <summary>The frame and the map.</summary>
    public BassProjection Projection { get; }

    /// <summary>The projected data the refit was run on.</summary>
    public BassData Data { get; }

    /// <summary>The refit, in the active variables.</summary>
    public BassFit Fit { get; }

    /// <summary>
    /// The refit of <paramref name="data"/>, the native fit's data or the same cases, projected by <paramref name="projection"/>
    /// (<see cref="BassProjection.Project(BassData)"/>), under <paramref name="config"/>, the native fit's by default being the
    /// caller's to pass.
    /// <para>The active variables' domains are the zonotope's bounding box, so every image of the input box lies inside, and in high
    /// dimension most of that box is far from the data. A hinge config that leaves <see cref="BassStructureSpec.HingeNormalization"/>
    /// unset is resolved to the training range (E8's slice 7d): the column-mass support and the raw ridge then read the data's
    /// scale, not the box's, and the refit's config records the choice. A stated normalization is kept as stated.</para>
    /// </summary>
    public static BassProjectedFit Of(BassProjection projection, BassData data, BassConfig config, BassRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(projection);
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(config);
        BassData projected = projection.Project(data);
        if (config.Structure is { Legendre: null, HingeNormalization: null } structure)
            config = config with { Structure = structure with { HingeNormalization = HingeNormalization.TrainingRange } };
        return new BassProjectedFit(projection, projected, Bass.Fit(config, projected, run));
    }

    /// <summary>The refit's posterior at <paramref name="cases"/>, given in the original inputs' units.</summary>
    public BassPrediction Predict(PredictionCases cases, BassPredictOptions? options = null)
        => Fit.Predict(Data, Projection.ProjectCases(cases), options);

    /// <summary>
    /// The refit's partial derivatives at <paramref name="cases"/> in the original continuous inputs, in their own units, by the chain
    /// rule through the active variables (<see cref="BassFit.Gradient"/> on the refit). <see cref="BassGradient.CaseOutsideRange"/> is
    /// the refit's: a case outside some active variable's training range.
    /// </summary>
    public BassGradient Gradient(PredictionCases cases, IReadOnlyList<int>? draws = null)
    {
        BassGradient g = Fit.Gradient(Data, Projection.ProjectCases(cases), draws);
        double[,] a = Projection.Map;
        IReadOnlyList<double> widths = Projection.Widths;
        int p = Projection.Inputs.Count, d = Projection.Dimension;
        int[] active = [.. Projection.Names.Select(name => IndexOf(g.Inputs, name))];
        double[,,] first = g.Partials[active[0]];
        int s = first.GetLength(0), n = first.GetLength(1), m = first.GetLength(2);
        var partials = new double[p][,,];
        for (int i = 0; i < p; i++)
        {
            var values = new double[s, n, m];
            for (int j = 0; j < d; j++)
            {
                double factor = a[i, j] / widths[i];
                if (factor == 0.0) continue;
                double[,,] du = g.Partials[active[j]];
                for (int b = 0; b < s; b++)
                    for (int r = 0; r < n; r++)
                        for (int l = 0; l < m; l++) values[b, r, l] += factor * du[b, r, l];
            }
            partials[i] = values;
        }
        return new BassGradient(Projection.Inputs, partials, g.Weights, g.CaseOutsideRange);
    }

    /// <summary>
    /// E8's readout of the refit (<see cref="ProjectedActiveSubspace"/>): in the active variables on the sampled route at a pushforward
    /// joint sample, and mapped to the original inputs in the analysis's metric. The spec is the refit's, with these rules: its
    /// <see cref="ActiveSubspaceSpec.Sample"/>, if given, is a joint sample of the original inputs, pushed forward; otherwise
    /// <see cref="ActiveSubspaceSpec.Points"/> points are drawn from the analysis's law with <see cref="ActiveSubspaceSpec.Seed"/>.
    /// The metric is u's box, and the route the sampled one; a declared product law (<see cref="ActiveSubspaceSpec.Measures"/>), another
    /// metric, or the exact route is refused. Its dimension k, at most d, is the mapped subspaces' too.
    /// </summary>
    public ProjectedActiveSubspace ActiveSubspace(ActiveSubspaceSpec? spec = null)
    {
        spec ??= new ActiveSubspaceSpec();
        if (spec.Measures is not null)
            throw new ArgumentException("The refit reads under the analysis's law pushed forward; a declared product law on the active variables is another law.", nameof(spec));
        if (spec.Metric != InputMetric.Box) throw new ArgumentException("The refit is read in its active variables' box; the original inputs' metric is the analysis's.", nameof(spec));
        if (spec.Route == ActiveSubspaceRoute.Exact) throw new ArgumentException("The pushforward is a joint sample, which the exact route does not take.", nameof(spec));
        if (spec.Dimension > Projection.Dimension)
            throw new ArgumentOutOfRangeException(nameof(spec), $"The refit has {Projection.Dimension} active variables, so k is at most that.");
        JointSample sample = spec.Sample is { } given ? Projection.Pushforward(given) : Projection.Pushforward(spec.Points, spec.Seed);
        ActiveSubspaceAnalysis active = Fit.ActiveSubspace(Data, spec with { Sample = sample, Route = ActiveSubspaceRoute.Sampled });

        // C = W₁D·C_u·DW₁ᵀ per draw, D = diag(1/width of u's box), u's box being its domain.
        int p = Projection.Inputs.Count, d = Projection.Dimension;
        double[,] frame = Projection.Frame;
        var b = new double[p, d];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < d; j++) b[i, j] = frame[i, j] / (Projection.Domains[j].Upper - Projection.Domains[j].Lower);
        var c = new double[active.Posterior.Draws.Count][,];
        for (int s = 0; s < c.Length; s++) c[s] = Sandwich(b, active.Posterior.Draws[s].C);
        SubspacePosterior inputs = SubspacePosterior.Of(c, active.Posterior.Weights, spec.Dimension);
        return new ProjectedActiveSubspace(active, inputs, Projection.Inputs, Projection.Metric);
    }

    /// <summary>
    /// This refit beside <paramref name="native"/>, the fit <paramref name="analysis"/> read the frame from, on its
    /// <paramref name="data"/> (slice 7c, <see cref="BassProjectionComparison"/>): each fit's size, PSIS-LOO, evidence and held-out
    /// error, the refit's subspace against the frame, the native draws' spread about it, and cross-fitting when the spec asks. Every
    /// report carries the data-reuse caveat.
    /// </summary>
    public BassProjectionComparison Compare(BassFit native, BassData data, ActiveSubspaceAnalysis analysis, ProjectionComparisonSpec? spec = null)
        => BassProjectionComparisons.Of(this, native, data, analysis, spec ?? new ProjectionComparisonSpec());

    private static int IndexOf(IReadOnlyList<string> names, string name)
    {
        for (int i = 0; i < names.Count; i++)
            if (names[i] == name) return i;
        throw new InvalidOperationException($"The refit has no input '{name}'.");
    }

    // B·M·Bᵀ, symmetrized.
    private static double[,] Sandwich(double[,] b, double[,] m)
    {
        int p = b.GetLength(0), d = b.GetLength(1);
        var bm = new double[p, d];
        for (int i = 0; i < p; i++)
            for (int k = 0; k < d; k++)
            {
                double s = 0.0;
                for (int j = 0; j < d; j++) s += b[i, j] * m[j, k];
                bm[i, k] = s;
            }
        var c = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int l = 0; l < p; l++)
            {
                double s = 0.0;
                for (int k = 0; k < d; k++) s += bm[i, k] * b[l, k];
                c[i, l] = s;
            }
        for (int i = 0; i < p; i++)
            for (int l = i + 1; l < p; l++) c[i, l] = c[l, i] = 0.5 * (c[i, l] + c[l, i]);
        return c;
    }
}

/// <summary>
/// E8's readout of a projected refit (slice 7b): <paramref name="Active"/> in the active variables, read on the sampled route at the
/// pushforward, and <paramref name="InputSpace"/>, each draw's C mapped to the original continuous inputs in the analysis's metric,
/// with its two summaries (E[C] and the projection mean). Every mapped C lies in span(W₁).
/// </summary>
/// <param name="Active">The refit's analysis in its active variables' box.</param>
/// <param name="InputSpace">The mapped draws and their summaries, in <paramref name="Metric"/>.</param>
/// <param name="Inputs">The original continuous inputs, in <paramref name="InputSpace"/>'s order.</param>
/// <param name="Metric">The analysis's metric.</param>
public sealed record ProjectedActiveSubspace(ActiveSubspaceAnalysis Active, SubspacePosterior InputSpace, IReadOnlyList<string> Inputs, InputMetric Metric);
