using System;
using System.Collections.Generic;
using Maths.Regression.Subspaces;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>How C per draw is computed (E8, decision 1).</summary>
public enum ActiveSubspaceRoute
{
    /// <summary>The exact route where it holds (a product law over terms that have factors), the sampled route otherwise.</summary>
    Auto,

    /// <summary>The moment service's derivative moments (<see cref="BassMoments.GradientOuterProduct"/>): exact.</summary>
    Exact,

    /// <summary>The gradient readout averaged over points drawn from the law, shared by every draw, with its Monte Carlo error.</summary>
    Sampled,
}

/// <summary>
/// The input metric C is read in (E8, decision 2). For metric coordinates z with box coordinates x = Lz + c, C_z = LᵀC_xL, so the
/// subspace moves under any non-orthogonal change of units; every readout states its metric.
/// </summary>
public enum InputMetric
{
    /// <summary>Each continuous input scaled to its box, the declared domain or else the training range, as the fit's transform does.</summary>
    Box,

    /// <summary>The inputs' own units.</summary>
    Raw,

    /// <summary>Whitened by the law's covariance, affine-equivariant: under a product law, each input scaled by its standard deviation,
    /// so C_z is the generalized problem C v = λ Σ⁻¹ v in z's coordinates. Under a product uniform law it is the box scaling up to a
    /// constant.</summary>
    Whitened,
}

/// <summary>
/// What E8's readout reads (<see cref="BassFit.ActiveSubspace"/>): the law, the metric, the subspace's dimension k, the route and the
/// draws. A declarative record, JSON-serializable.
/// </summary>
public sealed record ActiveSubspaceSpec
{
    /// <summary>The predictors' measures by name, in their own units, as <see cref="SensitivitySpec.Measures"/> takes them; a predictor
    /// left out takes the uniform law on its box, or uniform level probabilities.</summary>
    public IReadOnlyDictionary<string, InputMeasure>? Measures { get; init; }

    /// <summary>The metric C and the subspaces are read in.</summary>
    public InputMetric Metric { get; init; } = InputMetric.Box;

    /// <summary>k, the subspaces' dimension, from 1 to the number of continuous predictors.</summary>
    public int Dimension { get; init; } = 1;

    /// <summary>The route.</summary>
    public ActiveSubspaceRoute Route { get; init; } = ActiveSubspaceRoute.Auto;

    /// <summary>The sampled route's points N, drawn once from the law and shared by every draw.</summary>
    public int Points { get; init; } = 4096;

    /// <summary>The seed of the sampled route's points.</summary>
    public int Seed { get; init; } = 1;

    /// <summary>The draws, as indices into the fit's draws pooled in chain order (a basis fit's joint draws), repeats allowed; every
    /// draw when null.</summary>
    public IReadOnlyList<int>? Draws { get; init; }

    /// <summary>A functional fit's output metric W_y over its training output locations (slice 3): C = E[JᵀW_yJ], J the outputs × inputs
    /// Jacobian. The identity, F-X6's Q, when null; refused for a scalar fit.</summary>
    public OutputMetric? OutputMetric { get; init; }
}

/// <summary>
/// The output metric W_y of a functional fit's C = E[JᵀW_yJ] (E8, slice 3; F-X6's Q): a symmetric positive semidefinite m × m matrix over
/// the training output locations, saying which output units and which correlations count. The identity sums the locations' C;
/// quadrature weights on the diagonal integrate over the output grid. A basis fit's C is then Σ_kk′ (GᵀW_yG)_kk′ E[∇w_k∇w_k′ᵀ], with
/// G = diag(y_s)K the reconstruction's map.
/// </summary>
public sealed record OutputMetric(double[,] Matrix)
{
    /// <summary>W_y = diag(<paramref name="weights"/>), nonnegative.</summary>
    public static OutputMetric Diagonal(IReadOnlyList<double> weights)
    {
        ArgumentNullException.ThrowIfNull(weights);
        var w = new double[weights.Count, weights.Count];
        for (int l = 0; l < weights.Count; l++) w[l, l] = weights[l];
        return new OutputMetric(w);
    }
}

/// <summary>
/// What E8's co-active readout reads (slice 3; <see cref="BassFit.CoActiveSubspace(BassData, BassFit, BassData, CoActiveSpec?)"/>):
/// the law, the metric, the route, the draw pairs and, for a functional fit, the output each side reads. Two fits' draws pair by
/// position, <see cref="Draws"/> with <see cref="OtherDraws"/>, each pair weighted by the product of the two fits' weights, which samples
/// their posteriors' product; two outputs of one fit pair each draw with itself.
/// </summary>
public sealed record CoActiveSpec
{
    /// <summary>The predictors' measures by name, as <see cref="ActiveSubspaceSpec.Measures"/> takes them.</summary>
    public IReadOnlyDictionary<string, InputMeasure>? Measures { get; init; }

    /// <summary>The metric the matrices are read in.</summary>
    public InputMetric Metric { get; init; } = InputMetric.Box;

    /// <summary>The route.</summary>
    public ActiveSubspaceRoute Route { get; init; } = ActiveSubspaceRoute.Auto;

    /// <summary>The sampled route's points N, shared by every pair.</summary>
    public int Points { get; init; } = 4096;

    /// <summary>The seed of the sampled route's points.</summary>
    public int Seed { get; init; } = 1;

    /// <summary>The first fit's draws (or the one fit's), as <see cref="ActiveSubspaceSpec.Draws"/>; every draw when null.</summary>
    public IReadOnlyList<int>? Draws { get; init; }

    /// <summary>The second fit's draws, paired by position with <see cref="Draws"/>; every draw when null, which needs the two fits'
    /// selections to have one length.</summary>
    public IReadOnlyList<int>? OtherDraws { get; init; }

    /// <summary>f's output functional a_f, one weight per training output location, f = a_fᵀy: required where f is read from a functional
    /// fit, refused for a scalar one. For two outputs of one fit, the first output.</summary>
    public IReadOnlyList<double>? Output { get; init; }

    /// <summary>g's output functional a_g, as <see cref="Output"/>, over the second fit's training output locations, or for two outputs of
    /// one fit, the second output.</summary>
    public IReadOnlyList<double>? OtherOutput { get; init; }
}

/// <summary>
/// E8's co-active readout (slice 3): per draw pair, C_f, C_g and C_fg = E[∇f∇gᵀ] in the declared metric, and their posterior
/// (<see cref="CoActivePosterior"/>): V's co-active directions, the contributions and the concordances.
/// </summary>
/// <param name="Inputs">The continuous predictors, in the matrices' order.</param>
/// <param name="Excluded">The categorical predictors, which take no part (decision 3).</param>
/// <param name="Route">The route that ran.</param>
/// <param name="Metric">The metric.</param>
/// <param name="Map">L, as <see cref="ActiveSubspaceAnalysis.Map"/>: every matrix here is LᵀC_boxL.</param>
/// <param name="Measure">The law, in the transformed units.</param>
/// <param name="Points">The sampled route's points N; 0 on the exact route.</param>
/// <param name="Posterior">The pairs and their summaries.</param>
/// <param name="MonteCarloError">The sampled route's standard error of each entry of each pair's matrices, in the metric; null on the
/// exact route.</param>
/// <param name="DrawIndex">Each pair's draw of the first fit (or of the one fit).</param>
/// <param name="OtherDrawIndex">Each pair's draw of the second fit; null for two outputs of one fit.</param>
public sealed record CoActiveAnalysis(
    IReadOnlyList<string> Inputs, IReadOnlyList<string> Excluded, ActiveSubspaceRoute Route, InputMetric Metric, double[,] Map,
    ProductMeasure Measure, int Points, CoActivePosterior Posterior, IReadOnlyList<CoActiveError>? MonteCarloError,
    IReadOnlyList<int> DrawIndex, IReadOnlyList<int>? OtherDrawIndex);

/// <summary>The sampled route's standard error of each entry of one pair's C_f, C_g and C_fg (E8, slice 3).</summary>
public sealed record CoActiveError(double[,] Cf, double[,] Cg, double[,] Cfg);

/// <summary>
/// Decision 6's per-draw consistency checks, reported and not enforced, for each draw and continuous predictor in the readout's
/// order. Poincaré on a uniform input on [a, b]: total(v) ≤ ((b − a)²/π²)·C_vv in box units (NaN where the input's law is not
/// uniform); for a functional fit, with W_y = RᵀR, total(v) is Σ_r of the outputs Ry's totals, which each meet it, so the sum meets
/// C = E[JᵀW_yJ]'s. On a scalar Legendre fit under the product uniform law on the boxes, SG-4(b)'s total(v) ≤ ½·E[x_v(1 − x_v)(∂_v f)²]
/// as well. A ratio above one beyond rounding, or beyond three Monte Carlo errors on the sampled route, is a violation, which is a defect.
/// </summary>
/// <param name="Total">total(v), the input's unnormalized total variance, per draw (the moment route's, under W_y for a functional fit).</param>
/// <param name="Bound">The Poincaré bound per draw.</param>
/// <param name="Ratio">Total over Bound.</param>
/// <param name="JacobiRatio">Total over ½·the Jacobi measure, for a Legendre fit under the product uniform law; null otherwise.</param>
/// <param name="Violations">The ratios above one beyond their tolerance.</param>
public sealed record PoincareCheck(double[,] Total, double[,] Bound, double[,] Ratio, double[,]? JacobiRatio, int Violations);

/// <summary>
/// E8's readout (<see cref="BassFit.ActiveSubspace"/>): the active subspaces of a fit's draws under a law, in a declared metric.
/// </summary>
/// <param name="Inputs">The continuous predictors, in C's order.</param>
/// <param name="Excluded">The categorical predictors, which take no part in C (decision 3): their effect is Sobol's.</param>
/// <param name="Route">The route that ran.</param>
/// <param name="Metric">The metric C and the subspaces are in.</param>
/// <param name="Map">L, p × p: box coordinates are L times metric coordinates (plus a shift), so C here is LᵀC_boxL.</param>
/// <param name="Measure">The law, in the fit's transformed units.</param>
/// <param name="Points">The sampled route's points N; 0 on the exact route.</param>
/// <param name="Posterior">The draws' C, subspaces and summaries, in the metric.</param>
/// <param name="MonteCarloError">The sampled route's standard error of each entry of each draw's C, in the metric; null on the exact route.</param>
/// <param name="Poincare">The per-draw consistency checks.</param>
/// <param name="DrawIndex">Each draw's index among the fit's pooled draws.</param>
public sealed record ActiveSubspaceAnalysis(
    IReadOnlyList<string> Inputs, IReadOnlyList<string> Excluded, ActiveSubspaceRoute Route, InputMetric Metric, double[,] Map,
    ProductMeasure Measure, int Points, SubspacePosterior Posterior, double[][,]? MonteCarloError, PoincareCheck Poincare,
    IReadOnlyList<int> DrawIndex);
