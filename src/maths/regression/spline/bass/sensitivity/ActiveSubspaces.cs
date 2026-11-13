using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Regression.Subspaces;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>How C per draw is computed (E8, decision 1).</summary>
public enum ActiveSubspaceRoute
{
    /// <summary>The exact route where it holds (a product law over terms that have factors), the sampled route otherwise.</summary>
    Auto,

    /// <summary>The moment service's derivative moments (<see cref="BassMoments.GradientOuterProduct"/>): exact.</summary>
    Exact,

    /// <summary>The gradient readout averaged over points shared by every draw, drawn from a product law or given as a joint sample
    /// (<see cref="JointSample"/>), with its Monte Carlo error.</summary>
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

    /// <summary>Whitened by the law's covariance Σ_ρ in box units, L = Σ_ρ^{1/2} (the symmetric root, <see cref="Whitening"/>), so
    /// C_z = LᵀCL is the generalized problem C v = λ Σ_ρ⁻¹ v: affine-equivariant. Under a product law Σ_ρ is diagonal and L scales each
    /// input by its standard deviation, as <see cref="Standardized"/> does; under a joint sample, Σ_ρ is the sample's full covariance and
    /// L mixes the inputs. Under a product uniform law it is the box scaling up to a constant.</summary>
    Whitened,

    /// <summary>Each input scaled by its standard deviation under the law, L = diag(σ) in box units: whitening's diagonal, which keeps
    /// the input axes. It is <see cref="Whitened"/> under a product law; under a joint sample it leaves the inputs' correlations in C's
    /// metric, and it is not affine-equivariant.</summary>
    Standardized,
}

/// <summary>
/// A dependent law over the inputs given as a joint sample of points (E8, slice 4): the law is the sample's own, each point weighted by
/// <paramref name="Weights"/> (nonnegative, of positive sum, normalized; equal when null). <paramref name="Points"/> are in the inputs'
/// own units, one column per predictor of the fit, continuous and categorical, as prediction takes them, and lie in each input's box
/// (its declared domain, or else its training range), as a product law's readout does. A calibration's posterior draws or an SMC run's weighted particles are such samples. Under one, the readout
/// runs the sampled route at its points, the draws sharing them, and whitening takes its full covariance.
/// </summary>
public sealed record JointSample(PredictionCases Points, IReadOnlyList<double>? Weights = null);

/// <summary>
/// What E8's readout reads (<see cref="BassFit.ActiveSubspace"/>): the law, the metric, the subspace's dimension k, the route and the
/// draws. A declarative record, JSON-serializable.
/// </summary>
public sealed record ActiveSubspaceSpec
{
    /// <summary>The predictors' measures by name, in their own units, as <see cref="SensitivitySpec.Measures"/> takes them; a predictor
    /// left out takes the uniform law on its box, or uniform level probabilities. Together they are a product law.</summary>
    public IReadOnlyDictionary<string, InputMeasure>? Measures { get; init; }

    /// <summary>A dependent law as a joint sample (slice 4), in place of <see cref="Measures"/>: the sampled route at its points, which
    /// <see cref="Points"/> and <see cref="Seed"/> then do not draw; the exact route is refused.</summary>
    public JointSample? Sample { get; init; }

    /// <summary>The metric C and the subspaces are read in.</summary>
    public InputMetric Metric { get; init; } = InputMetric.Box;

    /// <summary>k, the subspaces' dimension, from 1 to the number of continuous predictors.</summary>
    public int Dimension { get; init; } = 1;

    /// <summary>The route.</summary>
    public ActiveSubspaceRoute Route { get; init; } = ActiveSubspaceRoute.Auto;

    /// <summary>The sampled route's points N, drawn once from the product law and shared by every draw.</summary>
    public int Points { get; init; } = 4096;

    /// <summary>The seed of the sampled route's points.</summary>
    public int Seed { get; init; } = 1;

    /// <summary>The draws, as indices into the fit's draws pooled in chain order (a basis fit's joint draws), repeats allowed; every
    /// draw when null.</summary>
    public IReadOnlyList<int>? Draws { get; init; }

    /// <summary>A functional fit's output metric W_y over its training output locations (slice 3): C = E[JᵀW_yJ], J the outputs × inputs
    /// Jacobian. The identity, F-X6's Q, when null; refused for a scalar fit.</summary>
    public OutputMetric? OutputMetric { get; init; }

    /// <summary>A declared rule for each draw's own k_s (slice 5), read beside the fixed-k readout as <see cref="RankPosterior"/>. There
    /// is no default: without a rule no per-draw rank is read.</summary>
    public RankRule? Rank { get; init; }

    /// <summary>The output patterns to read (slice 5): for a functional fit, G = E[J W_x Jᵀ]'s leading patterns per draw and for E[G],
    /// in the input metric W_x = LLᵀ and the output metric's symmetric root. 0, the default, reads none; refused for a scalar fit.</summary>
    public int OutputPatterns { get; init; }

    /// <summary>Whether to read the matrix Poincaré check (slice 5) where it applies. Off by default: it about doubles the exact
    /// readout's cost per draw on scalar, functional and basis fits alike (the slice-5 cost row), so it is read when asked.</summary>
    public bool MatrixPoincare { get; init; }
}

/// <summary>Whether the matrix Poincaré check was read (E8, slice 5).</summary>
public enum MatrixPoincareStatus
{
    /// <summary>Read: <see cref="ActiveSubspaceAnalysis.MatrixPoincare"/> holds it.</summary>
    Computed,

    /// <summary>Not asked: <see cref="ActiveSubspaceSpec.MatrixPoincare"/> is off, its default.</summary>
    NotAsked,

    /// <summary>Asked, but the law is not a product over continuous inputs alone of laws with a Poincaré constant here: a uniform law,
    /// or a truncated normal (<see cref="MatrixPoincareCheck"/>).</summary>
    NotApplicable,
}

/// <summary>
/// The matrix Poincaré check per draw (E8, slice 5): the largest generalized eigenvalue of the outputs' covariance Cov_ρ(f) against
/// E[J D Jᵀ], D the inputs' Poincaré constants in box units, both under the output metric's root; at most 1 when Efron–Stein and each
/// input's Poincaré inequality hold. Read from the moment route under a product law whose every input is continuous with a known
/// constant (slice 6): (b − a)²/π² for a uniform law on [a, b]; min(σ², (b − a)²/π²) for a normal N(μ, σ²) truncated to [a, b], a
/// mixture of one normal component, since it is log-concave with potential curvature 1/σ² (Bakry–Émery) on an interval (the weighted
/// Payne–Weinberger bound), both cited from memory. A ratio above 1 + 10⁻⁹ is a violation, which is a defect. The term-pair Grams it
/// reads are cached per structure within a readout (slice 6).
/// </summary>
/// <param name="Ratios">Each draw's ratio.</param>
/// <param name="Violations">The ratios above one beyond rounding.</param>
public sealed record MatrixPoincareCheck(IReadOnlyList<double> Ratios, int Violations);

/// <summary>
/// The output metric W_y of a functional fit's C = E[JᵀW_yJ] (E8, slice 3; F-X6's Q): a symmetric positive semidefinite m × m matrix over
/// the output locations, saying which output units and which correlations count. The identity sums the locations' C; quadrature weights
/// on the diagonal integrate over the output grid. A basis fit's C is then Σ_kk′ (GᵀW_yG)_kk′ E[∇w_k∇w_k′ᵀ], with G = diag(y_s)K the
/// reconstruction's map.
/// <para>The locations are the training ones unless <see cref="Locations"/> declares others (slice 6): a direct fit's output-side factors
/// are hinges in the output variables, defined across their box, so its J reads at any location there. A basis exists only at its
/// training locations, so a basis fit refuses them.</para>
/// </summary>
public sealed record OutputMetric(double[,] Matrix)
{
    /// <summary>The locations <see cref="Matrix"/> is over, one column per output variable matched by name, in the variables' own units,
    /// as prediction's <c>OutputLocations</c> takes them, each inside the variable's box (its declared domain, else its training range);
    /// the training locations when null.</summary>
    public IReadOnlyList<ContinuousColumn>? Locations { get; init; }

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

    /// <summary>A dependent law as a joint sample, as <see cref="ActiveSubspaceSpec.Sample"/> takes it.</summary>
    public JointSample? Sample { get; init; }

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
/// <param name="Measure">The product law, in the transformed units; null under a joint sample.</param>
/// <param name="Points">The sampled route's points N; 0 on the exact route.</param>
/// <param name="Posterior">The pairs and their summaries.</param>
/// <param name="MonteCarloError">The sampled route's standard error of each entry of each pair's matrices, in the metric; null on the
/// exact route.</param>
/// <param name="DrawIndex">Each pair's draw of the first fit (or of the one fit).</param>
/// <param name="OtherDrawIndex">Each pair's draw of the second fit; null for two outputs of one fit.</param>
public sealed record CoActiveAnalysis(
    IReadOnlyList<string> Inputs, IReadOnlyList<string> Excluded, ActiveSubspaceRoute Route, InputMetric Metric, double[,] Map,
    ProductMeasure? Measure, int Points, CoActivePosterior Posterior, IReadOnlyList<CoActiveError>? MonteCarloError,
    IReadOnlyList<int> DrawIndex, IReadOnlyList<int>? OtherDrawIndex);

/// <summary>The sampled route's standard error of each entry of one pair's C_f, C_g and C_fg (E8, slice 3).</summary>
public sealed record CoActiveError(double[,] Cf, double[,] Cg, double[,] Cfg);

/// <summary>
/// Decision 6's per-draw consistency checks, reported and not enforced, for each draw and continuous predictor in the readout's
/// order. Poincaré on a uniform input on [a, b]: total(v) ≤ ((b − a)²/π²)·C_vv in box units, and since slice 6 on a normal truncated to
/// [a, b] with min(σ², (b − a)²/π²) as the constant (<see cref="MatrixPoincareCheck"/>); NaN where no constant is known here; for a functional fit, with W_y = RᵀR, total(v) is Σ_r of the outputs Ry's totals, which each meet it, so the sum meets
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
/// <param name="Map">L, p × p: box coordinates are L times metric coordinates (plus a shift), so C here is LᵀC_boxL; diagonal but for
/// whitening under a joint sample.</param>
/// <param name="Measure">The product law, in the fit's transformed units; null under a joint sample.</param>
/// <param name="Points">The sampled route's points N, drawn or the joint sample's; 0 on the exact route.</param>
/// <param name="Posterior">The draws' C, subspaces and summaries, in the metric.</param>
/// <param name="MonteCarloError">The sampled route's standard error of each entry of each draw's C, in the metric; null on the exact route.</param>
/// <param name="Poincare">The per-draw consistency checks; null under a joint sample, which has no marginal bound to check.</param>
/// <param name="DrawIndex">Each draw's index among the fit's pooled draws.</param>
/// <param name="Ranks">Each draw at its own rank by the spec's rule (slice 5); null without a rule.</param>
/// <param name="Outputs">G = E[J W_x Jᵀ]'s patterns (slice 5); null unless asked of a functional fit.</param>
/// <param name="MatrixPoincare">The matrix Poincaré check (slice 5); null unless asked and the law is a product of uniform laws over
/// continuous inputs alone.</param>
/// <param name="MatrixPoincareStatus">Whether the matrix check was read, not asked, or not applicable.</param>
/// <param name="Boxes">Each continuous predictor's box in its own units, in <paramref name="Inputs"/>' order: its declared domain, or else
/// its training range, the box units C and <paramref name="Map"/> are read in (slice 7a). A frame taken from C holds only in these
/// boxes, so a projection refuses a fit whose boxes differ.</param>
public sealed record ActiveSubspaceAnalysis(
    IReadOnlyList<string> Inputs, IReadOnlyList<string> Excluded, ActiveSubspaceRoute Route, InputMetric Metric, double[,] Map,
    ProductMeasure? Measure, int Points, SubspacePosterior Posterior, double[][,]? MonteCarloError, PoincareCheck? Poincare,
    IReadOnlyList<int> DrawIndex, RankPosterior? Ranks = null, OutputPosterior? Outputs = null, MatrixPoincareCheck? MatrixPoincare = null,
    MatrixPoincareStatus MatrixPoincareStatus = MatrixPoincareStatus.NotApplicable, IReadOnlyList<InputDomain>? Boxes = null)
{
    /// <summary>The activity scores (slice 6), in the metric: C's diagonal always, and its truncation at k, and at each draw's k_s under a
    /// rank rule (<see cref="ActivityScores"/>). Read from <see cref="Posterior"/> and <see cref="Ranks"/> when asked.</summary>
    public ActivityScores Activity => ActivityScores.Of(Posterior, Ranks);

    private readonly (double[,] Map, double[,] Inverse)? _mapInverse;

    /// <summary>L⁻¹, which is L⁻ᵀ for the diagonal and symmetric maps a readout builds (spectral-solver brief, C5). The readout sets it
    /// as it built L: the diagonal's reciprocals, or under whitening the inverse root from the covariance's own eigensystem. It is kept
    /// with the map it inverts, so a copy with another <see cref="Map"/>, or an analysis built by hand, derives it from its own map.</summary>
    internal double[,] MapInverse
    {
        get => _mapInverse is { } kept && ReferenceEquals(kept.Map, Map) ? kept.Inverse : InverseOf(Map);
        init => _mapInverse = (Map, value);
    }

    // L⁻¹ for a diagonal map, or a symmetric positive definite one as V·Λ⁻¹·Vᵀ, symmetrized.
    private static double[,] InverseOf(double[,] l)
    {
        int p = l.GetLength(0);
        bool diagonal = true;
        for (int i = 0; i < p && diagonal; i++)
            for (int j = 0; j < p; j++)
                if (i != j && l[i, j] != 0.0) { diagonal = false; break; }
        if (!diagonal) return MatrixFunctions.Symmetric(DenseEigen.DecomposeSymmetric(l, tol: 1e-14), x => 1.0 / x);
        var inverse = new double[p, p];
        for (int i = 0; i < p; i++) inverse[i, i] = 1.0 / l[i, i];
        return inverse;
    }
}
