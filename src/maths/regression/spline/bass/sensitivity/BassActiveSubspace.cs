using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Subspaces;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// E8's readouts (the active-subspace brief, slices 2 to 4): C per draw on one of two routes, read in a declared metric, then the
/// posterior's summaries and decision 6's checks; and the co-active pairs.
/// <list type="bullet">
/// <item><b>Exact</b> (decision 1): the moment service's derivative moments under the product law, for terms that have factors, which
/// every term type built today has. <see cref="ActiveSubspaceRoute.Auto"/> takes it under a product law.</item>
/// <item><b>Sampled:</b> the gradient readout (<see cref="BassFit.Gradient"/>) at N points shared by every draw, so the draws' spread is
/// the posterior's and not the Monte Carlo error's (decision 5): drawn once from the product law, or a joint sample's own, weighted
/// (slice 4). Each draw's subspace comes from the thin SVD of its gradient factor, not from C, and each entry of C carries its standard
/// error. <see cref="ActiveSubspaceRoute.Auto"/> takes it under a joint sample.</item>
/// </list>
/// A draw is read as its terms B_a and H, the coefficient of each term at each output location (one row for a scalar fit; c_a·h_a(ℓ)
/// for a direct functional fit, h_a its output-side factors; G_ℓk·c_ka for a basis fit, G = diag(y_s)K). Then C = E[JᵀW_yJ] is
/// Σ_ab M_ab E[∇B_a∇B_bᵀ] with M = HᵀW_yH, a co-active pair's C_fg is the same form with the two functionals' rows, and the totals behind
/// the Poincaré check are Σ_ab M_ab·T_ab(v) over the term pairs' totals. The metric (decision 2) is a map L from metric to box
/// coordinates, C = LᵀC_boxL: diagonal under a product law (1 for the box, 1/width for raw units, the law's standard deviation in box
/// units for whitening), and the symmetric root of the sample's covariance for whitening under a joint sample.
/// </summary>
internal static partial class BassActiveSubspaces
{
    private const int Batch = 64;                  // a scalar fit's draws read from the gradient readout at once
    private const long BatchValues = 1L << 22;     // a functional or basis fit's batch: about this many partials at once

    // What a readout shares: the inputs, the law (a product measure, or a joint sample's cases and weights), the continuous predictors,
    // the metric's map and the route.
    private sealed record Frame(BassInputs Inputs, ProductMeasure? Measure, BassInputTransform Transform, int[] Variables, string[] Names,
                                string[] Excluded, double[] Widths, double[,] L, bool Diagonal, ActiveSubspaceRoute Route, int Degree,
                                PredictionCases? SampleCases, double[]? SampleWeights)
    {
        public int P => Variables.Length;

        public double[,] Map => (double[,])L.Clone();
    }

    // The sampled route's points, as cases in the inputs' own units, with their weights, normalized.
    private sealed record PointSet(PredictionCases Cases, double[] Weights)
    {
        public int N => Weights.Length;
    }

    // A draw: its terms and the coefficient of each term at each output location.
    private sealed record Surface(BassTerm[] Terms, double[,] H);

    // ── Active subspaces ──────────────────────────────────────────────────────────────────────────────────────────────

    public static ActiveSubspaceAnalysis Of(BassFit fit, BassData data, ActiveSubspaceSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        return Read(fit, data, spec, spec.OutputPatterns).Analysis;
    }

    // The readout, with G's per-draw factors when patterns are asked (slice 5's output side and basis adequacy).
    private static (ActiveSubspaceAnalysis Analysis, double[][,]? Factors) Read(BassFit fit, BassData data, ActiveSubspaceSpec spec, int patterns)
    {
        Fingerprint(fit, data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        Frame frame = FrameOf(inputs, fit.Transform, fit.Config.Structure.Degree, spec.Measures, spec.Sample, spec.Metric, spec.Route, spec.Points, spec.Dimension);
        bool functional = data.Response is FunctionalResponse;
        int m = data.Response.Locations;
        if (!functional && spec.OutputMetric is not null) throw new ArgumentException("An output metric applies to a functional fit's outputs; this fit's response is scalar.", nameof(spec));
        Patterns(patterns, functional, m);
        double[,] w = Metric(spec.OutputMetric, m);
        double[,]? s0 = patterns > 0 ? PsdRoot(w) : null;
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        (int[] selected, double[] weights) = Draws3.Select(spec.Draws, pooled.Length, fit.DrawWeights, nameof(spec));
        Func<int, Surface> surface = functional ? DirectFunctional(fit, data, inputs, pooled) : d => Scalar(fit, pooled[d]);
        BassMoments? moments = frame.Measure is null ? null : BassMoments.On(inputs, frame.Degree, frame.Measure);

        var cMetric = new double[selected.Length][,];
        var draws = new DrawSubspace[selected.Length];
        double[][,]? factors = patterns > 0 ? new double[selected.Length][,] : null;
        // The terms, H and M = HᵀW_yH, for the exact route and the checks, which a product law allows.
        var forms = moments is null ? null : new (BassTerm[] Terms, double[,] M)[selected.Length];
        var hs = moments is null ? null : new double[selected.Length][,];
        if (forms is not null)
            for (int d = 0; d < selected.Length; d++)
            {
                Surface s = surface(selected[d]);
                forms[d] = (s.Terms, functional ? Quadratic(s.H, w) : Outer(Row(s.H, 0)));
                hs![d] = s.H;
            }
        double[][,]? error = null;
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            for (int d = 0; d < selected.Length; d++)
            {
                (BassTerm[] terms, double[,] mm) = forms![d];
                double[,] cBox = functional ? Symmetric(moments!.GradientCross(terms, terms, mm))
                                            : moments!.GradientOuterProduct(terms, Coefficients(pooled[selected[d]]));
                cMetric[d] = Congruence(cBox, frame.L);
                draws[d] = DrawSubspace.FromC(cMetric[d], spec.Dimension);
                if (factors is not null) factors[d] = OutputFactor(moments, frame, terms, hs![d], s0!);
            }
        }
        else
        {
            error = new double[selected.Length][,];
            double[,] root = Root(w);
            PointSet points = PointsOf(frame, spec.Points, spec.Seed);
            int batch = functional ? (int)Math.Max(1, BatchValues / ((long)points.N * m * frame.P)) : Batch;
            for (int start = 0; start < selected.Length; start += batch)
            {
                int[] chunk = selected[start..Math.Min(selected.Length, start + batch)];
                BassGradient g = fit.Gradient(data, points.Cases, chunk);
                for (int b = 0; b < chunk.Length; b++)
                {
                    int bb = b;
                    (draws[start + b], cMetric[start + b], error[start + b], double[,]? gf) =
                        SampledDraw(frame, points, root, m, (r, l, i) => g.Partials[i][bb, r, l], spec.Dimension, s0);
                    if (factors is not null) factors[start + b] = gf!;
                }
            }
        }
        var posterior = SubspacePosterior.Of(draws, weights);
        PoincareCheck? poincare = null;
        if (forms is not null)
        {
            double[,]? jacobi = !functional && LegendreOnBoxes(fit.Terms, pooled, frame, selected) ? Jacobi(fit, data, frame, selected) : null;
            poincare = Poincare(frame, moments!, forms, cMetric, error, jacobi);
        }
        (MatrixPoincareCheck? matrix, MatrixPoincareStatus status) = MatrixCheckIf(spec.MatrixPoincare, frame, moments, forms, hs, w);
        var analysis = new ActiveSubspaceAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map,
                                                  frame.Measure, PointCount(frame, spec.Points), posterior, error, poincare, Array.AsReadOnly(selected),
                                                  spec.Rank is null ? null : RankPosterior.Of(draws, weights, spec.Rank),
                                                  factors is null ? null : OutputPosterior.Of(factors, weights, patterns), matrix, status);
        return (analysis, factors);
    }

    public static ActiveSubspaceAnalysis Of(BassBasisFit fit, ActiveSubspaceSpec spec)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(spec);
        (BassInputs inputs, BassFit first, int degree) = BasisInputs(fit);
        Frame frame = FrameOf(inputs, first.Transform, degree, spec.Measures, spec.Sample, spec.Metric, spec.Route, spec.Points, spec.Dimension);
        int m = fit.Decomposition.Locations;
        Patterns(spec.OutputPatterns, true, m);
        double[,] w = Metric(spec.OutputMetric, m);
        double[,]? s0 = spec.OutputPatterns > 0 ? PsdRoot(w) : null;
        (int[] selected, Func<int, Surface> surface) = BasisSurfaces(fit, spec.Draws, nameof(spec));
        BassMoments? moments = frame.Measure is null ? null : BassMoments.On(inputs, degree, frame.Measure);
        var cMetric = new double[selected.Length][,];
        var draws = new DrawSubspace[selected.Length];
        double[][,]? factors = spec.OutputPatterns > 0 ? new double[selected.Length][,] : null;
        var forms = moments is null ? null : new (BassTerm[] Terms, double[,] M)[selected.Length];
        var hs = moments is null ? null : new double[selected.Length][,];
        if (forms is not null)
            for (int d = 0; d < selected.Length; d++)
            {
                Surface s = surface(selected[d]);
                forms[d] = (s.Terms, Quadratic(s.H, w));
                hs![d] = s.H;
            }
        double[][,]? error = null;
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            for (int d = 0; d < selected.Length; d++)
            {
                cMetric[d] = Congruence(Symmetric(moments!.GradientCross(forms![d].Terms, forms[d].Terms, forms[d].M)), frame.L);
                draws[d] = DrawSubspace.FromC(cMetric[d], spec.Dimension);
                if (factors is not null) factors[d] = OutputFactor(moments, frame, forms[d].Terms, hs![d], s0!);
            }
        }
        else
        {
            // J = G·J_w with J_w the components' gradients, so JᵀW_yJ = J_wᵀ(GᵀW_yG)J_w: the output metric's root in component space.
            error = new double[selected.Length][,];
            double[,] root = Root(ComponentMetric(fit, w));
            int r = fit.Components.Count;
            double[,]? outputMap = s0 is null ? null : Multiply(s0, ReconstructionMap(fit));   // S·G, m × K
            PointSet points = PointsOf(frame, spec.Points, spec.Seed);
            int batch = (int)Math.Max(1, BatchValues / ((long)points.N * r * frame.P));
            for (int start = 0; start < selected.Length; start += batch)
            {
                int[] chunk = selected[start..Math.Min(selected.Length, start + batch)];
                BassGradient[] g = ComponentGradients(fit, points, chunk);
                for (int b = 0; b < chunk.Length; b++)
                {
                    int bb = b;
                    (draws[start + b], cMetric[start + b], error[start + b], double[,]? gf) =
                        SampledDraw(frame, points, root, r, (q, k, i) => g[k].Partials[i][bb, q, 0], spec.Dimension, outputMap);
                    if (factors is not null) factors[start + b] = gf!;
                }
            }
        }
        double[] weights = [.. selected.Select(_ => 1.0)];
        var posterior = SubspacePosterior.Of(draws, weights);
        PoincareCheck? poincare = forms is null ? null : Poincare(frame, moments!, forms, cMetric, error, null);
        (MatrixPoincareCheck? matrix, MatrixPoincareStatus status) = MatrixCheckIf(spec.MatrixPoincare, frame, moments, forms, hs, w);
        return new ActiveSubspaceAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map,
                                          frame.Measure, PointCount(frame, spec.Points), posterior, error, poincare, Array.AsReadOnly(selected),
                                          spec.Rank is null ? null : RankPosterior.Of(draws, weights, spec.Rank),
                                          factors is null ? null : OutputPosterior.Of(factors, weights, spec.OutputPatterns), matrix, status);
    }

    // Basis adequacy (slice 5): the direct fit's G against the basis's span, both under the spec's output metric.
    public static BasisAdequacy Adequacy(BassBasisFit fit, BassFit direct, BassData directData, ActiveSubspaceSpec spec)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(direct);
        ArgumentNullException.ThrowIfNull(spec);
        Fingerprint(direct, directData);
        if (directData.Response is not FunctionalResponse)
            throw new ArgumentException("Adequacy reads a direct-mode functional fit beside the basis fit; this fit's response is scalar.", nameof(directData));
        int m = fit.Decomposition.Locations;
        if (directData.Response.Locations != m)
            throw new ArgumentException($"The direct fit has {directData.Response.Locations} output locations and the basis {m}: adequacy needs the same cases.", nameof(directData));
        (BassInputs basisInputs, BassFit first, _) = BasisInputs(fit);
        SameInputs(first.Transform, basisInputs, direct.Transform, BassInputs.Prepare(directData, direct.Transform));
        int k = Math.Max(1, spec.OutputPatterns);
        (ActiveSubspaceAnalysis analysis, double[][,]? factors) = Read(direct, directData, spec, k);
        double[,] span = Multiply(PsdRoot(Metric(spec.OutputMetric, m)), ReconstructionMap(fit));
        return BasisAdequacy.Of(factors!, analysis.Posterior.Weights, span, k);
    }

    // One draw on the sampled route. At each point the outputs' Jacobian in box units (∂/∂x = width·∂/∂z) is mapped by the output
    // metric's root, over its rows, and into the input metric by L; each mapped row over √w joins the draw's factor, and its Gram
    // matrix C's weighted running sums. jacobian(r, l, i) is output (or component) l's partial in input i, in the input's own units.
    // With an output map O (slice 5), G's factor stacks each point's O·J·L·√w side by side, m_O × Np, so G = FFᵀ is never formed.
    private static (DrawSubspace Draw, double[,] C, double[,] Error, double[,]? G) SampledDraw(Frame frame, PointSet points, double[,] root, int m,
                                                                                            Func<int, int, int, double> jacobian, int dimension,
                                                                                            double[,]? outputMap = null)
    {
        int p = frame.P, rows = root.GetLength(0), n = points.N;
        var box = new double[m, p];
        var mapped = new double[p];
        var metric = new double[rows, p];
        var factor = new double[n * rows, p];
        var sums = new WeightedSums(p, symmetric: true);
        int outputs = outputMap?.GetLength(0) ?? 0;
        double[,]? g = outputMap is null ? null : new double[outputs, n * p];
        var outputRow = new double[1, p];
        for (int r = 0; r < n; r++)
        {
            double weight = points.Weights[r];
            if (weight == 0.0) continue;
            for (int l = 0; l < m; l++)
                for (int i = 0; i < p; i++) box[l, i] = frame.Widths[i] * jacobian(r, l, i);
            double rootWeight = Math.Sqrt(weight);
            for (int q = 0; q < rows; q++)
            {
                for (int i = 0; i < p; i++)
                {
                    double v = 0.0;
                    for (int l = 0; l < m; l++)
                        if (root[q, l] != 0.0) v += root[q, l] * box[l, i];
                    mapped[i] = v;
                }
                ToMetric(frame, mapped, metric, q);
                for (int i = 0; i < p; i++) factor[r * rows + q, i] = metric[q, i] * rootWeight;
            }
            for (int i = 0; i < p; i++)
                for (int j = i; j < p; j++)
                {
                    double v = 0.0;
                    for (int q = 0; q < rows; q++) v += metric[q, i] * metric[q, j];
                    sums.Add(i, j, v, weight);
                }
            sums.AddWeight(weight);
            for (int o = 0; o < outputs; o++)
            {
                for (int i = 0; i < p; i++)
                {
                    double v = 0.0;
                    for (int l = 0; l < m; l++)
                        if (outputMap![o, l] != 0.0) v += outputMap[o, l] * box[l, i];
                    mapped[i] = v;
                }
                ToMetric(frame, mapped, outputRow, 0);
                for (int j = 0; j < p; j++) g![o, r * p + j] = outputRow[0, j] * rootWeight;
            }
        }
        (double[,] c, double[,] e) = sums.Result();
        return (DrawSubspace.FromFactor(factor, dimension), c, e, g);
    }

    // G's factor on the exact route: G = E[J W_x Jᵀ] = S·H·K·Hᵀ·S with K_ab = Σ_ij (LLᵀ)_ij E[∂ᵢB_a ∂ⱼB_b], so F = S·H·F_K with K = F_KF_Kᵀ
    // from K's eigendecomposition, its positive part.
    private static double[,] OutputFactor(BassMoments moments, Frame frame, BassTerm[] terms, double[,] h, double[,] s)
    {
        int p = frame.P, t = terms.Length;
        var wx = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                for (int k = 0; k < p; k++) wx[i, j] += frame.L[i, k] * frame.L[j, k];
        double[,] gram = moments.GradientGram(terms, wx);
        EigenResult e = ScaledEigen.Of(gram);
        double top = Math.Max(e.Eigenvalues.Max(), 0.0);
        int[] kept = [.. Enumerable.Range(0, t).Where(q => e.Eigenvalues[q] > 1e-14 * top)];
        var fk = new double[t, Math.Max(1, kept.Length)];
        for (int c = 0; c < kept.Length; c++)
        {
            double root = Math.Sqrt(e.Eigenvalues[kept[c]]);
            for (int a = 0; a < t; a++) fk[a, c] = root * e.Eigenvectors[kept[c]][a];
        }
        return Multiply(s, Multiply(h, fk));
    }

    // The matrix check when asked and applicable (a product law, which gives the moment route's forms, of uniform laws over continuous
    // inputs alone); otherwise its status says why it was not read.
    private static (MatrixPoincareCheck? Check, MatrixPoincareStatus Status) MatrixCheckIf(bool asked, Frame frame, BassMoments? moments,
                                                                                          (BassTerm[] Terms, double[,] M)[]? forms, double[][,]? hs, double[,] w)
    {
        if (!asked) return (null, MatrixPoincareStatus.NotAsked);
        if (forms is null) return (null, MatrixPoincareStatus.NotApplicable);
        MatrixPoincareCheck? check = MatrixCheck(frame, moments!, forms, hs!, Root(w));
        return (check, check is null ? MatrixPoincareStatus.NotApplicable : MatrixPoincareStatus.Computed);
    }

    // The matrix Poincaré check per draw under a product law of uniform laws over continuous inputs alone: Cov(Rf) = (RH)Σ_B(RH)ᵀ against
    // E[RJ D JᵀRᵀ] = (RH)K_D(RH)ᵀ, Σ_B the terms' covariances, K_D their gradients' Gram under D = diag((b − a)²/π²) in box units.
    private static MatrixPoincareCheck? MatrixCheck(Frame frame, BassMoments moments, (BassTerm[] Terms, double[,] M)[] forms, double[][,] hs, double[,] root)
    {
        ProductMeasure measure = frame.Measure!;
        if (frame.Inputs.Cases.CategoricalCount > 0) return null;
        int p = frame.P;
        var d = new double[p, p];
        for (int i = 0; i < p; i++)
        {
            int v = frame.Variables[i];
            if (measure.Declared[v] is not UniformMeasure) return null;
            ContinuousMeasure law = measure.Continuous(v)!;
            d[i, i] = (law.Upper - law.Lower) * (law.Upper - law.Lower) / (Math.PI * Math.PI);
        }
        var ratios = new double[forms.Length];
        int violations = 0;
        for (int s = 0; s < forms.Length; s++)
        {
            BassTerm[] terms = forms[s].Terms;
            double[,] rh = Multiply(root, hs[s]);
            double[,] covariance = Sandwich(rh, moments.CovarianceGram(terms)), bound = Sandwich(rh, moments.GradientGram(terms, d));
            ratios[s] = MatrixPoincare.Ratio(covariance, bound);
            if (!(ratios[s] <= 1.0 + 1e-9)) violations++;
        }
        return new MatrixPoincareCheck(Array.AsReadOnly(ratios), violations);
    }

    // A·M·Aᵀ.
    private static double[,] Sandwich(double[,] a, double[,] m)
    {
        double[,] am = Multiply(a, m);
        int r = a.GetLength(0), t = a.GetLength(1);
        var s = new double[r, r];
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++)
            {
                double v = 0.0;
                for (int k = 0; k < t; k++) v += am[i, k] * a[j, k];
                s[i, j] = v;
            }
        for (int i = 0; i < r; i++)
            for (int j = i + 1; j < r; j++) s[i, j] = s[j, i] = 0.5 * (s[i, j] + s[j, i]);
        return s;
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), inner = a.GetLength(1), q = b.GetLength(1);
        if (b.GetLength(0) != inner) throw new ArgumentException("The shapes do not chain.");
        var c = new double[n, q];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < inner; k++)
            {
                double x = a[i, k];
                if (x == 0.0) continue;
                for (int j = 0; j < q; j++) c[i, j] += x * b[k, j];
            }
        return c;
    }

    // The reconstruction's map G = diag(y_s)K, m × K.
    private static double[,] ReconstructionMap(BassBasisFit fit)
    {
        LinearReconstruction map = fit.Decomposition.Map;
        int m = map.Locations, r = fit.Components.Count;
        var g = new double[m, r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++) g[l, k] = map.Scales[l] * map.Basis(l, k);
        return g;
    }

    // W_y's symmetric root, positive semidefinite: S = Σ √λ_q v_qv_qᵀ, so the output side keeps the locations' coordinates.
    private static double[,] PsdRoot(double[,] w)
    {
        int m = w.GetLength(0);
        EigenResult e = ScaledEigen.Of(w);
        var s = new double[m, m];
        for (int q = 0; q < m; q++)
        {
            double root = Math.Sqrt(Math.Max(0.0, e.Eigenvalues[q]));
            if (root == 0.0) continue;
            double[] v = e.Eigenvectors[q];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) s[i, j] += root * v[i] * v[j];
        }
        for (int i = 0; i < m; i++)
            for (int j = i + 1; j < m; j++) s[i, j] = s[j, i] = 0.5 * (s[i, j] + s[j, i]);
        return s;
    }

    // Output patterns: 0 reads none; otherwise a functional fit's, from 1 to its locations.
    private static void Patterns(int patterns, bool functional, int m)
    {
        if (patterns == 0) return;
        if (!functional) throw new ArgumentException("Output patterns read a functional fit's outputs; this fit's response is scalar.", nameof(patterns));
        if (patterns < 0 || patterns > m) throw new ArgumentOutOfRangeException(nameof(patterns), $"Output patterns run from 1 to the {m} output locations.");
    }

    // A row of box-unit partials into the metric: gᵀL, into row q of the target.
    private static void ToMetric(Frame frame, double[] boxRow, double[,] target, int q)
    {
        int p = frame.P;
        if (frame.Diagonal)
        {
            for (int j = 0; j < p; j++) target[q, j] = boxRow[j] * frame.L[j, j];
            return;
        }
        for (int j = 0; j < p; j++)
        {
            double v = 0.0;
            for (int i = 0; i < p; i++) v += boxRow[i] * frame.L[i, j];
            target[q, j] = v;
        }
    }

    // ── Co-active pairs ───────────────────────────────────────────────────────────────────────────────────────────────

    // Two fits on the same inputs, each read as itself (a scalar fit) or through its output functional (a direct functional fit); draws
    // pair by position, each pair weighted by the product of the two fits' weights.
    public static CoActiveAnalysis CoActive(BassFit fit, BassData data, BassFit other, BassData otherData, CoActiveSpec spec)
    {
        ArgumentNullException.ThrowIfNull(other);
        ArgumentNullException.ThrowIfNull(spec);
        Fingerprint(fit, data);
        Fingerprint(other, otherData);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform), otherInputs = BassInputs.Prepare(otherData, other.Transform);
        SameInputs(fit.Transform, inputs, other.Transform, otherInputs);
        if (other.Config.Structure.Degree != fit.Config.Structure.Degree)
            throw new ArgumentException("The two fits' hinges have different degrees; the cross moments need one.", nameof(other));
        double[]? a = OutputOf(spec.Output, data, "Output"), b = OutputOf(spec.OtherOutput, otherData, "OtherOutput");
        Frame frame = FrameOf(inputs, fit.Transform, fit.Config.Structure.Degree, spec.Measures, spec.Sample, spec.Metric, spec.Route, spec.Points, 1);
        BassDraw[] pooledF = [.. fit.Draws.SelectMany(c => c)], pooledG = [.. other.Draws.SelectMany(c => c)];
        (int[] selF, double[] wF) = Draws3.Select(spec.Draws, pooledF.Length, fit.DrawWeights, nameof(spec));
        (int[] selG, double[] wG) = Draws3.Select(spec.OtherDraws, pooledG.Length, other.DrawWeights, nameof(spec));
        if (selF.Length != selG.Length)
            throw new ArgumentException($"The two fits' draws pair by position, so their selections need one length: {selF.Length} and {selG.Length}.", nameof(spec));
        double[] weights = [.. wF.Select((x, s) => x * wG[s])];
        var pairs = new CoActivePair[selF.Length];
        CoActiveError[]? errors = null;
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            var moments = BassMoments.On(inputs, frame.Degree, frame.Measure!);
            Func<int, Surface> sf = Surfaces(fit, data, inputs, pooledF), sg = Surfaces(other, otherData, otherInputs, pooledG);
            for (int s = 0; s < pairs.Length; s++)
            {
                Surface f = sf(selF[s]), g = sg(selG[s]);
                pairs[s] = Pair(moments, frame, f.Terms, Read(a, f.H), g.Terms, Read(b, g.H));
            }
        }
        else
        {
            PointSet points = PointsOf(frame, spec.Points, spec.Seed);
            errors = new CoActiveError[pairs.Length];
            long m = Math.Max(data.Response.Locations, otherData.Response.Locations);
            int batch = (int)Math.Clamp(BatchValues / (points.N * m * frame.P), 1, Batch);
            for (int start = 0; start < pairs.Length; start += batch)
            {
                int end = Math.Min(pairs.Length, start + batch);
                BassGradient gf = fit.Gradient(data, points.Cases, selF[start..end]), gg = other.Gradient(otherData, points.Cases, selG[start..end]);
                for (int s = 0; s < end - start; s++)
                    (pairs[start + s], errors[start + s]) = SampledPair(frame, points, Directional(gf, s, a), Directional(gg, s, b));
            }
        }
        return new CoActiveAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map, frame.Measure,
                                    PointCount(frame, spec.Points), CoActivePosterior.Of(pairs, weights),
                                    errors, Array.AsReadOnly(selF), Array.AsReadOnly(selG));
    }

    // Two outputs of one functional fit, a_fᵀy and a_gᵀy; each draw pairs with itself.
    public static CoActiveAnalysis CoActive(BassFit fit, BassData data, CoActiveSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        Fingerprint(fit, data);
        if (data.Response is not FunctionalResponse)
            throw new ArgumentException("Two outputs of one fit need a functional response; this fit's is scalar: pair two fits instead.", nameof(data));
        if (spec.OtherDraws is not null) throw new ArgumentException("Two outputs of one fit pair each draw with itself; OtherDraws is for two fits.", nameof(spec));
        int m = data.Response.Locations;
        double[] a = OutputOf(spec.Output, data, "Output")!, c = OutputOf(spec.OtherOutput, data, "OtherOutput")!;
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        Frame frame = FrameOf(inputs, fit.Transform, fit.Config.Structure.Degree, spec.Measures, spec.Sample, spec.Metric, spec.Route, spec.Points, 1);
        BassDraw[] pooled = [.. fit.Draws.SelectMany(x => x)];
        (int[] selected, double[] weights) = Draws3.Select(spec.Draws, pooled.Length, fit.DrawWeights, nameof(spec));
        var pairs = new CoActivePair[selected.Length];
        CoActiveError[]? errors = null;
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            var moments = BassMoments.On(inputs, frame.Degree, frame.Measure!);
            Func<int, Surface> surface = DirectFunctional(fit, data, inputs, pooled);
            for (int s = 0; s < pairs.Length; s++)
            {
                Surface sf = surface(selected[s]);
                pairs[s] = Pair(moments, frame, sf.Terms, Contract(a, sf.H), sf.Terms, Contract(c, sf.H));
            }
        }
        else
        {
            PointSet points = PointsOf(frame, spec.Points, spec.Seed);
            errors = new CoActiveError[pairs.Length];
            int batch = (int)Math.Max(1, BatchValues / ((long)points.N * m * frame.P));
            for (int start = 0; start < pairs.Length; start += batch)
            {
                int end = Math.Min(pairs.Length, start + batch);
                BassGradient g = fit.Gradient(data, points.Cases, selected[start..end]);
                for (int s = 0; s < end - start; s++) (pairs[start + s], errors[start + s]) = SampledPair(frame, points, Directional(g, s, a), Directional(g, s, c));
            }
        }
        return new CoActiveAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map, frame.Measure,
                                    PointCount(frame, spec.Points), CoActivePosterior.Of(pairs, weights),
                                    errors, Array.AsReadOnly(selected), null);
    }

    // Two outputs of one basis fit, over its joint draws.
    public static CoActiveAnalysis CoActive(BassBasisFit fit, CoActiveSpec spec)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(spec);
        if (spec.OtherDraws is not null) throw new ArgumentException("Two outputs of one fit pair each draw with itself; OtherDraws is for two fits.", nameof(spec));
        int m = fit.Decomposition.Locations;
        double[] a = Functional(spec.Output ?? throw new ArgumentException("A basis fit's co-active pair reads two outputs: Output is required.", nameof(spec)), m, "Output");
        double[] c = Functional(spec.OtherOutput ?? throw new ArgumentException("A basis fit's co-active pair reads two outputs: OtherOutput is required.", nameof(spec)), m, "OtherOutput");
        (BassInputs inputs, BassFit first, int degree) = BasisInputs(fit);
        Frame frame = FrameOf(inputs, first.Transform, degree, spec.Measures, spec.Sample, spec.Metric, spec.Route, spec.Points, 1);
        (int[] selected, Func<int, Surface> surface) = BasisSurfaces(fit, spec.Draws, nameof(spec));
        var pairs = new CoActivePair[selected.Length];
        CoActiveError[]? errors = null;
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            var moments = BassMoments.On(inputs, degree, frame.Measure!);
            for (int s = 0; s < pairs.Length; s++)
            {
                Surface sf = surface(selected[s]);
                pairs[s] = Pair(moments, frame, sf.Terms, Contract(a, sf.H), sf.Terms, Contract(c, sf.H));
            }
        }
        else
        {
            // aᵀJ = (Gᵀa)ᵀJ_w: each output functional as weights on the components' gradients.
            double[] alpha = ComponentFunctional(fit, a), gamma = ComponentFunctional(fit, c);
            PointSet points = PointsOf(frame, spec.Points, spec.Seed);
            errors = new CoActiveError[pairs.Length];
            int batch = (int)Math.Max(1, BatchValues / ((long)points.N * fit.Components.Count * frame.P));
            for (int start = 0; start < pairs.Length; start += batch)
            {
                int[] chunk = selected[start..Math.Min(selected.Length, start + batch)];
                BassGradient[] g = ComponentGradients(fit, points, chunk);
                for (int s = 0; s < chunk.Length; s++)
                {
                    int ss = s;
                    (pairs[start + s], errors[start + s]) = SampledPair(frame, points, (r, i) => Along(alpha, g, ss, r, i), (r, i) => Along(gamma, g, ss, r, i));
                }
            }
        }
        return new CoActiveAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map, frame.Measure,
                                    PointCount(frame, spec.Points), CoActivePosterior.Of(pairs, [.. selected.Select(_ => 1.0)]), errors,
                                    Array.AsReadOnly(selected), null);
    }

    // C_f, C_g and C_fg on the exact route, in the metric.
    private static CoActivePair Pair(BassMoments moments, Frame frame, BassTerm[] termsF, double[] u, BassTerm[] termsG, double[] v)
        => new(Congruence(moments.GradientOuterProduct(termsF, u), frame.L), Congruence(moments.GradientOuterProduct(termsG, v), frame.L),
               Congruence(moments.GradientCross(termsF, u, termsG, v), frame.L));

    // C_f, C_g and C_fg from the two gradients at the shared points (box units after the width, then the metric), weighted, with each
    // entry's standard error.
    private static (CoActivePair Pair, CoActiveError Error) SampledPair(Frame frame, PointSet points, Func<int, int, double> gradF, Func<int, int, double> gradG)
    {
        int p = frame.P;
        var sf = new WeightedSums(p, symmetric: true);
        var sg = new WeightedSums(p, symmetric: true);
        var sfg = new WeightedSums(p, symmetric: false);
        var fb = new double[p];
        var gb = new double[p];
        var f = new double[1, p];
        var g = new double[1, p];
        for (int r = 0; r < points.N; r++)
        {
            double weight = points.Weights[r];
            if (weight == 0.0) continue;
            for (int i = 0; i < p; i++)
            {
                fb[i] = frame.Widths[i] * gradF(r, i);
                gb[i] = frame.Widths[i] * gradG(r, i);
            }
            ToMetric(frame, fb, f, 0);
            ToMetric(frame, gb, g, 0);
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    sfg.Add(i, j, f[0, i] * g[0, j], weight);
                    if (j < i) continue;
                    sf.Add(i, j, f[0, i] * f[0, j], weight);
                    sg.Add(i, j, g[0, i] * g[0, j], weight);
                }
            sf.AddWeight(weight);
            sg.AddWeight(weight);
            sfg.AddWeight(weight);
        }
        (double[,] cf, double[,] ef) = sf.Result();
        (double[,] cg, double[,] eg) = sg.Result();
        (double[,] cfg, double[,] efg) = sfg.Result();
        return (new CoActivePair(cf, cg, cfg), new CoActiveError(ef, eg, efg));
    }

    // A draw's partial ∂/∂z_i at point r: the scalar output's, or aᵀ∂y/∂z_i for an output functional a.
    private static Func<int, int, double> Directional(BassGradient g, int b, double[]? a)
        => a is null ? (r, i) => g.Partials[i][b, r, 0] : (r, i) => Along(a, g.Partials[i], b, r);

    private static double Along(double[] a, double[,,] partials, int b, int r)
    {
        double s = 0.0;
        for (int l = 0; l < a.Length; l++) s += a[l] * partials[b, r, l];
        return s;
    }

    // Σ_k α_k ∂w_k/∂z_i at point r of joint draw b, the components' gradients weighted.
    private static double Along(double[] alpha, BassGradient[] components, int b, int r, int i)
    {
        double s = 0.0;
        for (int k = 0; k < alpha.Length; k++) s += alpha[k] * components[k].Partials[i][b, r, 0];
        return s;
    }

    // A draw's coefficient per term: the scalar output's row, or aᵀH.
    private static double[] Read(double[]? a, double[,] h) => a is null ? Row(h, 0) : Contract(a, h);

    // ── Draws as surfaces ─────────────────────────────────────────────────────────────────────────────────────────────

    private static Func<int, Surface> Surfaces(BassFit fit, BassData data, BassInputs inputs, BassDraw[] pooled)
        => data.Response is FunctionalResponse ? DirectFunctional(fit, data, inputs, pooled) : d => Scalar(fit, pooled[d]);

    private static Surface Scalar(BassFit fit, BassDraw draw)
    {
        BassTerm[] terms = [.. draw.Structure.TermIds.Select(id => fit.Terms[id])];
        var h = new double[1, terms.Length];
        for (int a = 0; a < terms.Length; a++) h[0, a] = draw.Coefficients[1 + a];
        return new Surface(terms, h);
    }

    private static double[] Coefficients(BassDraw draw) => draw.Coefficients[1..];

    // A direct functional fit's draws: each term's output-side factors at the training output locations, its output normalizer
    // included, times its coefficient (BassSensitivity's conditioned form, at every training location).
    private static Func<int, Surface> DirectFunctional(BassFit fit, BassData data, BassInputs inputs, BassDraw[] pooled)
    {
        var response = (FunctionalResponse)data.Response;
        int m = response.Locations, degree = fit.Config.Structure.Degree;
        IReadOnlyList<RangeScale> outputs = fit.Transform.OutputVariables;
        var columns = new double[outputs.Count][];
        for (int q = 0; q < outputs.Count; q++)
        {
            columns[q] = new double[m];
            ReadOnlySpan<double> values = response.OutputVariables[q].Values.Span;
            for (int l = 0; l < m; l++) columns[q][l] = outputs[q].Scale(values[l]);
        }
        var kernels = new DesignKernels(inputs, degree);
        var sides = new Dictionary<int, double[]>();
        double[] Side(int id)
        {
            if (sides.TryGetValue(id, out double[]? side)) return side;
            BassTerm term = fit.Terms[id];
            double scale = 1.0;
            foreach (IDesignKernel<OutputInputs> kernel in kernels.Output) scale *= kernel.Scale(term);
            side = new double[m];
            Array.Fill(side, scale);
            foreach (HingeFactor h in term.Functional)
            {
                double knot = inputs.FunctionalKnots(h.Variable)[h.Knot];
                for (int l = 0; l < m; l++) side[l] *= Hinge.Value(degree, h.Sign, knot, columns[h.Variable][l]);
            }
            return sides[id] = side;
        }
        return d =>
        {
            BassDraw draw = pooled[d];
            int[] ids = draw.Structure.TermIds;
            var h = new double[m, ids.Length];
            for (int a = 0; a < ids.Length; a++)
            {
                double[] side = Side(ids[a]);
                for (int l = 0; l < m; l++) h[l, a] = draw.Coefficients[1 + a] * side[l];
            }
            return new Surface([.. ids.Select(id => fit.Terms[id])], h);
        };
    }

    private static (BassInputs Inputs, BassFit First, int Degree) BasisInputs(BassBasisFit fit)
    {
        int degree = fit.Components[0].Config.Structure.Degree;
        for (int k = 1; k < fit.Components.Count; k++)
            if (fit.Components[k].Config.Structure.Degree != degree)
                throw new ArgumentException($"Component {k} was fitted at hinge degree {fit.Components[k].Config.Structure.Degree} and component 0 at {degree}: the cross moments need one degree.", nameof(fit));
        BassFit first = fit.Components[0];
        return (BassInputs.Prepare(fit.ComponentData(0), first.Transform), first, degree);
    }

    // A basis fit's joint draws: every component's terms, each column G_ℓk·c_ka with G = diag(y_s)K, the reconstruction's map.
    private static (int[] Selected, Func<int, Surface> Surface) BasisSurfaces(BassBasisFit fit, IReadOnlyList<int>? draws, string name)
    {
        BassJointDraws joint = fit.JointDraws;
        (int[] selected, _) = Draws3.Select(draws, joint.Count, null, name);
        int r = fit.Components.Count;
        var pooled = new List<BassDraw>[r];
        for (int k = 0; k < r; k++) pooled[k] = BassDrawFit.Pool(fit.Components[k].Draws, fit.Components[k].Terms, false, name).Pooled;
        LinearReconstruction map = fit.Decomposition.Map;
        int m = map.Locations;
        return (selected, j =>
        {
            var terms = new List<BassTerm>();
            var columns = new List<double[]>();
            for (int k = 0; k < r; k++)
            {
                BassDraw draw = pooled[k][joint.Draw(k, j)];
                int[] ids = draw.Structure.TermIds;
                for (int a = 0; a < ids.Length; a++)
                {
                    terms.Add(fit.Components[k].Terms[ids[a]]);
                    var column = new double[m];
                    for (int l = 0; l < m; l++) column[l] = map.Scales[l] * map.Basis(l, k) * draw.Coefficients[1 + a];
                    columns.Add(column);
                }
            }
            var h = new double[m, terms.Count];
            for (int a = 0; a < terms.Count; a++)
                for (int l = 0; l < m; l++) h[l, a] = columns[a][l];
            return new Surface([.. terms], h);
        });
    }

    // GᵀW_yG over the components, G = diag(y_s)K.
    private static double[,] ComponentMetric(BassBasisFit fit, double[,] w)
    {
        LinearReconstruction map = fit.Decomposition.Map;
        int m = map.Locations, r = fit.Components.Count;
        var g = new double[m, r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++) g[l, k] = map.Scales[l] * map.Basis(l, k);
        return Quadratic(g, w);
    }

    // Gᵀa: an output functional's weight on each component.
    private static double[] ComponentFunctional(BassBasisFit fit, double[] a)
    {
        LinearReconstruction map = fit.Decomposition.Map;
        int r = fit.Components.Count;
        var alpha = new double[r];
        for (int k = 0; k < r; k++)
            for (int l = 0; l < map.Locations; l++) alpha[k] += a[l] * map.Scales[l] * map.Basis(l, k);
        return alpha;
    }

    // Each component's gradient readout at the points for a chunk of joint draws, each joint draw's own draw of the component.
    private static BassGradient[] ComponentGradients(BassBasisFit fit, PointSet points, int[] chunk)
    {
        int r = fit.Components.Count;
        var g = new BassGradient[r];
        for (int k = 0; k < r; k++)
        {
            int kk = k;
            g[k] = fit.Components[k].Gradient(fit.ComponentData(k), points.Cases, [.. chunk.Select(j => fit.JointDraws.Draw(kk, j))]);
        }
        return g;
    }

    // ── The frame ─────────────────────────────────────────────────────────────────────────────────────────────────────

    private static void Fingerprint(BassFit fit, BassData data)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(data);
        string fingerprint = Bass.Fingerprint(fit.Config, data);
        if (fingerprint != fit.Identity.Fingerprint)
            throw new ArgumentException($"The fit was run on other data: its fingerprint is {fit.Identity.Fingerprint}, theirs {fingerprint}.", nameof(data));
    }

    private static Frame FrameOf(BassInputs inputs, BassInputTransform transform, int degree, IReadOnlyDictionary<string, InputMeasure>? measures,
                                 JointSample? sample, InputMetric metric, ActiveSubspaceRoute route, int points, int dimension)
    {
        if (sample is not null && measures is not null)
            throw new ArgumentException("A joint sample is the law, and Measures declare a product law: give one of them.", nameof(sample));
        if (sample is null && points < 2) throw new ArgumentOutOfRangeException(nameof(points), "The sampled route needs at least two points.");
        // Under a joint sample the product measure, the uniform default, only names the inputs.
        var measure = ProductMeasure.Resolve(transform, inputs, measures);
        int p = inputs.Cases.ContinuousCount;
        if (p == 0) throw new ArgumentException("C needs a continuous predictor: a categorical one has no derivative.", nameof(inputs));
        if (dimension < 1 || dimension > p) throw new ArgumentOutOfRangeException(nameof(dimension), $"The dimension k is from 1 to the {p} continuous predictors.");
        int[] variables = [.. Enumerable.Range(0, p).Select(measure.OfContinuous)];
        string[] names = [.. variables.Select(v => measure.Names[v])];
        string[] excluded = [.. Enumerable.Range(0, inputs.Cases.CategoricalCount).Select(c => measure.Names[measure.OfCategorical(c)])];
        double[] widths = [.. variables.Select(v => { var s = (RangeScale)transform.Predictors[v]; return s.Max - s.Min; })];
        ActiveSubspaceRoute resolved = route switch
        {
            ActiveSubspaceRoute.Auto => sample is null ? ActiveSubspaceRoute.Exact : ActiveSubspaceRoute.Sampled,
            ActiveSubspaceRoute.Exact or ActiveSubspaceRoute.Sampled => route,
            _ => throw new ArgumentOutOfRangeException(nameof(route), $"No route {route}."),
        };
        if (sample is not null && resolved == ActiveSubspaceRoute.Exact)
            throw new ArgumentException("The exact route needs a product law; a joint sample runs on the sampled route.", nameof(route));
        (PredictionCases? cases, double[]? weights, double[][]? box) = sample is null ? (null, null, null) : SampleOf(sample, transform, inputs);
        double[,] l = new double[p, p];
        bool diagonal = true;
        switch (metric)
        {
            case InputMetric.Box:
                for (int i = 0; i < p; i++) l[i, i] = 1.0;
                break;
            case InputMetric.Raw:
                for (int i = 0; i < p; i++) l[i, i] = 1.0 / widths[i];
                break;
            case InputMetric.Standardized or InputMetric.Whitened when sample is null:
                for (int i = 0; i < p; i++) l[i, i] = StandardDeviation(measure.Continuous(variables[i])!);
                break;
            case InputMetric.Standardized:
            {
                double[,] sigma = Whitening.Moments(box!, weights).Covariance;
                for (int i = 0; i < p; i++)
                {
                    if (!(sigma[i, i] > 1e-24)) throw new ArgumentException($"Standardizing needs each input to vary over the sample, its sd above 1e-12 of the box; '{names[i]}' does not.", nameof(sample));
                    l[i, i] = Math.Sqrt(sigma[i, i]);
                }
                break;
            }
            case InputMetric.Whitened:
                l = Whitening.SymmetricRoot(Whitening.Moments(box!, weights).Covariance);
                diagonal = false;
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(metric), $"No metric {metric}.");
        }
        return new Frame(inputs, sample is null ? measure : null, transform, variables, names, excluded, widths, l, diagonal, resolved, degree,
                         cases, weights);
    }

    // A joint sample's cases, its normalized weights and its continuous columns in box units. Refused: points outside the boxes (the
    // declared domains, or the training ranges), which a product law's readout never reaches either; fewer than two points of positive
    // weight; and what prediction refuses (a missing column, an unseen level).
    private static (PredictionCases Cases, double[] Weights, double[][] Box) SampleOf(JointSample sample, BassInputTransform transform, BassInputs inputs)
    {
        PredictionCases cases = sample.Points ?? throw new ArgumentException("A joint sample needs its points.", nameof(sample));
        (CaseInputs prepared, _) = cases.Prepare(transform, inputs, UnseenLevel.Refuse);
        double[][] box = [.. Enumerable.Range(0, inputs.Cases.ContinuousCount).Select(k => prepared.Continuous(k).ToArray())];
        int beyond = Enumerable.Range(0, cases.Cases).Count(r => box.Any(column => !(column[r] >= 0.0 && column[r] <= 1.0)));
        if (beyond > 0)
            throw new ArgumentException($"{beyond} of the sample's {cases.Cases} points lie outside the inputs' boxes: truncate the sample, " +
                                        "or declare wider domains when fitting.", nameof(sample));
        double[] weights = Whitening.Normalized(sample.Weights, cases.Cases, nameof(sample));
        if (!(weights.Sum(x => x * x) < 1.0)) throw new ArgumentException("A joint sample needs two points of positive weight.", nameof(sample));
        return (cases, weights, box);
    }

    private static PointSet PointsOf(Frame frame, int n, int seed)
    {
        if (frame.SampleCases is not null) return new PointSet(frame.SampleCases, frame.SampleWeights!);
        return new PointSet(Points(frame.Transform, frame.Measure!, n, seed), [.. Enumerable.Repeat(1.0 / n, n)]);
    }

    private static int PointCount(Frame frame, int n)
        => frame.Route != ActiveSubspaceRoute.Sampled ? 0 : frame.SampleCases?.Cases ?? n;

    // Two fits on the same inputs: the same predictors, boxes, training values and levels.
    private static void SameInputs(BassInputTransform a, BassInputs ia, BassInputTransform b, BassInputs ib)
    {
        const string why = "Co-active pairs need two fits on the same inputs: the same predictors, boxes and training values.";
        if (a.Predictors.Count != b.Predictors.Count) throw new ArgumentException(why);
        for (int v = 0; v < a.Predictors.Count; v++)
        {
            PredictorScale x = a.Predictors[v], y = b.Predictors[v];
            if (x.Name != y.Name) throw new ArgumentException(why);
            if (x is RangeScale rx && (y is not RangeScale ry || rx.Min != ry.Min || rx.Max != ry.Max)) throw new ArgumentException(why);
            if (x is TrainingLevels lx && (y is not TrainingLevels ly || !lx.Levels.SequenceEqual(ly.Levels))) throw new ArgumentException(why);
        }
        if (ia.Cases.ContinuousCount != ib.Cases.ContinuousCount) throw new ArgumentException(why);
        for (int k = 0; k < ia.Cases.ContinuousCount; k++)
            if (!ia.ContinuousKnots(k).SequenceEqual(ib.ContinuousKnots(k)) || !ia.Cases.Continuous(k).SequenceEqual(ib.Cases.Continuous(k)))
                throw new ArgumentException(why);
    }

    // W_y over m locations: the identity when none is declared; refused unless m × m, symmetric and positive semidefinite.
    private static double[,] Metric(OutputMetric? declared, int m)
    {
        if (declared is null)
        {
            var identity = new double[m, m];
            for (int l = 0; l < m; l++) identity[l, l] = 1.0;
            return identity;
        }
        double[,] w = declared.Matrix ?? throw new ArgumentException("The output metric needs a matrix.");
        if (w.GetLength(0) != m || w.GetLength(1) != m) throw new ArgumentException($"The output metric is {m} × {m}, one row per training output location.");
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                if (!double.IsFinite(w[i, j]) || w[i, j] != w[j, i]) throw new ArgumentException("The output metric is finite and symmetric.");
        Root(w);
        return w;
    }

    // R with W_y = RᵀR, from W_y's eigendecomposition, one row per positive eigenvalue; refused where W_y has a negative eigenvalue
    // beyond rounding.
    private static double[,] Root(double[,] w)
    {
        int m = w.GetLength(0);
        EigenResult e = ScaledEigen.Of(w);
        double top = Math.Max(e.Eigenvalues.Max(x => Math.Abs(x)), double.Epsilon);
        var rows = new List<double[]>();
        for (int q = 0; q < m; q++)
        {
            double lambda = e.Eigenvalues[q];
            if (lambda < -1e-12 * top) throw new ArgumentException("The output metric is positive semidefinite.");
            if (!(lambda > 1e-14 * top)) continue;   // rounding-level directions add nothing
            double s = Math.Sqrt(lambda);
            rows.Add([.. Enumerable.Range(0, m).Select(l => s * e.Eigenvectors[q][l])]);
        }
        if (rows.Count == 0) rows.Add(new double[m]);   // W_y = 0: one zero row, so C is 0
        var root = new double[rows.Count, m];
        for (int q = 0; q < rows.Count; q++)
            for (int l = 0; l < m; l++) root[q, l] = rows[q][l];
        return root;
    }

    // A fit's output functional: required for a functional response, one weight per training output location; refused for a scalar one.
    private static double[]? OutputOf(IReadOnlyList<double>? declared, BassData data, string name)
    {
        if (data.Response is FunctionalResponse)
            return Functional(declared ?? throw new ArgumentException($"A functional fit is read through an output functional: {name} is required.", name),
                              data.Response.Locations, name);
        if (declared is not null) throw new ArgumentException($"{name} reads a functional fit's outputs; this fit's response is scalar.", name);
        return null;
    }

    private static double[] Functional(IReadOnlyList<double> a, int m, string name)
    {
        ArgumentNullException.ThrowIfNull(a, name);
        if (a.Count != m) throw new ArgumentException($"An output functional has one weight per training output location, {m}.", name);
        foreach (double x in a)
            if (!double.IsFinite(x)) throw new ArgumentException("An output functional's weights are finite.", name);
        return [.. a];
    }

    // ── Algebra ───────────────────────────────────────────────────────────────────────────────────────────────────────

    // M = HᵀW_yH.
    private static double[,] Quadratic(double[,] h, double[,] w)
    {
        int m = h.GetLength(0), t = h.GetLength(1);
        var wh = new double[m, t];
        for (int l = 0; l < m; l++)
            for (int q = 0; q < m; q++)
            {
                double x = w[l, q];
                if (x == 0.0) continue;
                for (int a = 0; a < t; a++) wh[l, a] += x * h[q, a];
            }
        var mm = new double[t, t];
        for (int a = 0; a < t; a++)
            for (int b = 0; b < t; b++)
            {
                double s = 0.0;
                for (int l = 0; l < m; l++) s += h[l, a] * wh[l, b];
                mm[a, b] = s;
            }
        return mm;
    }

    // bbᵀ.
    private static double[,] Outer(double[] b)
    {
        var o = new double[b.Length, b.Length];
        for (int i = 0; i < b.Length; i++)
            for (int j = 0; j < b.Length; j++) o[i, j] = b[i] * b[j];
        return o;
    }

    // aᵀH, an output functional's coefficient per term.
    private static double[] Contract(double[] a, double[,] h)
    {
        int m = h.GetLength(0), t = h.GetLength(1);
        var u = new double[t];
        for (int l = 0; l < m; l++)
            for (int b = 0; b < t; b++) u[b] += a[l] * h[l, b];
        return u;
    }

    private static double[] Row(double[,] h, int l) => [.. Enumerable.Range(0, h.GetLength(1)).Select(a => h[l, a])];

    private static double[,] Symmetric(double[,] c)
    {
        int p = c.GetLength(0);
        var s = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) s[i, j] = 0.5 * (c[i, j] + c[j, i]);
        return s;
    }

    // LᵀCL. A diagonal L takes only its own entries, each product in the order Lᵢᵢ·Cᵢⱼ·Lⱼⱼ.
    private static double[,] Congruence(double[,] c, double[,] l)
    {
        int p = l.GetLength(0);
        var z = new double[p, p];
        bool diagonal = true;
        for (int i = 0; i < p && diagonal; i++)
            for (int j = 0; j < p; j++)
                if (i != j && l[i, j] != 0.0) { diagonal = false; break; }
        if (diagonal)
        {
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) z[i, j] = l[i, i] * c[i, j] * l[j, j];
            return z;
        }
        var cl = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double s = 0.0;
                for (int k = 0; k < p; k++) s += c[i, k] * l[k, j];
                cl[i, j] = s;
            }
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double s = 0.0;
                for (int k = 0; k < p; k++) s += l[k, i] * cl[k, j];
                z[i, j] = s;
            }
        return z;
    }

    // Weighted running sums of a matrix over the points, weights summing to one: the mean Σwv, and each entry's standard error from
    // Σw²(v − m)²/(1 − Σw²), which with equal weights is the sample variance over N. A symmetric matrix keeps its upper triangle.
    private sealed class WeightedSums(int p, bool symmetric)
    {
        private readonly double[,] _first = new double[p, p], _square = new double[p, p], _cross = new double[p, p];
        private double _weights;

        public void Add(int i, int j, double v, double w)
        {
            _first[i, j] += w * v;
            double w2 = w * w;
            _square[i, j] += w2 * v * v;
            _cross[i, j] += w2 * v;
        }

        public void AddWeight(double w) => _weights += w * w;

        public (double[,] Mean, double[,] Error) Result()
        {
            var mean = new double[p, p];
            var error = new double[p, p];
            double denominator = 1.0 - _weights;
            for (int i = 0; i < p; i++)
                for (int j = symmetric ? i : 0; j < p; j++)
                {
                    double m = _first[i, j];
                    double variance = Math.Max(0.0, (_square[i, j] - 2.0 * m * _cross[i, j] + m * m * _weights) / denominator);
                    mean[i, j] = m;
                    error[i, j] = Math.Sqrt(variance);
                    if (!symmetric) continue;
                    mean[j, i] = m;
                    error[j, i] = error[i, j];
                }
            return (mean, error);
        }
    }

    // The law's standard deviation in box units, from its interval moments about the support's centre.
    private static double StandardDeviation(ContinuousMeasure m)
    {
        Span<double> moments = stackalloc double[3];
        double centre = 0.5 * (m.Lower + m.Upper);
        m.Measure.IntervalMoments(Math.BitDecrement(m.Lower), m.Upper, centre, 0.0, moments);
        double mean = moments[1] / moments[0], variance = moments[2] / moments[0] - mean * mean;
        if (!(variance > 0.0)) throw new ArgumentException("Whitening needs each input's law to have a positive variance.");
        return Math.Sqrt(variance);
    }

    // N points from the product law, as cases in the inputs' own units, drawn once from the seed: each continuous input from its law in
    // box units (an empirical marginal by its atoms), mapped back through the box; each categorical one by its level probabilities.
    private static PredictionCases Points(BassInputTransform transform, ProductMeasure measure, int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var columns = new List<PredictionColumn>();
        for (int v = 0; v < measure.Count; v++)
        {
            if (measure.OutputOf(v) >= 0) continue;
            if (measure.Continuous(v) is { } m)
            {
                var s = (RangeScale)transform.Predictors[v];
                var values = new double[n];
                for (int r = 0; r < n; r++)
                {
                    double x = m.Measure is Empirical e ? e.Values[rng.NextInt(e.Count)] : ((IContinuousDistribution)m.Measure).Sample(rng);
                    values[r] = s.Min + x * (s.Max - s.Min);
                }
                columns.Add(new ContinuousColumn(measure.Names[v], values));
            }
            else
            {
                IReadOnlyList<double> probabilities = measure.Probabilities(v)!;
                var levels = ((TrainingLevels)transform.Predictors[v]).Levels;
                var picks = new string[n];
                for (int r = 0; r < n; r++)
                {
                    double u = rng.NextDouble(), cumulative = 0.0;
                    int l = 0;
                    for (; l < probabilities.Count - 1; l++)
                    {
                        cumulative += probabilities[l];
                        if (u < cumulative) break;
                    }
                    while (probabilities[l] == 0.0 && l > 0) l--;
                    picks[r] = levels[l];
                }
                columns.Add(new CategoricalColumn(measure.Names[v], picks));
            }
        }
        return new PredictionCases(columns);
    }

    // ── Decision 6 ────────────────────────────────────────────────────────────────────────────────────────────────────

    // Poincaré per draw and continuous input with a uniform law, from totals Σ_ab M_ab·T_ab(v) over the term pairs' totals, the moment
    // route's on either route; and SG-4(b) where the Jacobi measure is given. C and its errors come in the metric, which a product law
    // makes diagonal, so the box's C_vv is C_vv/L_vv².
    private static PoincareCheck Poincare(Frame frame, BassMoments moments, (BassTerm[] Terms, double[,] M)[] forms, double[][,] cMetric,
                                          double[][,]? error, double[,]? jacobi)
    {
        int draws = forms.Length, p = frame.P;
        var total = new double[draws, p];
        var bound = new double[draws, p];
        var ratio = new double[draws, p];
        var position = new Dictionary<int, int>();
        for (int i = 0; i < p; i++) position[frame.Variables[i]] = i;
        for (int d = 0; d < draws; d++)
        {
            (BassTerm[] terms, double[,] mm) = forms[d];
            for (int a = 0; a < terms.Length; a++)
                for (int b = 0; b < terms.Length; b++)
                {
                    double w = mm[a, b];
                    if (w == 0.0) continue;
                    foreach (VariableCovariance t in moments.Pair(terms[a], terms[b]).Totals)
                        if (position.TryGetValue(t.Variable, out int i)) total[d, i] += w * t.Value;
                }
        }
        double[,]? jacobiRatio = jacobi is null ? null : new double[draws, p];
        int violations = 0;
        for (int i = 0; i < p; i++)
        {
            int v = frame.Variables[i];
            bool uniform = frame.Measure!.Declared[v] is UniformMeasure;
            ContinuousMeasure m = frame.Measure.Continuous(v)!;
            double constant = uniform ? (m.Upper - m.Lower) * (m.Upper - m.Lower) / (Math.PI * Math.PI) : double.NaN;
            double toBox = 1.0 / (frame.L[i, i] * frame.L[i, i]);
            for (int d = 0; d < draws; d++)
            {
                double t = total[d, i];
                bound[d, i] = constant * cMetric[d][i, i] * toBox;
                ratio[d, i] = uniform ? (bound[d, i] > 0.0 ? t / bound[d, i] : (t > 0.0 ? double.PositiveInfinity : 0.0)) : double.NaN;
                if (uniform)
                {
                    // Rounding on the exact route; three standard errors of the bound on the sampled route.
                    double slack = error is null ? 1e-9 * Math.Max(Math.Abs(t), Math.Abs(bound[d, i]))
                                                 : 3.0 * constant * error[d][i, i] * toBox;
                    if (t > bound[d, i] + slack) violations++;
                }
                if (jacobi is not null)
                {
                    double half = 0.5 * jacobi[d, i];
                    jacobiRatio![d, i] = half > 0.0 ? t / half : (t > 0.0 ? double.PositiveInfinity : 0.0);
                    if (t > half + 1e-9 * Math.Max(Math.Abs(t), Math.Abs(half))) violations++;
                }
            }
        }
        return new PoincareCheck(total, bound, ratio, jacobiRatio, violations);
    }

    // The Jacobi measure per selected draw and continuous input, from a scalar Legendre fit's coefficients.
    private static double[,] Jacobi(BassFit fit, BassData data, Frame frame, int[] selected)
    {
        double[,] all = BassSensitivity.FromCoefficients(fit, data).Jacobi;
        var jacobi = new double[selected.Length, frame.P];
        for (int d = 0; d < selected.Length; d++)
            for (int i = 0; i < frame.P; i++) jacobi[d, i] = all[selected[d], frame.Variables[i]];
        return jacobi;
    }

    // Whether SG-4(b)'s closed form applies: every term of the draws a Legendre product, and every input's law uniform on its box.
    private static bool LegendreOnBoxes(BassTermTable table, BassDraw[] pooled, Frame frame, int[] selected)
    {
        foreach (int s in selected)
            foreach (int id in pooled[s].Structure.TermIds)
            {
                BassTerm term = table[id];
                if (term.Legendre.Count == 0 || term.Continuous.Count > 0 || term.Categorical.Count > 0 || term.Functional.Count > 0) return false;
            }
        if (frame.Measure!.Count != frame.Variables.Length) return false;
        foreach (int v in frame.Variables)
        {
            ContinuousMeasure m = frame.Measure.Continuous(v)!;
            if (frame.Measure.Declared[v] is not UniformMeasure || m.Lower != 0.0 || m.Upper != 1.0) return false;   // the box, as FromCoefficients reads it
        }
        return true;
    }
}
