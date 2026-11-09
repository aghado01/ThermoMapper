using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Subspaces;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// E8's readouts (the active-subspace brief, slices 2 and 3): C per draw on one of two routes, read in a declared metric, then the
/// posterior's summaries and decision 6's checks; and the co-active pairs.
/// <list type="bullet">
/// <item><b>Exact</b> (decision 1): the moment service's derivative moments under the product law, for terms that have factors, which
/// every term type built today has. <see cref="ActiveSubspaceRoute.Auto"/> takes it.</item>
/// <item><b>Sampled:</b> the gradient readout (<see cref="BassFit.Gradient"/>) at N points drawn once from the law and shared by every
/// draw, so the draws' spread is the posterior's and not the Monte Carlo error's (decision 5). Each draw's subspace comes from the thin
/// SVD of its gradient factor, not from C, and each entry of C carries its standard error.</item>
/// </list>
/// A draw is read as its terms B_a and H, the coefficient of each term at each output location (one row for a scalar fit; c_a·h_a(ℓ)
/// for a direct functional fit, h_a its output-side factors; G_ℓk·c_ka for a basis fit, G = diag(y_s)K). Then C = E[JᵀW_yJ] is
/// Σ_ab M_ab E[∇B_a∇B_bᵀ] with M = HᵀW_yH, a co-active pair's C_fg is the same form with the two functionals' rows, and the totals behind
/// the Poincaré check are Σ_ab M_ab·T_ab(v) over the term pairs' totals. The metric (decision 2) is a diagonal map L under a product
/// law, C = LᵀC_boxL: 1 for the box, 1/width for raw units, the law's standard deviation in box units for whitening.
/// </summary>
internal static class BassActiveSubspaces
{
    private const int Batch = 64;                  // a scalar fit's draws read from the gradient readout at once
    private const long BatchValues = 1L << 22;     // a functional fit's batch: about this many partials at once

    // What a readout shares: the inputs and the law, the continuous predictors, the metric's diagonal and the route.
    private sealed record Frame(BassInputs Inputs, ProductMeasure Measure, BassInputTransform Transform, int[] Variables, string[] Names,
                                string[] Excluded, double[] Widths, double[] Scale, ActiveSubspaceRoute Route, int Degree)
    {
        public int P => Variables.Length;

        public double[,] Map
        {
            get
            {
                var map = new double[P, P];
                for (int i = 0; i < P; i++) map[i, i] = Scale[i];
                return map;
            }
        }
    }

    // A draw: its terms and the coefficient of each term at each output location.
    private sealed record Surface(BassTerm[] Terms, double[,] H);

    // ── Active subspaces ──────────────────────────────────────────────────────────────────────────────────────────────

    public static ActiveSubspaceAnalysis Of(BassFit fit, BassData data, ActiveSubspaceSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        Fingerprint(fit, data);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        Frame frame = FrameOf(inputs, fit.Transform, fit.Config.Structure.Degree, spec.Measures, spec.Metric, spec.Route, spec.Points, spec.Dimension);
        bool functional = data.Response is FunctionalResponse;
        int m = data.Response.Locations;
        if (!functional && spec.OutputMetric is not null) throw new ArgumentException("An output metric applies to a functional fit's outputs; this fit's response is scalar.", nameof(spec));
        double[,] w = Metric(spec.OutputMetric, m);
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        (int[] selected, double[] weights) = Draws3.Select(spec.Draws, pooled.Length, fit.DrawWeights, nameof(spec));
        Func<int, Surface> surface = functional ? DirectFunctional(fit, data, inputs, pooled) : d => Scalar(fit, pooled[d]);
        var moments = BassMoments.On(inputs, frame.Degree, frame.Measure);

        var cBox = new double[selected.Length][,];
        var draws = new DrawSubspace[selected.Length];
        var forms = new (BassTerm[] Terms, double[,] M)[selected.Length];
        var surfaces = new Surface[selected.Length];
        double[][,]? error = null;
        for (int d = 0; d < selected.Length; d++)
        {
            surfaces[d] = surface(selected[d]);
            forms[d] = (surfaces[d].Terms, Quadratic(surfaces[d].H, w));
        }
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            for (int d = 0; d < selected.Length; d++)
            {
                cBox[d] = functional ? Symmetric(moments.GradientCross(forms[d].Terms, forms[d].Terms, forms[d].M))
                                     : moments.GradientOuterProduct(forms[d].Terms, Row(surfaces[d].H, 0));
                draws[d] = DrawSubspace.FromC(Congruence(cBox[d], frame.Scale), spec.Dimension);
            }
        }
        else
        {
            error = new double[selected.Length][,];
            double[,] root = Root(w);
            int rows = root.GetLength(0), p = frame.P;
            (PredictionCases cases, int n) = Points(fit.Transform, frame.Measure, spec.Points, spec.Seed);
            int batch = functional ? (int)Math.Max(1, BatchValues / ((long)n * m * p)) : Batch;
            double sqrtN = Math.Sqrt(n);
            var rj = new double[rows, p];
            for (int start = 0; start < selected.Length; start += batch)
            {
                int[] chunk = selected[start..Math.Min(selected.Length, start + batch)];
                BassGradient g = fit.Gradient(data, cases, chunk);
                for (int b = 0; b < chunk.Length; b++)
                {
                    int d = start + b;
                    // Each point's Jacobian in box units (∂f/∂x = width·∂f/∂z), mapped by W_y's root: its rows join the factor, in the
                    // metric, over √N, and its Gram matrix C's running sums.
                    var factor = new double[n * rows, p];
                    var sum = new double[p, p];
                    var square = new double[p, p];
                    for (int r = 0; r < n; r++)
                    {
                        for (int q = 0; q < rows; q++)
                            for (int i = 0; i < p; i++)
                            {
                                double[,,] partials = g.Partials[i];
                                double v = 0.0;
                                for (int l = 0; l < m; l++)
                                    if (root[q, l] != 0.0) v += root[q, l] * partials[b, r, l];
                                rj[q, i] = frame.Widths[i] * v;
                                factor[r * rows + q, i] = rj[q, i] * frame.Scale[i] / sqrtN;
                            }
                        for (int i = 0; i < p; i++)
                            for (int j = i; j < p; j++)
                            {
                                double v = 0.0;
                                for (int q = 0; q < rows; q++) v += rj[q, i] * rj[q, j];
                                sum[i, j] += v;
                                square[i, j] += v * v;
                            }
                    }
                    draws[d] = DrawSubspace.FromFactor(factor, spec.Dimension);
                    (cBox[d], double[,] se) = MeanAndError(sum, square, n, symmetric: true);
                    error[d] = Congruence(se, frame.Scale);
                }
            }
        }
        var posterior = SubspacePosterior.Of(draws, weights);
        double[,]? jacobi = !functional && LegendreOnBoxes(fit.Terms, pooled, frame, selected) ? Jacobi(fit, data, frame, selected) : null;
        PoincareCheck poincare = Poincare(frame, moments, forms, cBox, error, jacobi);
        return new ActiveSubspaceAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map,
                                          frame.Measure, frame.Route == ActiveSubspaceRoute.Sampled ? spec.Points : 0, posterior, error, poincare,
                                          Array.AsReadOnly(selected));
    }

    public static ActiveSubspaceAnalysis Of(BassBasisFit fit, ActiveSubspaceSpec spec)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(spec);
        (BassInputs inputs, BassFit first, int degree) = BasisInputs(fit);
        Frame frame = FrameOf(inputs, first.Transform, degree, spec.Measures, spec.Metric, spec.Route, spec.Points, spec.Dimension);
        if (frame.Route == ActiveSubspaceRoute.Sampled)
            throw new ArgumentException("A basis fit's C runs on the exact route: its sampled route waits for slice 4's joint samples.", nameof(spec));
        double[,] w = Metric(spec.OutputMetric, fit.Decomposition.Locations);
        (int[] selected, Func<int, Surface> surface) = BasisSurfaces(fit, spec.Draws, nameof(spec));
        var moments = BassMoments.On(inputs, degree, frame.Measure);
        var cBox = new double[selected.Length][,];
        var draws = new DrawSubspace[selected.Length];
        var forms = new (BassTerm[] Terms, double[,] M)[selected.Length];
        for (int d = 0; d < selected.Length; d++)
        {
            Surface s = surface(selected[d]);
            forms[d] = (s.Terms, Quadratic(s.H, w));
            cBox[d] = Symmetric(moments.GradientCross(s.Terms, s.Terms, forms[d].M));
            draws[d] = DrawSubspace.FromC(Congruence(cBox[d], frame.Scale), spec.Dimension);
        }
        double[] weights = [.. selected.Select(_ => 1.0)];
        var posterior = SubspacePosterior.Of(draws, weights);
        PoincareCheck poincare = Poincare(frame, moments, forms, cBox, null, null);
        return new ActiveSubspaceAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map,
                                          frame.Measure, 0, posterior, null, poincare, Array.AsReadOnly(selected));
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
        Frame frame = FrameOf(inputs, fit.Transform, fit.Config.Structure.Degree, spec.Measures, spec.Metric, spec.Route, spec.Points, 1);
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
            var moments = BassMoments.On(inputs, frame.Degree, frame.Measure);
            Func<int, Surface> sf = Surfaces(fit, data, inputs, pooledF), sg = Surfaces(other, otherData, otherInputs, pooledG);
            for (int s = 0; s < pairs.Length; s++)
            {
                Surface f = sf(selF[s]), g = sg(selG[s]);
                pairs[s] = Pair(moments, frame, f.Terms, Read(a, f.H), g.Terms, Read(b, g.H));
            }
        }
        else
        {
            (PredictionCases cases, int n) = Points(fit.Transform, frame.Measure, spec.Points, spec.Seed);
            errors = new CoActiveError[pairs.Length];
            long m = Math.Max(data.Response.Locations, otherData.Response.Locations);
            int batch = (int)Math.Clamp(BatchValues / (n * m * frame.P), 1, Batch);
            for (int start = 0; start < pairs.Length; start += batch)
            {
                int end = Math.Min(pairs.Length, start + batch);
                BassGradient gf = fit.Gradient(data, cases, selF[start..end]), gg = other.Gradient(otherData, cases, selG[start..end]);
                for (int s = 0; s < end - start; s++)
                    (pairs[start + s], errors[start + s]) = SampledPair(frame, n, Directional(gf, s, a), Directional(gg, s, b));
            }
        }
        return new CoActiveAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map, frame.Measure,
                                    frame.Route == ActiveSubspaceRoute.Sampled ? spec.Points : 0, CoActivePosterior.Of(pairs, weights),
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
        Frame frame = FrameOf(inputs, fit.Transform, fit.Config.Structure.Degree, spec.Measures, spec.Metric, spec.Route, spec.Points, 1);
        BassDraw[] pooled = [.. fit.Draws.SelectMany(x => x)];
        (int[] selected, double[] weights) = Draws3.Select(spec.Draws, pooled.Length, fit.DrawWeights, nameof(spec));
        var pairs = new CoActivePair[selected.Length];
        CoActiveError[]? errors = null;
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            var moments = BassMoments.On(inputs, frame.Degree, frame.Measure);
            Func<int, Surface> surface = DirectFunctional(fit, data, inputs, pooled);
            for (int s = 0; s < pairs.Length; s++)
            {
                Surface sf = surface(selected[s]);
                pairs[s] = Pair(moments, frame, sf.Terms, Contract(a, sf.H), sf.Terms, Contract(c, sf.H));
            }
        }
        else
        {
            (PredictionCases cases, int n) = Points(fit.Transform, frame.Measure, spec.Points, spec.Seed);
            errors = new CoActiveError[pairs.Length];
            int batch = (int)Math.Max(1, BatchValues / ((long)n * m * frame.P));
            for (int start = 0; start < pairs.Length; start += batch)
            {
                int end = Math.Min(pairs.Length, start + batch);
                BassGradient g = fit.Gradient(data, cases, selected[start..end]);
                for (int s = 0; s < end - start; s++) (pairs[start + s], errors[start + s]) = SampledPair(frame, n, Directional(g, s, a), Directional(g, s, c));
            }
        }
        return new CoActiveAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map, frame.Measure,
                                    frame.Route == ActiveSubspaceRoute.Sampled ? spec.Points : 0, CoActivePosterior.Of(pairs, weights),
                                    errors, Array.AsReadOnly(selected), null);
    }

    // Two outputs of one basis fit, over its joint draws, on the exact route.
    public static CoActiveAnalysis CoActive(BassBasisFit fit, CoActiveSpec spec)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(spec);
        if (spec.OtherDraws is not null) throw new ArgumentException("Two outputs of one fit pair each draw with itself; OtherDraws is for two fits.", nameof(spec));
        int m = fit.Decomposition.Locations;
        double[] a = Functional(spec.Output ?? throw new ArgumentException("A basis fit's co-active pair reads two outputs: Output is required.", nameof(spec)), m, "Output");
        double[] c = Functional(spec.OtherOutput ?? throw new ArgumentException("A basis fit's co-active pair reads two outputs: OtherOutput is required.", nameof(spec)), m, "OtherOutput");
        (BassInputs inputs, BassFit first, int degree) = BasisInputs(fit);
        Frame frame = FrameOf(inputs, first.Transform, degree, spec.Measures, spec.Metric, spec.Route, spec.Points, 1);
        if (frame.Route == ActiveSubspaceRoute.Sampled)
            throw new ArgumentException("A basis fit's co-active pairs run on the exact route: its sampled route waits for slice 4's joint samples.", nameof(spec));
        (int[] selected, Func<int, Surface> surface) = BasisSurfaces(fit, spec.Draws, nameof(spec));
        var moments = BassMoments.On(inputs, degree, frame.Measure);
        var pairs = new CoActivePair[selected.Length];
        for (int s = 0; s < pairs.Length; s++)
        {
            Surface sf = surface(selected[s]);
            pairs[s] = Pair(moments, frame, sf.Terms, Contract(a, sf.H), sf.Terms, Contract(c, sf.H));
        }
        return new CoActiveAnalysis(Array.AsReadOnly(frame.Names), Array.AsReadOnly(frame.Excluded), frame.Route, spec.Metric, frame.Map, frame.Measure,
                                    0, CoActivePosterior.Of(pairs, [.. selected.Select(_ => 1.0)]), null, Array.AsReadOnly(selected), null);
    }

    // C_f, C_g and C_fg on the exact route, in the metric.
    private static CoActivePair Pair(BassMoments moments, Frame frame, BassTerm[] termsF, double[] u, BassTerm[] termsG, double[] v)
        => new(Congruence(moments.GradientOuterProduct(termsF, u), frame.Scale), Congruence(moments.GradientOuterProduct(termsG, v), frame.Scale),
               Congruence(moments.GradientCross(termsF, u, termsG, v), frame.Scale));

    // C_f, C_g and C_fg from the two gradients at N shared points (box units after the width), in the metric, with each entry's
    // standard error.
    private static (CoActivePair Pair, CoActiveError Error) SampledPair(Frame frame, int n, Func<int, int, double> gradF, Func<int, int, double> gradG)
    {
        int p = frame.P;
        double[,] sf = new double[p, p], qf = new double[p, p], sg = new double[p, p], qg = new double[p, p], sfg = new double[p, p], qfg = new double[p, p];
        var f = new double[p];
        var g = new double[p];
        for (int r = 0; r < n; r++)
        {
            for (int i = 0; i < p; i++)
            {
                f[i] = frame.Widths[i] * gradF(r, i) * frame.Scale[i];
                g[i] = frame.Widths[i] * gradG(r, i) * frame.Scale[i];
            }
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    double fg = f[i] * g[j];
                    sfg[i, j] += fg;
                    qfg[i, j] += fg * fg;
                    if (j < i) continue;
                    double ff = f[i] * f[j], gg = g[i] * g[j];
                    sf[i, j] += ff;
                    qf[i, j] += ff * ff;
                    sg[i, j] += gg;
                    qg[i, j] += gg * gg;
                }
        }
        (double[,] cf, double[,] ef) = MeanAndError(sf, qf, n, symmetric: true);
        (double[,] cg, double[,] eg) = MeanAndError(sg, qg, n, symmetric: true);
        (double[,] cfg, double[,] efg) = MeanAndError(sfg, qfg, n, symmetric: false);
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
                                 InputMetric metric, ActiveSubspaceRoute route, int points, int dimension)
    {
        if (points < 2) throw new ArgumentOutOfRangeException(nameof(points), "The sampled route needs at least two points.");
        var measure = ProductMeasure.Resolve(transform, inputs, measures);
        int p = inputs.Cases.ContinuousCount;
        if (p == 0) throw new ArgumentException("C needs a continuous predictor: a categorical one has no derivative.", nameof(inputs));
        if (dimension < 1 || dimension > p) throw new ArgumentOutOfRangeException(nameof(dimension), $"The dimension k is from 1 to the {p} continuous predictors.");
        int[] variables = [.. Enumerable.Range(0, p).Select(measure.OfContinuous)];
        string[] names = [.. variables.Select(v => measure.Names[v])];
        string[] excluded = [.. Enumerable.Range(0, inputs.Cases.CategoricalCount).Select(c => measure.Names[measure.OfCategorical(c)])];
        double[] widths = [.. variables.Select(v => { var s = (RangeScale)transform.Predictors[v]; return s.Max - s.Min; })];
        double[] scale = metric switch
        {
            InputMetric.Box => [.. variables.Select(_ => 1.0)],
            InputMetric.Raw => [.. widths.Select(x => 1.0 / x)],
            InputMetric.Whitened => [.. variables.Select(v => StandardDeviation(measure.Continuous(v)!))],
            _ => throw new ArgumentOutOfRangeException(nameof(metric), $"No metric {metric}."),
        };
        return new Frame(inputs, measure, transform, variables, names, excluded, widths, scale,
                         route == ActiveSubspaceRoute.Auto ? ActiveSubspaceRoute.Exact : route, degree);
    }

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

    // R with W_y = RᵀR, from W_y's eigendecomposition; refused where W_y has a negative eigenvalue beyond rounding.
    private static double[,] Root(double[,] w)
    {
        int m = w.GetLength(0);
        EigenResult e = DenseEigen.DecomposeSymmetric(w);
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

    // LᵀCL for a diagonal L.
    private static double[,] Congruence(double[,] c, double[] scale)
    {
        int p = scale.Length;
        var z = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) z[i, j] = scale[i] * c[i, j] * scale[j];
        return z;
    }

    // The mean over N points of a matrix from its running sum and sum of squares, and each entry's standard error; a symmetric one's
    // upper triangle mirrored.
    private static (double[,] Mean, double[,] Error) MeanAndError(double[,] sum, double[,] square, int n, bool symmetric)
    {
        int p = sum.GetLength(0);
        var mean = new double[p, p];
        var error = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = symmetric ? i : 0; j < p; j++)
            {
                double m = sum[i, j] / n, variance = Math.Max(0.0, (square[i, j] - n * m * m) / (n - 1));
                mean[i, j] = m;
                error[i, j] = Math.Sqrt(variance / n);
                if (!symmetric) continue;
                mean[j, i] = m;
                error[j, i] = error[i, j];
            }
        return (mean, error);
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
    private static (PredictionCases Cases, int N) Points(BassInputTransform transform, ProductMeasure measure, int n, int seed)
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
        return (new PredictionCases(columns), n);
    }

    // ── Decision 6 ────────────────────────────────────────────────────────────────────────────────────────────────────

    // Poincaré per draw and continuous input with a uniform law, from totals Σ_ab M_ab·T_ab(v) over the term pairs' totals, the moment
    // route's on either route; and SG-4(b) where the Jacobi measure is given.
    private static PoincareCheck Poincare(Frame frame, BassMoments moments, (BassTerm[] Terms, double[,] M)[] forms, double[][,] cBox,
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
            bool uniform = frame.Measure.Declared[v] is UniformMeasure;
            ContinuousMeasure m = frame.Measure.Continuous(v)!;
            double constant = uniform ? (m.Upper - m.Lower) * (m.Upper - m.Lower) / (Math.PI * Math.PI) : double.NaN;
            for (int d = 0; d < draws; d++)
            {
                double t = total[d, i];
                bound[d, i] = constant * cBox[d][i, i];
                ratio[d, i] = uniform ? (bound[d, i] > 0.0 ? t / bound[d, i] : (t > 0.0 ? double.PositiveInfinity : 0.0)) : double.NaN;
                if (uniform)
                {
                    // Rounding on the exact route; three standard errors of the bound on the sampled route (the error is in the metric).
                    double slack = error is null ? 1e-9 * Math.Max(Math.Abs(t), Math.Abs(bound[d, i]))
                                                 : 3.0 * constant * error[d][i, i] / (frame.Scale[i] * frame.Scale[i]);
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
        if (frame.Measure.Count != frame.Variables.Length) return false;
        foreach (int v in frame.Variables)
        {
            ContinuousMeasure m = frame.Measure.Continuous(v)!;
            if (frame.Measure.Declared[v] is not UniformMeasure || m.Lower != 0.0 || m.Upper != 1.0) return false;   // the box, as FromCoefficients reads it
        }
        return true;
    }
}
