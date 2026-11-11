using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Rotate-and-refit's frame (E8, slice 7; 2307.11241 projects and refits in §4.1.3 and §4.3): the d leading eigenvectors W₁ of a
/// posterior's E[C], read in the analysis's metric, and the active variables they define,
/// u = W₁ᵀL⁻¹x_box = Aᵀx_box with A = L⁻ᵀW₁ (p × d), x_box each continuous input in its box's units and L the analysis's map. For the
/// box metric u = W₁ᵀx_box; for whitening, the paper's Σ^(−1/2)x up to a shift, which BASS's rescaling makes immaterial.
/// <list type="bullet">
/// <item><b>The frame</b> is E[C]'s, by the user's decision of 2026-11-10: for any frame, the posterior mean of each draw's energy along
/// the dropped directions is Σ_{j>d} w_jᵀE[C]w_j, which E[C]'s trailing eigenvectors minimize. BASS ignores an input's sign and the
/// inputs' order, so W₁ fixes the refit; no rotation within the span is taken. The eigengaps within the active block are reported, since
/// near-equal eigenvalues leave their eigenvectors arbitrary within their span.</item>
/// <item><b>The domains.</b> The input box [0, 1]^p maps onto a zonotope whose bounding box is
/// [Σᵢ min(0, aᵢⱼ), Σᵢ max(0, aᵢⱼ)] per active variable, each end padded by 10⁻¹⁴ of Σᵢ|aᵢⱼ| so that a case's rounding never leaves
/// it. Each active variable takes it as its declared domain (the user's decision 4), so every image of the input box lies in the
/// refit's box. In high p it can be far wider than the projected data's spread, which <see cref="Spreads"/> reports.</item>
/// <item><b>The law.</b> The analysis's law pushes forward to the active variables, in general not as a product, so the refit's
/// readouts take E8's sampled route on a pushforward joint sample (<see cref="Pushforward(int, int)"/>).</item>
/// </list>
/// The categorical predictors are carried unchanged; the refit's continuous inputs are the d active variables alone (the user's decision
/// 3; d = p rotates without dropping).
/// </summary>
public sealed class BassProjection
{
    private readonly BassInputTransform _transform;
    private readonly RangeScale[] _scales;
    private readonly double[,] _a;

    private BassProjection(BassInputTransform transform, RangeScale[] scales, string[] inputs, string[] categorical, string[] names, InputMetric metric,
                           double[,] frame, double[,] a, double[] eigenvalues, double[] gaps, InputDomain[] domains, ProductMeasure? measure)
    {
        _transform = transform;
        _scales = scales;
        _a = a;
        Inputs = Array.AsReadOnly(inputs);
        Categorical = Array.AsReadOnly(categorical);
        Names = Array.AsReadOnly(names);
        Metric = metric;
        Frame = frame;
        Eigenvalues = Array.AsReadOnly(eigenvalues);
        Gaps = Array.AsReadOnly(gaps);
        Domains = Array.AsReadOnly(domains);
        Measure = measure;
    }

    /// <summary>The continuous predictors projected, in the analysis's order.</summary>
    public IReadOnlyList<string> Inputs { get; }

    /// <summary>The categorical predictors, carried unchanged.</summary>
    public IReadOnlyList<string> Categorical { get; }

    /// <summary>The active variables' names, u1 to ud unless declared.</summary>
    public IReadOnlyList<string> Names { get; }

    /// <summary>The metric the frame was read in.</summary>
    public InputMetric Metric { get; }

    /// <summary>d.</summary>
    public int Dimension => Names.Count;

    /// <summary>W₁, p × d: E[C]'s leading eigenvectors as columns, in the metric.</summary>
    public double[,] Frame { get; }

    /// <summary>A = L⁻ᵀW₁, p × d: u = Aᵀx_box.</summary>
    public double[,] Map => (double[,])_a.Clone();

    /// <summary>E[C]'s eigenvalues, descending, all p.</summary>
    public IReadOnlyList<double> Eigenvalues { get; }

    /// <summary>The eigengaps (λ_j − λ_{j+1})/λ₁ for j = 1 to d: within the active block for j &lt; d, and at its edge for j = d
    /// (λ_{p+1} = 0). A small gap within the block leaves the two eigenvectors arbitrary within their span.</summary>
    public IReadOnlyList<double> Gaps { get; }

    /// <summary>Each active variable's declared domain, the zonotope's bounding box.</summary>
    public IReadOnlyList<InputDomain> Domains { get; }

    /// <summary>The analysis's product law, in the native fit's transformed units; null when it read a joint sample.</summary>
    public ProductMeasure? Measure { get; }

    /// <summary>
    /// The frame of <paramref name="analysis"/>'s E[C] at dimension <paramref name="dimension"/>, for the fit whose transform is
    /// <paramref name="transform"/> (a fit's <see cref="BassFit.Transform"/>, or a basis fit's component's). Refused: an analysis
    /// whose continuous and categorical predictors are not the transform's, in its order, d outside 1 to p, and names that are not
    /// d distinct new ones.
    /// </summary>
    public static BassProjection Of(BassInputTransform transform, ActiveSubspaceAnalysis analysis, int dimension, IReadOnlyList<string>? names = null)
    {
        ArgumentNullException.ThrowIfNull(transform);
        ArgumentNullException.ThrowIfNull(analysis);
        string[] continuous = [.. transform.Predictors.OfType<RangeScale>().Select(s => s.Name)];
        string[] categorical = [.. transform.Predictors.OfType<TrainingLevels>().Select(s => s.Name)];
        if (!continuous.SequenceEqual(analysis.Inputs) || !categorical.SequenceEqual(analysis.Excluded))
            throw new ArgumentException("The analysis was read on other inputs: its continuous and categorical predictors must be the fit's, in its order.", nameof(analysis));
        int p = continuous.Length;
        if (dimension < 1 || dimension > p) throw new ArgumentOutOfRangeException(nameof(dimension), $"d is from 1 to the {p} continuous predictors.");
        string[] active = names is null ? [.. Enumerable.Range(1, dimension).Select(j => $"u{j}")] : [.. names];
        if (active.Length != dimension || active.Distinct(StringComparer.Ordinal).Count() != dimension || active.Any(string.IsNullOrWhiteSpace))
            throw new ArgumentException($"The active variables need {dimension} distinct names.", nameof(names));
        if (active.Intersect(categorical, StringComparer.Ordinal).Any())
            throw new ArgumentException("An active variable's name is a categorical predictor's.", nameof(names));

        SubspacePosterior posterior = analysis.Posterior;
        double[] eigenvalues = [.. posterior.MeanCEigenvalues];
        var frame = new double[p, dimension];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < dimension; j++) frame[i, j] = posterior.MeanCVectors[i, j];
        double top = eigenvalues[0];
        double[] gaps = [.. Enumerable.Range(0, dimension).Select(j => top > 0.0 ? (eigenvalues[j] - (j + 1 < p ? eigenvalues[j + 1] : 0.0)) / top : 0.0)];
        double[,] a = Multiply(InverseTranspose(analysis.Map), frame);

        var domains = new InputDomain[dimension];
        for (int j = 0; j < dimension; j++)
        {
            double lo = 0.0, hi = 0.0, size = 0.0;
            for (int i = 0; i < p; i++)
            {
                lo += Math.Min(0.0, a[i, j]);
                hi += Math.Max(0.0, a[i, j]);
                size += Math.Abs(a[i, j]);
            }
            double pad = 1e-14 * size;
            domains[j] = new InputDomain(lo - pad, hi + pad);
        }
        RangeScale[] scales = [.. transform.Predictors.OfType<RangeScale>()];
        return new BassProjection(transform, scales, continuous, categorical, active, analysis.Metric, frame, a, eigenvalues, gaps, domains, analysis.Measure);
    }

    /// <summary>The active variables' values at <paramref name="cases"/>: one row per case, one column per active variable.</summary>
    public double[,] Project(PredictionCases cases)
    {
        ArgumentNullException.ThrowIfNull(cases);
        double[][] columns = ContinuousColumns(cases);
        return Apply(columns, cases.Cases);
    }

    /// <summary>
    /// <paramref name="cases"/> in the refit's inputs: the active variables, then the categorical columns unchanged; every other
    /// column must be a projected continuous predictor's.
    /// </summary>
    public PredictionCases ProjectCases(PredictionCases cases)
    {
        double[,] u = Project(cases);
        var columns = new List<PredictionColumn>();
        for (int j = 0; j < Dimension; j++) columns.Add(new ContinuousColumn(Names[j], Column(u, j)));
        foreach (PredictionColumn c in cases.Columns)
            if (c is CategoricalColumn) columns.Add(c);
        return new PredictionCases(columns);
    }

    /// <summary>
    /// <paramref name="data"/> in the refit's inputs (its training data, or the same cases): the active variables, each with its
    /// bounding box as its declared domain, then the categorical predictors, the response, the groupings and the case weights unchanged.
    /// </summary>
    public BassData Project(BassData data)
    {
        ArgumentNullException.ThrowIfNull(data);
        string[] continuous = [.. data.Predictors.OfType<ContinuousPredictor>().Select(c => c.Name)];
        string[] categorical = [.. data.Predictors.OfType<CategoricalPredictor>().Select(c => c.Name)];
        if (!continuous.SequenceEqual(Inputs) || !categorical.SequenceEqual(Categorical))
            throw new ArgumentException("The data's predictors are not the projection's: its continuous and categorical predictors, in its order.", nameof(data));
        double[][] columns = [.. data.Predictors.OfType<ContinuousPredictor>().Select(c => c.Values.ToArray())];
        double[,] u = Apply(columns, columns[0].Length);
        var predictors = new List<BassPredictor>();
        for (int j = 0; j < Dimension; j++) predictors.Add(new ContinuousPredictor(Names[j], Column(u, j), Domains[j]));
        predictors.AddRange(data.Predictors.OfType<CategoricalPredictor>());
        return new BassData(predictors, data.Response, data.Groupings, data.Weights);
    }

    /// <summary>
    /// <paramref name="points"/> points of the analysis's product law, drawn once from <paramref name="seed"/> as E8's sampled route
    /// draws them, pushed forward to the refit's inputs as a joint sample of equal weights. Refused when the analysis read a joint
    /// sample: push that sample forward instead.
    /// </summary>
    public JointSample Pushforward(int points, int seed)
    {
        if (Measure is null) throw new InvalidOperationException("The analysis read a joint sample, which has no product law to draw from: push the sample forward.");
        if (points < 2) throw new ArgumentOutOfRangeException(nameof(points), "A joint sample needs at least two points.");
        return new JointSample(ProjectCases(BassActiveSubspaces.Points(_transform, Measure, points, seed)));
    }

    /// <summary><paramref name="sample"/>, a joint sample of the native inputs, pushed forward to the refit's inputs with its weights.</summary>
    public JointSample Pushforward(JointSample sample)
    {
        ArgumentNullException.ThrowIfNull(sample);
        return new JointSample(ProjectCases(sample.Points ?? throw new ArgumentException("A joint sample needs its points.", nameof(sample))), sample.Weights);
    }

    /// <summary>
    /// Each active variable's projected training range against its domain: the domain's width over the range's, at least 1. In high p
    /// a dense direction's zonotope reaches far past the cases, so the ratio grows (decision 4's gate reads it).
    /// </summary>
    public double[] Spreads(BassData data)
    {
        BassData projected = Project(data);
        var spreads = new double[Dimension];
        for (int j = 0; j < Dimension; j++)
        {
            double[] u = ((ContinuousPredictor)projected.Predictors[j]).Values.ToArray();
            spreads[j] = (Domains[j].Upper - Domains[j].Lower) / (u.Max() - u.Min());
        }
        return spreads;
    }

    // Each projected continuous predictor's column of the cases, in raw units; refused when one is missing.
    private double[][] ContinuousColumns(PredictionCases cases)
    {
        var byName = new Dictionary<string, PredictionColumn>(StringComparer.Ordinal);
        foreach (PredictionColumn c in cases.Columns) byName[c.Name] = c;
        var columns = new double[Inputs.Count][];
        for (int i = 0; i < Inputs.Count; i++)
        {
            if (!byName.TryGetValue(Inputs[i], out PredictionColumn? c) || c is not ContinuousColumn continuous)
                throw new ArgumentException($"The cases have no continuous column '{Inputs[i]}'.", nameof(cases));
            columns[i] = continuous.Values.ToArray();
        }
        foreach (PredictionColumn c in cases.Columns)
            if (c is ContinuousColumn && !Inputs.Contains(c.Name))
                throw new ArgumentException($"The column '{c.Name}' is no projected predictor.", nameof(cases));
        return columns;
    }

    // u = Aᵀx_box per case, x_box each raw value through its predictor's box.
    private double[,] Apply(double[][] columns, int n)
    {
        int p = Inputs.Count, d = Dimension;
        var u = new double[n, d];
        var x = new double[p];
        for (int r = 0; r < n; r++)
        {
            for (int i = 0; i < p; i++) x[i] = _scales[i].Scale(columns[i][r]);
            for (int j = 0; j < d; j++)
            {
                double s = 0.0;
                for (int i = 0; i < p; i++) s += _a[i, j] * x[i];
                u[r, j] = s;
            }
        }
        return u;
    }

    private static double[] Column(double[,] m, int j) => [.. Enumerable.Range(0, m.GetLength(0)).Select(r => m[r, j])];

    // L⁻ᵀ for the analysis's map: diagonal, or a symmetric positive definite root (whitening under a joint sample), by its eigenvalues.
    private static double[,] InverseTranspose(double[,] l)
    {
        int p = l.GetLength(0);
        var inverse = new double[p, p];
        bool diagonal = true;
        for (int i = 0; i < p && diagonal; i++)
            for (int j = 0; j < p; j++)
                if (i != j && l[i, j] != 0.0) { diagonal = false; break; }
        if (diagonal)
        {
            for (int i = 0; i < p; i++) inverse[i, i] = 1.0 / l[i, i];
            return inverse;
        }
        EigenResult e = ScaledEigen.Of(l);
        for (int q = 0; q < p; q++)
        {
            double s = 1.0 / e.Eigenvalues[q];
            double[] v = e.Eigenvectors[q];
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) inverse[i, j] += s * v[i] * v[j];
        }
        return inverse;
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), inner = a.GetLength(1), q = b.GetLength(1);
        var c = new double[n, q];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < inner; k++)
                for (int j = 0; j < q; j++) c[i, j] += a[i, k] * b[k, j];
        return c;
    }
}
