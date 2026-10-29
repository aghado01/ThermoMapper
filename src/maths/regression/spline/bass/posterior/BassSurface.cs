using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A fit's draws evaluated at any cases and, for a direct functional fit, any output locations: f = a₀ + Σⱼ aⱼBⱼ, with
/// each term's case and output factors evaluated once per id through the training basis (its knots, levels and degree)
/// and shared by the draws. The accumulation is <see cref="BassDrawFit.Evaluate"/>'s, so at the training inputs the
/// values are the stage-4 fitted values bit for bit. Not thread-safe: one per call.
/// </summary>
internal sealed class BassSurface
{
    private readonly BassBasis _basis;
    private readonly BassTermTable _terms;
    private readonly CaseInputs _cases;
    private readonly OutputInputs? _outputs;
    private readonly Dictionary<int, (double[] Case, double[]? Output)> _columns = [];

    /// <summary>The surface of <paramref name="basis"/>'s fit, whose draws index <paramref name="terms"/>, at
    /// <paramref name="cases"/> and, for a direct functional fit, <paramref name="outputs"/>.</summary>
    public BassSurface(BassBasis basis, BassTermTable terms, CaseInputs cases, OutputInputs? outputs)
    {
        if ((basis.Inputs.Outputs is null) != (outputs is null))
            throw new ArgumentException("A direct functional fit is evaluated at output locations, a scalar one at none.", nameof(outputs));
        _basis = basis;
        _terms = terms;
        _cases = cases;
        _outputs = outputs;
    }

    /// <summary>The cases.</summary>
    public int Cases => _cases.Cases;

    /// <summary>The locations: 1 for a scalar fit.</summary>
    public int Locations => _outputs?.Locations ?? 1;

    /// <summary>Term <paramref name="id"/>'s case factor and, for a direct functional fit, its output factor.</summary>
    public (double[] Case, double[]? Output) Columns(int id)
    {
        if (_columns.TryGetValue(id, out var cached)) return cached;
        BassTerm term = _terms[id];
        try
        {
            _basis.Validate(term);
        }
        catch (ArgumentException e)
        {
            throw new ArgumentException($"Term {id} ({term}) does not fit the model's inputs: {e.Message}", "terms", e);
        }
        var a = new double[_cases.Cases];
        _basis.CaseColumn(term, _cases, a);
        double[]? h = null;
        if (_outputs is { } outputs)
        {
            h = new double[outputs.Locations];
            _basis.OutputColumn(term, outputs, h);
        }
        _columns.Add(id, (a, h));
        return (a, h);
    }

    /// <summary>The function of <paramref name="coefficients"/> (the intercept, then one per id of
    /// <paramref name="ids"/>) at every case and location, case by case, into <paramref name="f"/>.</summary>
    public void Evaluate(ReadOnlySpan<int> ids, ReadOnlySpan<double> coefficients, Span<double> f)
    {
        int cases = Cases, locations = Locations;
        f.Fill(coefficients[0]);
        for (int q = 0; q < ids.Length; q++)
        {
            (double[] a, double[]? h) = Columns(ids[q]);
            double c = coefficients[q + 1];
            if (h is null)
            {
                for (int i = 0; i < cases; i++) f[i] += c * a[i];
                continue;
            }
            for (int i = 0; i < cases; i++)
            {
                double ca = c * a[i];
                if (ca == 0.0) continue;
                Span<double> row = f.Slice(i * locations, locations);
                for (int r = 0; r < locations; r++) row[r] += ca * h[r];
            }
        }
    }
}

/// <summary>
/// The Rao–Blackwellized coefficients of a draw (D07, *(2026-10-29)*): E[a | structure, τ, y], the g-prior's posterior
/// mean given the draw's structure and τ, in the uncentred basis the draws use. The design and block are built once per
/// structure from the training data, and the mean does not depend on σ². Not thread-safe: one per call.
/// </summary>
internal sealed class BassConditionalMeans
{
    private readonly BassModel _model;
    private readonly BassDrawFit _fit;
    private readonly Dictionary<int[], (CollapsedGPrior Block, double[] ColumnMeans)> _blocks = new(IdsComparer.Instance);

    public BassConditionalMeans(BassModel model, BassTermTable terms)
    {
        _model = model;
        _fit = new BassDrawFit(model, terms);
    }

    /// <summary>a₀ and one coefficient per id of <paramref name="draw"/>'s structure, at its τ.</summary>
    public double[] Coefficients(BassDraw draw)
    {
        int[] ids = draw.Structure.TermIds;
        if (!_blocks.TryGetValue(ids, out var entry))
        {
            var structure = new BassTerm[ids.Length];
            for (int j = 0; j < ids.Length; j++) structure[j] = _fit.Term(ids[j]);
            BassDesign design = _model.NewDesign();
            BassRefusal? refusal = design.Load(structure) ?? _model.Target.Refusal(BassStructureFit.Of(design, 0.0));
            if (refusal is { } reason) throw new ArgumentException($"A draw's structure lies outside F ({reason}).", "draws");
            CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                             design.FittedSquares, design.Observations, design.ResponseMean,
                                                                             _model.Target.Intercept);
            var means = new double[ids.Length];
            for (int j = 0; j < ids.Length; j++)
            {
                (double[] a, double[]? h) = _fit.Columns(ids[j]);
                double mean = PairwiseSum.Of(a) / a.Length;
                means[j] = h is null ? mean : mean * (PairwiseSum.Of(h) / h.Length);
            }
            entry = (block, means);
            _blocks.Add(ids, entry);
        }
        ReadOnlySpan<double> centred = entry.Block.At(draw.Tau).Mean;   // (α, β_c) in the ids' order
        var coefficients = new double[ids.Length + 1];
        double intercept = centred[0];
        for (int j = 0; j < ids.Length; j++)
        {
            coefficients[j + 1] = centred[j + 1];
            intercept -= centred[j + 1] * entry.ColumnMeans[j];
        }
        coefficients[0] = intercept;
        return coefficients;
    }

    private sealed class IdsComparer : IEqualityComparer<int[]>
    {
        public static readonly IdsComparer Instance = new();

        public bool Equals(int[]? a, int[]? b) => a is not null && b is not null && a.AsSpan().SequenceEqual(b);

        public int GetHashCode(int[] ids)
        {
            var hash = new HashCode();
            foreach (int id in ids) hash.Add(id);
            return hash.ToHashCode();
        }
    }
}
