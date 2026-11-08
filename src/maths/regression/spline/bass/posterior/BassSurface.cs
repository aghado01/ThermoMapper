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
    private readonly Dictionary<(int Id, int Variable), double[]> _derivatives = [];

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

    /// <summary>Term <paramref name="id"/>'s case factor's partial derivative in continuous predictor <paramref name="variable"/>
    /// (its position among the continuous predictors), in the transformed units (<see cref="BassBasis.CaseDerivative"/>).</summary>
    public double[] DerivativeColumn(int id, int variable)
    {
        if (_derivatives.TryGetValue((id, variable), out double[]? cached)) return cached;
        Columns(id);   // validates the term
        var a = new double[_cases.Cases];
        _basis.CaseDerivative(_terms[id], variable, _cases, a);
        _derivatives.Add((id, variable), a);
        return a;
    }

    /// <summary>The function of <paramref name="coefficients"/> (the intercept, then one per id of
    /// <paramref name="ids"/>) at every case and location, case by case, into <paramref name="f"/>.</summary>
    public void Evaluate(ReadOnlySpan<int> ids, ReadOnlySpan<double> coefficients, Span<double> f) => Accumulate(ids, coefficients, -1, f);

    /// <summary>The function's partial derivative in continuous predictor <paramref name="variable"/>, in the transformed units,
    /// at every case and location, into <paramref name="f"/>: the intercept drops, and each term contributes its derivative
    /// column times its output factor.</summary>
    public void EvaluateDerivative(ReadOnlySpan<int> ids, ReadOnlySpan<double> coefficients, int variable, Span<double> f)
    {
        if (variable < 0) throw new ArgumentOutOfRangeException(nameof(variable), "A continuous predictor's position is nonnegative.");
        Accumulate(ids, coefficients, variable, f);
    }

    // Σ coefficients × columns: the values' (variable −1, from the intercept) or one variable's derivatives (from 0).
    private void Accumulate(ReadOnlySpan<int> ids, ReadOnlySpan<double> coefficients, int variable, Span<double> f)
    {
        int cases = Cases, locations = Locations;
        f.Fill(variable < 0 ? coefficients[0] : 0.0);
        for (int q = 0; q < ids.Length; q++)
        {
            (double[] values, double[]? h) = Columns(ids[q]);
            double[] a = variable < 0 ? values : DerivativeColumn(ids[q], variable);
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
    private readonly Dictionary<(int[] Ids, WorkingResponse Working), (Func<double, double, GPriorConditional> Block, double[] ColumnMeans)> _blocks = new(BassDrawGroupComparer.Instance);

    public BassConditionalMeans(BassModel model, BassTermTable terms)
    {
        _model = model;
        _fit = new BassDrawFit(model, terms);
    }

    /// <summary>a₀ and one coefficient per id of <paramref name="draw"/>'s structure, at its τ.</summary>
    public double[] Coefficients(BassDraw draw)
    {
        int[] ids = draw.Structure.TermIds;
        WorkingResponse working = draw.WorkingOn(_model);
        if (!_blocks.TryGetValue((ids, working), out var entry))
        {
            var structure = new BassTerm[ids.Length];
            for (int j = 0; j < ids.Length; j++) structure[j] = _fit.Term(ids[j]);
            BassDesign design = _model.NewDesign(working);
            BassRefusal? refusal = design.Load(structure) ?? _model.Target.Refusal(BassStructureFit.Of(design, 0.0));
            if (refusal is { } reason) throw new ArgumentException($"A draw's structure lies outside F ({reason}).", "draws");
            Func<double, double, GPriorConditional> block = _model.Target.CoefficientBlock(design);
            var means = new double[ids.Length];
            for (int j = 0; j < ids.Length; j++)
            {
                (double[] a, double[]? h) = _fit.Columns(ids[j]);
                means[j] = working.ColumnMean(a, h);
            }
            entry = (block, means);
            _blocks.Add((ids, working), entry);
        }
        ReadOnlySpan<double> centred = entry.Block(draw.Tau, 1.0).Mean;   // (α, β_c) in the ids' order
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
}
