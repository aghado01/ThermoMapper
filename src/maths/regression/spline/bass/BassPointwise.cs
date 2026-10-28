using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's pointwise log density (OQ17; D04; stage 4's S9): log p(y_unit | draw) for every retained draw and every unit
/// of the training observations, a view over the draws that reduces them to the shared <see cref="PointwiseLogDensity"/>
/// the <see cref="PredictiveChecks"/> read. The draws come per chain in draw order, so PSIS-LOO reads each unit's
/// relative efficiency from their chains.
/// <list type="bullet">
/// <item><b>Units</b> partition the N observations, case by case and, within a case, location by location. By default
/// (<see cref="DefaultUnits"/>) each observation of a scalar response is a unit, and each case of a direct functional
/// response, its curve of n_r values, the hub's unit for functional data.</item>
/// <item><b>Retained parameters</b> (<see cref="PointwiseContent.RetainedParameters"/>, D04's default): the Normal
/// density of the unit's observations at the draw's coefficients and σ², the sum of their log terms, since the noise is
/// independent across observations. The structure, τ, the coefficients and σ² are global, so every unit conditions on
/// them.</item>
/// <item><b>The collapsed predictive</b> (<see cref="PointwiseContent.CollapsedPredictive"/>): the draw's structure and
/// τ kept, and the coefficients and σ² integrated out given the observations outside the unit, ln p(y_U | y_rest,
/// structure, τ), a multivariate Student t. It is the g-prior block's <see cref="GPriorConditional.LeaveOutLogDensity"/>
/// on the structure's explicit design [1, Z_c], under the fit's intercept convention and σ² prior, with the prior (its
/// centre, the Gram and τ) held at the full data's. It costs one explicit N × (M + 1) design and one factor per distinct
/// structure, which the draws at that structure share, then N solves per draw.</item>
/// </list>
/// λ enters neither content: it is the count prior's rate. Thread-safe: it reads the model and the table and writes
/// only its own buffers.
/// </summary>
public static class BassPointwise
{
    /// <summary>The default units: each observation of a scalar response, each case (its curve) of a direct functional
    /// one.</summary>
    public static UnitPartition DefaultUnits(BassModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        BassResponse response = model.Data.Response;
        return response is FunctionalResponse
            ? UnitPartition.ConsecutiveBlocks(response.Cases, response.Locations)
            : UnitPartition.EachObservation(response.Observations);
    }

    /// <summary>
    /// log p(y_unit | draw) for the retained draws <paramref name="draws"/> of <paramref name="model"/>, per chain in draw
    /// order, whose structures index <paramref name="terms"/>; the rows pool the chains in chain order.
    /// </summary>
    /// <param name="units">A partition of the N observations; null takes <see cref="DefaultUnits"/>.</param>
    /// <param name="content">The draw content; the retained parameters by default (D04).</param>
    public static PointwiseLogDensity LogDensity(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws,
                                                 UnitPartition? units = null,
                                                 PointwiseContent content = PointwiseContent.RetainedParameters)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(terms);
        ArgumentNullException.ThrowIfNull(draws);
        if (content is not (PointwiseContent.RetainedParameters or PointwiseContent.CollapsedPredictive))
            throw new ArgumentOutOfRangeException(nameof(content));
        UnitPartition partition = units ?? DefaultUnits(model);
        if (partition.Observations != model.Data.Response.Observations)
            throw new ArgumentException("The partition must cover the model's N observations, case by case and location by location.", nameof(units));
        bool collapsed = content == PointwiseContent.CollapsedPredictive;
        if (collapsed)
        {
            // The predictive conditions on the observations outside the unit, beyond the one a flat intercept's integral
            // spends; under Jeffreys that is σ²'s shape given the rest, (N − spent − m)/2, staying positive.
            int spent = model.Target.Intercept == GPriorCentre.Flat ? 1 : 0;
            for (int u = 0; u < partition.Units; u++)
                if (partition.Observations - spent - partition.Members[u].Length < 1)
                    throw new ArgumentException($"Unit {u} leaves no observation outside it to condition on"
                                                + (spent == 1 ? ", beyond the one the flat intercept spends." : "."), nameof(units));
        }
        (List<BassDraw> pooled, int[] chains) = BassDrawFit.Pool(draws, terms, collapsed, nameof(draws));
        var fit = new BassDrawFit(model, terms);
        double[,] values = collapsed ? Collapsed(model, fit, pooled, partition, nameof(draws)) : Retained(model, fit, pooled, partition);

        LatentTreatment treatment = collapsed ? LatentTreatment.IntegratedGivenOthers : LatentTreatment.Conditioned;
        var latents = new List<LatentReport> { new("structure", LatentScope.Global, LatentTreatment.Conditioned) };
        if (model.Target.SamplesTau) latents.Add(new LatentReport("prior precision τ", LatentScope.Global, LatentTreatment.Conditioned));
        latents.Add(new LatentReport("coefficients a", LatentScope.Global, treatment));
        latents.Add(new LatentReport("noise variance σ²", LatentScope.Global, treatment));
        return new PointwiseLogDensity(values, partition, latents, chains);
    }

    // Each draw's Normal log terms at its fit and σ², summed over each unit's observations.
    private static double[,] Retained(BassModel model, BassDrawFit fit, List<BassDraw> pooled, UnitPartition partition)
    {
        ReadOnlySpan<double> y = model.Data.Response.Values;
        int n = y.Length;
        var values = new double[pooled.Count, partition.Units];
        var f = new double[n];
        var log = new double[n];
        for (int s = 0; s < pooled.Count; s++)
        {
            BassDraw draw = pooled[s];
            fit.Evaluate(draw, f);
            double sigma2 = draw.Sigma2, normalizer = -0.5 * Math.Log(2.0 * Math.PI * sigma2);
            for (int i = 0; i < n; i++)
            {
                double e = y[i] - f[i];
                log[i] = normalizer - 0.5 * e * e / sigma2;
            }
            for (int u = 0; u < partition.Units; u++)
            {
                double sum = 0.0;
                foreach (int i in partition.Members[u]) sum += log[i];
                values[s, u] = sum;
            }
        }
        return values;
    }

    // The draws grouped by structure: one design, block and explicit design per structure, then each draw's leave-unit-out
    // predictive at its τ.
    private static double[,] Collapsed(BassModel model, BassDrawFit fit, List<BassDraw> pooled, UnitPartition partition, string parameter)
    {
        BassTarget target = model.Target;
        ReadOnlySpan<double> y = model.Data.Response.Values;
        var byStructure = new Dictionary<int[], List<int>>(StructureComparer.Instance);
        for (int s = 0; s < pooled.Count; s++)
        {
            int[] ids = pooled[s].Structure.TermIds;
            if (!byStructure.TryGetValue(ids, out List<int>? rows)) byStructure.Add(ids, rows = []);
            rows.Add(s);
        }

        var values = new double[pooled.Count, partition.Units];
        foreach ((int[] ids, List<int> rows) in byStructure)
        {
            var structure = new BassTerm[ids.Length];
            for (int j = 0; j < ids.Length; j++) structure[j] = fit.Term(ids[j]);
            BassDesign design = model.NewDesign();
            BassRefusal? refusal = design.Load(structure) ?? target.Refusal(BassStructureFit.Of(design, 0.0));
            if (refusal is { } reason) throw new ArgumentException($"A draw's structure lies outside F ({reason}).", parameter);
            CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                             design.FittedSquares, design.Observations, design.ResponseMean, target.Intercept);
            double[,] x = fit.CentredDesign(ids);
            foreach (int s in rows)
            {
                GPriorConditional conditional = block.At(pooled[s].Tau);
                for (int u = 0; u < partition.Units; u++)
                    values[s, u] = conditional.LeaveOutLogDensity(target.Sigma2, x, y, null, partition.Members[u]);
            }
        }
        return values;
    }

    private sealed class StructureComparer : IEqualityComparer<int[]>
    {
        public static readonly StructureComparer Instance = new();

        public bool Equals(int[]? a, int[]? b) => a is not null && b is not null && a.AsSpan().SequenceEqual(b);

        public int GetHashCode(int[] ids)
        {
            var hash = new HashCode();
            foreach (int id in ids) hash.Add(id);
            return hash.ToHashCode();
        }
    }
}

/// <summary>
/// Retained draws' training fits, f = a₀ + Σⱼ aⱼBⱼ at every observation, case by case and location by location, with
/// each term's case and output factors evaluated once per id and shared by the draws. Not thread-safe: one per call.
/// </summary>
internal sealed class BassDrawFit
{
    private readonly BassModel _model;
    private readonly BassTermTable _terms;
    private readonly Dictionary<int, (double[] Case, double[]? Output)> _columns = [];
    private readonly int _cases, _locations;

    public BassDrawFit(BassModel model, BassTermTable terms)
    {
        _model = model;
        _terms = terms;
        _cases = model.Data.Response.Cases;
        _locations = model.Data.Response.Locations;
    }

    /// <summary>
    /// The draws of every chain, pooled in chain order, with each chain's length; refused: no chain, an empty chain, a
    /// null draw, ids outside the table or not ascending and distinct (an unordered collection, D01), other than one
    /// finite coefficient per term plus the intercept's, a σ² that is not positive and finite and, where τ is read, a τ
    /// that is not.
    /// </summary>
    public static (List<BassDraw> Pooled, int[] Chains) Pool(IReadOnlyList<IReadOnlyList<BassDraw>> draws, BassTermTable terms,
                                                            bool readsTau, string parameter)
    {
        if (draws is null) throw new ArgumentNullException(parameter);
        if (draws.Count == 0) throw new ArgumentException("At least one chain of draws is required.", parameter);
        var pooled = new List<BassDraw>();
        var chains = new int[draws.Count];
        for (int c = 0; c < draws.Count; c++)
        {
            IReadOnlyList<BassDraw> chain = draws[c] ?? throw new ArgumentException($"Chain {c} is null.", parameter);
            if (chain.Count == 0) throw new ArgumentException($"Chain {c} holds no draw; every chain must hold at least one.", parameter);
            chains[c] = chain.Count;
            foreach (BassDraw draw in chain)
            {
                Check(draw, terms, readsTau, parameter);
                pooled.Add(draw);
            }
        }
        return (pooled, chains);
    }

    private static void Check(BassDraw? draw, BassTermTable terms, bool readsTau, string parameter)
    {
        if (draw is null) throw new ArgumentException("A draw is null.", parameter);
        int[] ids = draw.Structure.TermIds ?? throw new ArgumentException("A draw's structure holds no ids.", parameter);
        for (int j = 0; j < ids.Length; j++)
        {
            if ((uint)ids[j] >= (uint)terms.Count) throw new ArgumentException($"Term id {ids[j]} is not in the table.", parameter);
            if (j > 0 && ids[j] <= ids[j - 1]) throw new ArgumentException("A draw's term ids must be ascending and distinct.", parameter);
        }
        if (draw.Coefficients is null || draw.Coefficients.Length != ids.Length + 1)
            throw new ArgumentException("A draw needs one coefficient per term and the intercept's.", parameter);
        foreach (double a in draw.Coefficients)
            if (!double.IsFinite(a)) throw new ArgumentException("A draw's coefficients must be finite.", parameter);
        if (!(draw.Sigma2 > 0.0) || double.IsPositiveInfinity(draw.Sigma2))
            throw new ArgumentException("A draw's σ² must be positive and finite.", parameter);
        if (readsTau && (!(draw.Tau > 0.0) || double.IsPositiveInfinity(draw.Tau)))
            throw new ArgumentException("A draw's τ must be positive and finite.", parameter);
    }

    /// <summary>The term with id <paramref name="id"/>, checked against the model's inputs.</summary>
    public BassTerm Term(int id)
    {
        Columns(id);
        return _terms[id];
    }

    /// <summary>Term <paramref name="id"/>'s case factor and, for a direct functional response, its output factor.</summary>
    public (double[] Case, double[]? Output) Columns(int id)
    {
        if (_columns.TryGetValue(id, out var cached)) return cached;
        BassTerm term = _terms[id];
        BassBasis basis = _model.Basis;
        try
        {
            basis.Validate(term);
        }
        catch (ArgumentException e)
        {
            throw new ArgumentException($"Term {id} ({term}) does not fit the model's inputs: {e.Message}", "terms", e);
        }
        var a = new double[_cases];
        basis.CaseColumn(term, basis.Inputs.Cases, a);
        double[]? h = null;
        if (basis.Inputs.Outputs is { } outputs)
        {
            h = new double[_locations];
            basis.OutputColumn(term, outputs, h);
        }
        _columns.Add(id, (a, h));
        return (a, h);
    }

    /// <summary>The draw's fit at every observation, into <paramref name="fit"/>.</summary>
    public void Evaluate(BassDraw draw, Span<double> fit)
    {
        int[] ids = draw.Structure.TermIds;
        fit.Fill(draw.Coefficients[0]);
        for (int q = 0; q < ids.Length; q++)
        {
            (double[] a, double[]? h) = Columns(ids[q]);
            double c = draw.Coefficients[q + 1];
            if (h is null)
            {
                for (int i = 0; i < _cases; i++) fit[i] += c * a[i];
                continue;
            }
            for (int i = 0; i < _cases; i++)
            {
                double ca = c * a[i];
                if (ca == 0.0) continue;
                Span<double> row = fit.Slice(i * _locations, _locations);
                for (int r = 0; r < _locations; r++) row[r] += ca * h[r];
            }
        }
    }

    /// <summary>
    /// The explicit design [1, Z_c] of the structure <paramref name="ids"/>, in their order: the intercept's column of
    /// ones, then each term's column centred over all N observations, a − ā for a scalar response and
    /// a(xᵢ)h(t_r) − āh̄ for a direct functional one, as the design layer centres them.
    /// </summary>
    public double[,] CentredDesign(int[] ids)
    {
        int n = _cases * _locations;
        var x = new double[n, ids.Length + 1];
        for (int i = 0; i < n; i++) x[i, 0] = 1.0;
        for (int j = 0; j < ids.Length; j++)
        {
            (double[] a, double[]? h) = Columns(ids[j]);
            double aBar = PairwiseSum.Of(a) / _cases;
            if (h is null)
            {
                for (int i = 0; i < _cases; i++) x[i, j + 1] = a[i] - aBar;
                continue;
            }
            double level = aBar * (PairwiseSum.Of(h) / _locations);
            for (int i = 0; i < _cases; i++)
                for (int r = 0; r < _locations; r++) x[i * _locations + r, j + 1] = a[i] * h[r] - level;
        }
        return x;
    }
}
