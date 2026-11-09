using System;
using System.Collections.Generic;
using Maths.Distributions;
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
/// structure, τ), a multivariate Student t. It is the coefficient block's <see cref="GPriorConditional.LeaveOutLogDensity"/>
/// on the structure's explicit design [1, Z_c], under the fit's coefficient form, intercept convention and σ² prior, with
/// the prior (its centre, the Gram or the ridge's scaling, and τ) held at the full data's. It costs one explicit
/// N × (M + 1) design and one block per distinct structure, which the draws at that structure share, then N solves per
/// draw. The g-form's precision is the Gram at every τ, so the units' quadratic forms are taken once per structure; the
/// ridge's, bG_c + τD, moves with τ, so they are taken once per distinct τ.</item>
/// </list>
/// λ enters neither content: it is the count prior's rate. A draw composed on a working response reads it, the moved ones of a chain
/// each their own (<see cref="BassDraw.Working"/>), and the offset and weights are then reported as latents conditioned on. A model with
/// a noise law (<see cref="IBassObservationModel.Law"/>) takes the law's density instead, each observation's latent integrated out
/// (OQ17's scope rule), at the draw's σ² and its latents' globals; its collapsed predictive is refused until it has a check of its own.
/// Thread-safe: it reads the model and the table and writes only its own buffers.
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
    /// order, whose structures index <paramref name="terms"/>; the rows pool the chains in chain order. Weighted draws, an
    /// SMC fit's replicate populations, carry <paramref name="weights"/> (<see cref="PointwiseLogDensity.Weights"/>) in
    /// place of chain lengths: their particles share ancestry, so the checks' Monte Carlo errors, which treat them as
    /// independent, read optimistic.
    /// </summary>
    /// <param name="units">A partition of the N observations; null takes <see cref="DefaultUnits"/>.</param>
    /// <param name="content">The draw content; the retained parameters by default (D04).</param>
    /// <param name="weights">One weight per draw, pooled in chain order; null for equal weights (MCMC).</param>
    /// <param name="quadrature">Takes a noise law whose density is an integral over each observation's latent
    /// (<see cref="IBassNoiseLaw.DensityByQuadrature"/>), N × draws integrals of hundreds of nodes each; refused without it.</param>
    public static PointwiseLogDensity LogDensity(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws,
                                                 UnitPartition? units = null,
                                                 PointwiseContent content = PointwiseContent.RetainedParameters,
                                                 IReadOnlyList<double>? weights = null, bool quadrature = false)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(terms);
        ArgumentNullException.ThrowIfNull(draws);
        if (content is not (PointwiseContent.RetainedParameters or PointwiseContent.CollapsedPredictive))
            throw new ArgumentOutOfRangeException(nameof(content));
        if (content == PointwiseContent.RetainedParameters && model.Observation.Law is { DensityByQuadrature: true } && !quadrature)
            throw new NotSupportedException("This noise's density is an integral over each observation's latent, N × draws integrals of hundreds "
                                            + "of nodes each; pass quadrature: true to take it.");
        UnitPartition partition = units ?? DefaultUnits(model);
        if (partition.Observations != model.Data.Response.Observations)
            throw new ArgumentException("The partition must cover the model's N observations, case by case and location by location.", nameof(units));
        bool collapsed = content == PointwiseContent.CollapsedPredictive;
        if (collapsed)
        {
            if (model.Observation.Law is not null)
                throw new NotSupportedException("The collapsed predictive integrates a Gaussian noise; under a noise law it waits for a check of its own. Take the retained content, which integrates each observation's latent.");
            if (model.RandomIntercept is not null)
                throw new NotSupportedException("The collapsed predictive under random intercepts would leave a unit out of the block as well, whose whitening couples every case of its level; not built (E21). Take the retained content, which reads each draw's group effects.");
            if (model.Working is { IsUnweighted: false, Meaning: WeightMeaning.Replication })
                throw new NotSupportedException("The collapsed predictive leaves a row out whole; under weights that count rows it would leave out all its copies at once. Take the retained content.");
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

        List<LatentReport> latents = Latents(model, pooled, collapsed);
        return weights is null
            ? new PointwiseLogDensity(values, partition, latents, chains)
            : new PointwiseLogDensity(values, partition, latents) { Weights = weights };
    }

    /// <summary>
    /// How the terms treat each latent, for <paramref name="pooled"/>'s draws: the structure, τ, σ²'s rate, the coefficients and σ² global
    /// and conditioned on, or the last two integrated given the rest for the collapsed content; a noise law's latent per observation
    /// integrated out, its prior being independent of every other observation's given the globals; a moved working response's offset and
    /// weights conditioned on. The candidate density (<see cref="BassCandidates"/>) reports the same, since every unit is read by one
    /// density (CB-4).
    /// </summary>
    internal static List<LatentReport> Latents(BassModel model, List<BassDraw> pooled, bool collapsed)
    {
        LatentTreatment treatment = collapsed ? LatentTreatment.IntegratedGivenOthers : LatentTreatment.Conditioned;
        var latents = new List<LatentReport> { new("structure", LatentScope.Global, LatentTreatment.Conditioned) };
        if (model.Target.SamplesTau) latents.Add(new LatentReport("prior precision τ", LatentScope.Global, LatentTreatment.Conditioned));
        if (model.Target.Sigma2Rate is not null) latents.Add(new LatentReport("σ²'s rate", LatentScope.Global, LatentTreatment.Conditioned));
        if (model.RandomIntercept is not null)   // the whitening is the block's integral, not an offset conditioned on
            latents.Add(new LatentReport("the group effects u_g", LatentScope.Group, LatentTreatment.Conditioned));
        else if (model.Observation.Law is not null)
            latents.Add(new LatentReport("the noise's latent per observation", LatentScope.Observation, LatentTreatment.Integrated));
        else if (pooled.Exists(d => !ReferenceEquals(d.WorkingOn(model), model.Working)))
            latents.Add(new LatentReport("the working response's offset and weights", LatentScope.Observation, LatentTreatment.Conditioned));
        latents.Add(new LatentReport("coefficients a", LatentScope.Global, treatment));
        latents.Add(new LatentReport("noise variance σ²", LatentScope.Global, treatment));
        return latents;
    }

    // Each draw's Normal log terms at its fit and σ², summed over each unit's observations: on the working response it was composed on,
    // under W at its meaning, N(y | f, σ²/w) or w copies of N(y | f, σ²). Under a noise law, the law's density of y − f, each
    // observation's latent integrated out.
    private static double[,] Retained(BassModel model, BassDrawFit fit, List<BassDraw> pooled, UnitPartition partition)
    {
        int n = model.Data.Response.Observations;
        var values = new double[pooled.Count, partition.Units];
        var f = new double[n];
        var log = new double[n];
        IBassNoiseLaw? law = model.Observation.Law;
        for (int s = 0; s < pooled.Count; s++)
        {
            BassDraw draw = pooled[s];
            if (law is not null)
            {
                ReadOnlySpan<double> observed = model.Data.Response.Values;
                ObservationLatents? globals = draw.LatentsOn(model);
                fit.Evaluate(draw, f);
                for (int i = 0; i < n; i++) log[i] = law.LogDensity(observed[i] - f[i], draw.Sigma2, globals);
                Sum(log, partition, values, s);
                continue;
            }
            WorkingResponse working = draw.WorkingOn(model);
            // Under random intercepts the working response is the whitened one, and the fit reads f + u_g at the observations (E21).
            ReadOnlySpan<double> y = model.RandomIntercept is null ? working.Values : model.Data.Response.Values;
            double[]? weights = working.WeightArray;
            bool counts = working.Meaning == WeightMeaning.Replication;
            fit.Evaluate(draw, f);
            double sigma2 = draw.Sigma2, normalizer = -0.5 * Math.Log(2.0 * Math.PI * sigma2);
            for (int i = 0; i < n; i++)
            {
                double e = y[i] - f[i];
                if (weights is null) log[i] = normalizer - 0.5 * e * e / sigma2;
                else if (counts) log[i] = weights[i] * (normalizer - 0.5 * e * e / sigma2);
                else log[i] = normalizer + 0.5 * Math.Log(weights[i]) - 0.5 * weights[i] * e * e / sigma2;
            }
            Sum(log, partition, values, s);
        }
        return values;
    }

    // Draw s's log terms summed over each unit's observations.
    private static void Sum(double[] log, UnitPartition partition, double[,] values, int s)
    {
        for (int u = 0; u < partition.Units; u++)
        {
            double sum = 0.0;
            foreach (int i in partition.Members[u]) sum += log[i];
            values[s, u] = sum;
        }
    }

    // The draws grouped by structure and working response: one design, block and explicit design per group, then each draw's
    // leave-unit-out predictive at its τ.
    private static double[,] Collapsed(BassModel model, BassDrawFit fit, List<BassDraw> pooled, UnitPartition partition, string parameter)
    {
        BassTarget target = model.Target;
        var groups = new Dictionary<(int[] Ids, WorkingResponse Working), List<int>>(BassDrawGroupComparer.Instance);
        for (int s = 0; s < pooled.Count; s++)
        {
            (int[], WorkingResponse) key = (pooled[s].Structure.TermIds, pooled[s].WorkingOn(model));
            if (!groups.TryGetValue(key, out List<int>? rows)) groups.Add(key, rows = []);
            rows.Add(s);
        }

        var values = new double[pooled.Count, partition.Units];
        foreach (((int[] ids, WorkingResponse working), List<int> rows) in groups)
        {
            ReadOnlySpan<double> y = working.Values;
            double[]? weights = working.WeightArray;
            var structure = new BassTerm[ids.Length];
            for (int j = 0; j < ids.Length; j++) structure[j] = fit.Term(ids[j]);
            BassDesign design = model.NewDesign(working);
            BassRefusal? refusal = design.Load(structure) ?? target.Refusal(BassStructureFit.Of(design, 0.0));
            if (refusal is { } reason) throw new ArgumentException($"A draw's structure lies outside F ({reason}).", parameter);
            Func<double, double, GPriorConditional> block = target.CoefficientBlock(design);
            bool perTau = design.PriorTerms is not null;
            double[,] x = fit.CentredDesign(ids, working);
            // The g-form's units' quadratic forms are τ-free, taken once per structure, those of a prior apart from its Gram once
            // per τ; a draw at a τ already seen copies its row.
            LeaveOutQuadratics[]? quadratics = null;
            var seen = new Dictionary<(double Tau, double? Rate), int>();
            foreach (int s in rows)
            {
                if (seen.TryGetValue((pooled[s].Tau, pooled[s].Sigma2Rate), out int first))
                {
                    for (int u = 0; u < partition.Units; u++) values[s, u] = values[first, u];
                    continue;
                }
                seen.Add((pooled[s].Tau, pooled[s].Sigma2Rate), s);
                InverseGamma? sigma2 = pooled[s].Sigma2Rate is double rate ? InverseGamma.ShapeScale(target.Sigma2Rate!.Shape, rate) : target.Sigma2;
                GPriorConditional conditional = block(pooled[s].Tau, 1.0);
                if (quadratics is null || perTau)
                {
                    quadratics = new LeaveOutQuadratics[partition.Units];
                    for (int u = 0; u < partition.Units; u++) quadratics[u] = conditional.LeaveOutQuadratics(x, partition.Members[u]);
                }
                for (int u = 0; u < partition.Units; u++)
                    values[s, u] = target.GivenSigma2 is double given
                        ? conditional.LeaveOutLogDensityGiven(given, x, y, weights, quadratics[u])
                        : conditional.LeaveOutLogDensity(sigma2, x, y, weights, quadratics[u]);
            }
        }
        return values;
    }
}

/// <summary>Draws that share a block: the same structure's ids on the same working response.</summary>
internal sealed class BassDrawGroupComparer : IEqualityComparer<(int[] Ids, WorkingResponse Working)>
{
    public static readonly BassDrawGroupComparer Instance = new();

    public bool Equals((int[] Ids, WorkingResponse Working) a, (int[] Ids, WorkingResponse Working) b)
        => ReferenceEquals(a.Working, b.Working) && a.Ids.AsSpan().SequenceEqual(b.Ids);

    public int GetHashCode((int[] Ids, WorkingResponse Working) key)
    {
        var hash = new HashCode();
        hash.Add(key.Working);
        foreach (int id in key.Ids) hash.Add(id);
        return hash.ToHashCode();
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
        _block = model.RandomIntercept;
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
        if (draw.OutputCoefficients is not null)
            throw new NotSupportedException("A draw of several outputs is read by E34's readouts (its slice 3); the readers of one response refuse it.");
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

    private readonly RandomInterceptBlock? _block;

    /// <summary>The draw's fit at every observation, into <paramref name="fit"/>: f, and under random intercepts f + u_g, each case reading
    /// its level's effect (E21).</summary>
    public void Evaluate(BassDraw draw, Span<double> fit)
    {
        EvaluateTerms(draw, fit);
        if (_block is null) return;
        double[] u = draw.GroupEffects ?? throw new ArgumentException("A random-intercept fit's draw carries its group effects.", nameof(draw));
        for (int i = 0; i < _cases; i++) fit[i] += u[_block.GroupOf(i)];
    }

    private void EvaluateTerms(BassDraw draw, Span<double> fit)
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
    /// The posterior mean of the offsets of the working responses <paramref name="read"/> gives the draws, at every observation under
    /// <paramref name="weights"/> (null weighs them equally), a draw with none counting as zero; null when no draw reads one.
    /// </summary>
    public static double[]? MeanOffset(IReadOnlyList<BassDraw> pooled, IReadOnlyList<double>? weights, int observations, Func<BassDraw, WorkingResponse?> read)
    {
        double[]? sum = null;
        double total = 0.0;
        for (int s = 0; s < pooled.Count; s++)
        {
            double weight = weights?[s] ?? 1.0;
            total += weight;
            ReadOnlySpan<double> offset = read(pooled[s]) is { } working ? working.Offset : default;
            if (offset.IsEmpty) continue;
            sum ??= new double[observations];
            for (int i = 0; i < observations; i++) sum[i] += weight * offset[i];
        }
        if (sum is null) return null;
        for (int i = 0; i < observations; i++) sum[i] /= total;
        return sum;
    }

    /// <summary>
    /// The explicit design [1, Z_c] of the structure <paramref name="ids"/>, in their order: the intercept's column of
    /// ones, then each term's column centred over all N observations under W, a − ā for a scalar response and
    /// a(xᵢ)h(t_r) − level for a direct functional one, as the design layer centres them. W is <paramref name="working"/>'s, the model's
    /// when null.
    /// </summary>
    public double[,] CentredDesign(int[] ids, WorkingResponse? working = null)
    {
        WorkingResponse centring = working ?? _model.Working;
        int n = _cases * _locations;
        var x = new double[n, ids.Length + 1];
        for (int i = 0; i < n; i++) x[i, 0] = 1.0;
        for (int j = 0; j < ids.Length; j++)
        {
            (double[] a, double[]? h) = Columns(ids[j]);
            double level = centring.ColumnMean(a, h);
            if (h is null)
            {
                for (int i = 0; i < _cases; i++) x[i, j + 1] = a[i] - level;
                continue;
            }
            for (int i = 0; i < _cases; i++)
                for (int r = 0; r < _locations; r++) x[i * _locations + r, j + 1] = a[i] * h[r] - level;
        }
        return x;
    }
}
