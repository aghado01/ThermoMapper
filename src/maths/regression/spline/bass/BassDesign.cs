using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline.Bass;

/// <summary>Why a structure lies outside the support F (D02's global truncation; D06's refusal reasons).</summary>
public enum BassRefusal
{
    /// <summary>M would exceed the basis cap.</summary>
    Cap,

    /// <summary>A term falls short of the support threshold in a part (D01).</summary>
    Support,

    /// <summary>A column's angle to the span of the others falls below s_min, a repeated term included (D06).</summary>
    Numerical,

    /// <summary>A column's variance inflation exceeds the declared cap, a statistical restriction (D06's option).</summary>
    VarianceInflation,

    /// <summary>Under the Jeffreys σ² prior, the structure fits the response exactly, which would make the posterior
    /// improper (D03's propriety rule); the collapsed target refuses it.</summary>
    ExactFit,
}

/// <summary>A design's restart state (<see cref="BassDesign.Snapshot"/>): its terms in column order, a copy of its
/// factor, and its coefficients, residual and statistics; for an in-memory hand-over, also its evaluated columns, shared
/// (null in a checkpoint, which restores them from the terms); and under the ridge its spectrum when it has one, or under a
/// prior apart from its Gram its pencil, immutable values the restored design shares.</summary>
internal sealed record BassDesignSnapshot(BassTerm[] Terms, BorderedCholesky Factor, double[] Beta, double[] Residual,
                                          double ResidualSquares, double FittedSquares, BassDesignColumns? Columns = null,
                                          RidgeSpectrum? Spectrum = null, PriorPencil? Pencil = null);

/// <summary>A design's evaluated columns in column order: each term's case side, its output side for a direct response held
/// as the Khatri–Rao pair, and its cross product with the response; under a prior apart from its Gram, each term's prior
/// side too. A design never writes into a column once evaluated, so designs over the same data can share them.</summary>
internal sealed record BassDesignColumns(double[][] Cases, double[][]? Outputs, double[] CrossResponses, BassPriorColumn[]? Priors = null);

/// <summary>A column on the prior's side: centred over N rows free of the weights, with δ = z̄ − z̄_w, the two centres' gap,
/// which couples the intercept to the column in the coordinates the design holds.</summary>
internal sealed record BassPriorColumn(double[] Centred, double Shift);

/// <summary>The admission rules a design applies, resolved from the config on the data.</summary>
public sealed record BassAdmission(
    int MaxBasis, SupportKind Support, double SupportMinimum, double? SupportMinimumFunctional, double AngleFloor,
    double? MaxVarianceInflation)
{
    /// <summary>The rules <paramref name="config"/> sets, with the data's thresholds from <paramref name="resolved"/>.</summary>
    public static BassAdmission From(BassConfig config, BassResolved resolved)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(resolved);
        return new(resolved.MaxBasis, config.Structure.Support.Kind, resolved.SupportMinimum, resolved.SupportMinimumFunctional,
                   config.Structure.AngleFloor, config.Structure.MaxVarianceInflation);
    }
}

/// <summary>
/// A proposed structure's account: admitted, or refused with its reason; and, when admitted, its column count, its
/// residual and fitted sums of squares in the centred basis, its smallest leave-one-out angle (D06's gate), the new
/// column's angle to the span of the others (NaN for a death), and under a prior the g-form's closed form does not serve
/// its columns' terms at any (τ, b) (<see cref="RidgeColumns"/>; null under the g-form on the design's own Gram).
/// </summary>
public readonly record struct DesignProposal(
    BassRefusal? Refusal, int Columns, double ResidualSquares, double FittedSquares, double MinAngle, double NewestAngle,
    RidgeColumns? Ridge = null)
{
    /// <summary>Whether the structure lies in F.</summary>
    public bool Admitted => Refusal is null;

    internal static DesignProposal Refused(BassRefusal reason)
        => new(reason, -1, double.NaN, double.NaN, double.NaN, double.NaN);
}

/// <summary>
/// One chain's design: its terms' columns in the centred basis, the factor of their Gram, and the least-squares fit of
/// the centred working response (the response less the observation model's offset, <see cref="WorkingResponse"/>),
/// maintained through births, deaths and changes (stage 4's S4; D06; the numerics segment's dense factor maintenance).
/// <para><b>The centred basis.</b> The marginal depends only on the column space of [1, Z] under every intercept
/// convention, so the design holds Z_c, the columns centred over all N rows, and ỹ = y − ȳ: the intercept is
/// orthogonal to them, and its Gram entry N stands apart (D03). Its statistics are the fit's residual ‖r‖² and fitted
/// sum Q_c = b_cᵀG_c⁻¹b_c, with b_c = Z_cᵀỹ; a convention adds its own pieces (the zero centre's Nȳ²) on top.
/// Scalar fits hold explicitly centred columns; direct functional fits hold the case-side and output-side factors,
/// and their Gram entries come from the centred Hadamard form (<see cref="KhatriRao"/>), a sum of three positive
/// semidefinite terms. Either way every Gram entry is formed without cancellation, and every inner product is a
/// <see cref="PairwiseSum"/>.</para>
/// <para><b>Admission</b> (D06), checked in this order on a birth or a change: the cap; a repeated term, whose column
/// duplicates one already there (s = 0); the support threshold in each part, by count or by mass; the order-free
/// angle gate, every column's leave-one-out angle s_j = 1/(G_jj (G⁻¹)_jj) at least s_min; and the declared variance
/// inflation cap. A death needs no gate, since the gate is downward closed. When the bordered factor fails the gate or
/// comes within a factor of two of s_min, a strict factorization of the candidate's Gram decides, and an admitted
/// candidate is then accepted from that factorization.</para>
/// <para><b>Statistics.</b> A candidate's ‖r′‖² updates the current residual's: a birth's removes (zᵀr)²/δ, with δ
/// the new column's pivot; a death's adds β̂_m²/(G⁻¹)_mm; a change does both, through the reduced model. These lose
/// digits only when one move explains nearly all the residual, where the candidate's residual is formed explicitly
/// instead. Q′ is the candidate's quadratic form. On every acceptance β̂ is solved and the residual recomputed
/// explicitly, r = ỹ − Z_cβ̂, so it never drifts; the factor rebuilds strictly on its cadence.</para>
/// <para><b>Under the normalized ridge</b> (D03, decided 2026-10-29) the design also keeps its current structure's
/// <see cref="RidgeSpectrum"/>, decomposed once per structure when first read, and each admitted proposal carries its
/// columns' ridge terms: a bordered candidate of the current spectrum, O(M²) to start and O(M) per (τ, b), or, where the
/// strict factorization decided admission near s_min, the candidate's own fresh decomposition, which an acceptance then
/// keeps. Admission is the same under either form, so F is unchanged.</para>
/// <para><b>Under weights</b> (<see cref="WorkingResponse.Weights"/>) the rows are whitened: a column is held as
/// W^½(z − z̄_w), centred under W, and the response as W^½(y − ȳ_w). Every inner product above then reads the weighted form,
/// so bordering, rebuilds, the gate and the statistics carry the weights with no branch of their own: G_c = Z_cᵀWZ_c,
/// b_c = Z_cᵀWỹ, ‖r‖² = ‖r‖²_W, and the intercept's Gram entry is Σwᵢ. A direct functional design holds its stacked columns,
/// since per-row weights do not contract through the Khatri–Rao pair. The support threshold reads the raw parts, or a scalar
/// response's rows as the weights count them under <see cref="WeightMeaning.Replication"/>.</para>
/// <para><b>With the prior free of the weights</b> (<see cref="WeightedPriorGram.Unweighted"/>, D03's robust row) the design also
/// holds each column on the prior's side and carries a <see cref="PriorPencil"/> for each admitted proposal, in the coordinates
/// (α, β_c) with the columns centred under W. The prior's own Gram must clear the angle gate too, since it is a precision; so
/// F is unchanged by the weights when they are unit.</para>
/// Not thread-safe: one per chain.
/// </summary>
public sealed class BassDesign
{
    private const int InitialCapacity = 32;
    private const double ExplicitBelow = 1e-8;   // ‖r′‖²/‖r_reduced‖² under which the candidate's residual is formed

    private enum Move { Birth, Death, Change }

    private sealed record Pending(Move Kind, int Index, BassTerm? Term, double[]? Column, double[]? Output, double BEntry,
                                  double[,]? Strict, DesignProposal Account, RidgeSpectrum? Fresh = null, BassPriorColumn? Prior = null,
                                  double[,]? PriorGram = null, PriorPencil? Pencil = null);

    private readonly BassBasis _basis;
    private readonly BassAdmission _rules;
    private readonly bool _functional, _ridge;
    private readonly bool _khatriRao;               // a direct response held as case and output factors, its rows unweighted
    private readonly int _cases, _locations, _rows;
    private readonly double[]? _weights, _roots;    // W and W^½; null at unit weights
    private readonly double[]? _counted;            // a scalar response's weights, when they count rows for the support threshold
    private readonly GPriorCentre _intercept;       // what a prior apart from the Gram does with the intercept
    private readonly List<BassPriorColumn>? _priorColumns;   // the prior's side, when it stands apart from the Gram
    private double[,]? _priorGram;                  // G_c on the prior's side
    private PriorPencil? _pencil;                   // the current structure's; null until read
    private readonly IInnerProduct _output;         // Q on the output index; the identity for independent noise
    private readonly double[] _y;                   // ỹ, case by case and location by location
    private readonly double[] _residual;
    private readonly List<BassTerm> _terms = [];
    private readonly List<double[]> _a = [];        // scalar: centred columns; direct: case-side factors
    private readonly List<double[]> _h = [];        // direct: output-side factors
    private readonly List<double> _b = [];          // Z_cᵀỹ
    private double[] _beta = [];
    private BorderedCholesky _factor;
    private BorderedCholesky? _strict;
    private Pending? _pending;
    private RidgeSpectrum? _spectrum;   // the current structure's, under the ridge; null until read

    /// <summary>The intercept-only design over <paramref name="basis"/>'s inputs for <paramref name="working"/>, the
    /// observation model's working response for <paramref name="data"/>'s (the response itself when null);
    /// <paramref name="ridge"/> keeps the normalized ridge's spectrum and proposal terms beside the factor; and
    /// <paramref name="intercept"/> places the intercept for a prior apart from the Gram, which alone reads it.</summary>
    public BassDesign(BassBasis basis, BassData data, BassAdmission rules, bool ridge = false, WorkingResponse? working = null,
                      GPriorCentre intercept = GPriorCentre.Mean)
    {
        ArgumentNullException.ThrowIfNull(basis);
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(rules);
        _basis = basis;
        _rules = rules;
        _ridge = ridge;
        _functional = basis.Inputs.Outputs is not null;
        if (_functional != data.Response is FunctionalResponse)
            throw new ArgumentException("The basis and the data disagree on the response's kind.", nameof(data));
        _cases = data.Response.Cases;
        _locations = data.Response.Locations;
        _rows = data.Response.Observations;
        if (basis.Inputs.Cases.Cases != _cases) throw new ArgumentException("The basis was not prepared from these data.", nameof(basis));
        working ??= new WorkingResponse(data.Response);
        if (!ReferenceEquals(working.Response, data.Response)) throw new ArgumentException("The working response is of another response.", nameof(working));
        _intercept = intercept;
        _output = working.Output;
        if (!_output.IsIdentity)
            throw new NotSupportedException("The core holds the independent Gram; a Q other than the identity arrives with the extension that centres under it (F-X6).");
        _weights = working.WeightArray;
        _khatriRao = _functional && _weights is null;
        ReadOnlySpan<double> y = working.Values;
        ResponseMean = working.ResponseMean;
        Weighting = working.Weighting;
        _y = new double[_rows];
        for (int i = 0; i < _rows; i++) _y[i] = y[i] - ResponseMean;
        if (_weights is not null)
        {
            _roots = new double[_rows];
            for (int i = 0; i < _rows; i++)
            {
                _roots[i] = Math.Sqrt(_weights[i]);
                _y[i] *= _roots[i];
            }
            if (!_functional && working.Meaning == WeightMeaning.Replication) _counted = _weights;
            if (working.PriorStandsApart)
            {
                if (intercept is not (GPriorCentre.Mean or GPriorCentre.Zero or GPriorCentre.Flat)) throw new ArgumentOutOfRangeException(nameof(intercept));
                _priorColumns = [];
                _priorGram = new double[0, 0];
                PriorLevel = intercept == GPriorCentre.Zero ? 0.0 : PairwiseSum.Of(y) / _rows;
            }
        }
        ResponseSquares = _output.Dot(_y, _y);
        _residual = (double[])_y.Clone();
        ResidualSquares = ResponseSquares;
        _factor = new BorderedCholesky(Math.Min(InitialCapacity, Math.Max(rules.MaxBasis, 1)));
    }

    /// <summary>M, the terms.</summary>
    public int Count => _terms.Count;

    /// <summary>The terms, in the factor's column order.</summary>
    public IReadOnlyList<BassTerm> Terms => _terms;

    /// <summary>N, the observations.</summary>
    public int Observations => _rows;

    /// <summary>The count, the mass and Σ ln wᵢ the marginal reads beyond these statistics; (N, N, 0) at unit weights.</summary>
    public GPriorWeights Weighting { get; }

    /// <summary>ȳ, the response's mean over all N observations under W.</summary>
    public double ResponseMean { get; }

    /// <summary>The level a prior apart from the Gram centres the intercept on, free of the weights; 0 otherwise.</summary>
    public double PriorLevel { get; }

    /// <summary>‖ỹ‖²_W = ‖y − ȳ‖²_W.</summary>
    public double ResponseSquares { get; }

    /// <summary>‖r‖²_W, the fit's residual sum of squares, from the explicit residual.</summary>
    public double ResidualSquares { get; private set; }

    /// <summary>Q_c = b_cᵀG_c⁻¹b_c, the fitted sum of squares in the centred basis.</summary>
    public double FittedSquares { get; private set; }

    /// <summary>β̂_c, the least-squares coefficients of the centred columns.</summary>
    public ReadOnlySpan<double> Coefficients => _beta;

    /// <summary>The explicit residual r = ỹ − Z_cβ̂, whitened (W^½r) under weights, case by case and location by location.</summary>
    public ReadOnlySpan<double> Residual => _residual;

    /// <summary>The factor of G_c, in the terms' order.</summary>
    public ISpdFactor Factor => _factor;

    /// <summary>Whether the coefficient prior's form is the normalized ridge.</summary>
    public bool IsRidge => _ridge;

    /// <summary>The current structure's terms for a prior the g-form's closed form does not serve: the normalized ridge's
    /// <see cref="Spectrum"/>, or the <see cref="Pencil"/> of a prior apart from the Gram; null otherwise.</summary>
    public RidgeColumns? PriorTerms => _priorColumns is not null ? Pencil : _ridge ? Spectrum : null;

    /// <summary>The current structure's pencil, under weights with the prior free of them; null otherwise.</summary>
    public PriorPencil? Pencil
    {
        get
        {
            if (_priorColumns is null) return null;
            if (_pencil is { } cached) return cached;
            return _pencil = PencilOf(CurrentGram(), _priorGram!, CurrentShifts(), [.. _b])
                ?? throw new InvalidOperationException("The current structure's prior is not positive definite.");
        }
    }

    /// <summary>The current structure's spectrum under the normalized ridge, decomposed from the Gram when first read after
    /// the structure changed; a ridge design only.</summary>
    public RidgeSpectrum Spectrum
    {
        get
        {
            if (!_ridge) throw new InvalidOperationException("The design holds the g-form; it keeps no ridge spectrum.");
            if (_priorColumns is not null) throw new InvalidOperationException("The prior stands apart from the Gram; the design keeps a pencil, not a spectrum.");
            if (_spectrum is { } cached) return cached;
            int n = _terms.Count;
            var g = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) g[i, j] = _factor.Entry(i, j);
            return _spectrum = RidgeSpectrum.Of(g, [.. _b], _rows, Weighting.Mass);
        }
    }

    /// <summary>The smallest leave-one-out angle of the current columns; +∞ with none.</summary>
    public double MinAngle => _factor.MinLeaveOneOutPivot();

    /// <summary>Term <paramref name="index"/>'s leave-one-out angle s_j: its column's squared sine to the span of the others
    /// and the intercept, the reciprocal of its centred variance inflation factor (D06).</summary>
    public double Angle(int index) => _factor.LeaveOneOutPivot(index);

    /// <summary>Proposes adding <paramref name="term"/>.</summary>
    public DesignProposal ProposeBirth(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        _pending = null;
        if (_terms.Count >= _rules.MaxBasis) return DesignProposal.Refused(BassRefusal.Cap);
        if (_terms.Contains(term)) return DesignProposal.Refused(BassRefusal.Numerical);
        if (Evaluate(term) is not var (column, output, prior)) return DesignProposal.Refused(BassRefusal.Support);
        EnsureCapacity(_terms.Count + 1);
        double[] border = Border(column, output, -1, out double diagonal);
        double bEntry = CrossResponse(column, output);
        FactorProposal p = _factor.ProposeBirth(border, diagonal);
        return Decide(Move.Birth, _terms.Count, term, column, output, border, border, diagonal, bEntry, p,
                      ResidualSquares, CrossResidual(column, output), prior);
    }

    /// <summary>Proposes removing term <paramref name="index"/>. A death never leaves F, since the gate is downward closed.</summary>
    public DesignProposal ProposeDeath(int index)
    {
        if ((uint)index >= (uint)_terms.Count) throw new ArgumentOutOfRangeException(nameof(index));
        _pending = null;
        FactorProposal p = _factor.ProposeDeath(index);
        double[] b = Without(_b, index);
        double q = _factor.ProposedInverseQuadraticForm(b);
        double rss = ResidualSquares + _beta[index] * _beta[index] / _factor.InverseDiagonal(index);
        if (_priorColumns is not null)
        {
            // Downward closed: a principal block of a positive definite prior Gram stays so.
            double[,] priorGram = PriorGramWithout(index);
            PriorPencil pencil = PencilOf(Without(CurrentGram(), index), priorGram, Without(CurrentShifts(), index), b)!;
            var apart = new DesignProposal(null, _terms.Count - 1, rss, q, p.MinLeaveOneOutPivot, double.NaN, pencil);
            _pending = new Pending(Move.Death, index, null, null, null, 0.0, null, apart, null, null, priorGram, pencil);
            return apart;
        }
        var account = new DesignProposal(null, _terms.Count - 1, rss, q, p.MinLeaveOneOutPivot, double.NaN, _ridge ? Spectrum.Death(index) : null);
        _pending = new Pending(Move.Death, index, null, null, null, 0.0, null, account);
        return account;
    }

    /// <summary>Proposes replacing term <paramref name="index"/> with <paramref name="term"/>, which goes last.</summary>
    public DesignProposal ProposeChange(int index, BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        if ((uint)index >= (uint)_terms.Count) throw new ArgumentOutOfRangeException(nameof(index));
        _pending = null;
        if (term.Equals(_terms[index]))   // the same state
        {
            var same = new DesignProposal(null, _terms.Count, ResidualSquares, FittedSquares, MinAngle, _factor.LeaveOneOutPivot(index), PriorTerms);
            _pending = new Pending(Move.Change, index, term, null, null, 0.0, null, same);
            return same;
        }
        if (_terms.Contains(term)) return DesignProposal.Refused(BassRefusal.Numerical);
        if (Evaluate(term) is not var (column, output, prior)) return DesignProposal.Refused(BassRefusal.Support);
        double[] full = Border(column, output, -1, out double diagonal);
        double bEntry = CrossResponse(column, output);

        // The reduced model, term index removed: ‖r‖² gains β̂_m²/(G⁻¹)_mm, and zᵀr gains β̂_m·uᵀG⁻¹e_m/(G⁻¹)_mm.
        var unit = new double[_terms.Count];
        unit[index] = 1.0;
        _factor.SolveInPlace(unit);
        double inverseDiagonal = unit[index], beta = _beta[index];
        double rssReduced = ResidualSquares + beta * beta / inverseDiagonal;
        double zrReduced = CrossResidual(column, output) + beta * PairwiseSum.Dot(full, unit) / inverseDiagonal;

        double[] border = Without(full, index);
        FactorProposal p = _factor.ProposeChange(index, border, diagonal);
        return Decide(Move.Change, index, term, column, output, border, full, diagonal, bEntry, p, rssReduced, zrReduced, prior);
    }

    /// <summary>Commits the last proposal, which must have been admitted.</summary>
    public void Accept()
    {
        Pending p = _pending ?? throw new InvalidOperationException("Nothing admitted is proposed.");
        _pending = null;
        if (p.Kind == Move.Change && p.Column is null) return;   // a change to the same term
        if (p.Strict is not null) _factor.Decompose(p.Strict);
        else _factor.Accept();
        _spectrum = p.Fresh;
        if (p.Kind is Move.Death or Move.Change)
        {
            _terms.RemoveAt(p.Index);
            _a.RemoveAt(p.Index);
            if (_khatriRao) _h.RemoveAt(p.Index);
            _b.RemoveAt(p.Index);
            _priorColumns?.RemoveAt(p.Index);
        }
        if (p.Kind is Move.Birth or Move.Change)
        {
            _terms.Add(p.Term!);
            _a.Add(p.Column!);
            if (_khatriRao) _h.Add(p.Output!);
            _b.Add(p.BEntry);
            if (p.Prior is { } prior) _priorColumns!.Add(prior);
        }
        if (_priorColumns is not null)
        {
            _priorGram = p.PriorGram;
            _pencil = p.Pencil;
        }
        Refit();
    }

    /// <summary>
    /// Replaces the structure with <paramref name="terms"/>, built and factored strictly from scratch, when it lies in
    /// F; otherwise returns the refusal and leaves the design as it was. Starts and restores use it.
    /// </summary>
    public BassRefusal? Load(IReadOnlyList<BassTerm> terms)
    {
        ArgumentNullException.ThrowIfNull(terms);
        _pending = null;
        if (terms.Count > _rules.MaxBasis) return BassRefusal.Cap;
        var seen = new HashSet<BassTerm>();
        var columns = new List<double[]>(terms.Count);
        var outputs = new List<double[]>(terms.Count);
        var priors = new List<BassPriorColumn>(terms.Count);
        foreach (BassTerm term in terms)
        {
            if (!seen.Add(term)) return BassRefusal.Numerical;
            if (Evaluate(term) is not var (column, output, prior)) return BassRefusal.Support;
            columns.Add(column);
            if (_khatriRao) outputs.Add(output!);
            if (prior is not null) priors.Add(prior);
        }
        double[,] g = Gram(columns, outputs);
        var fresh = new BorderedCholesky(Math.Max(_factor.Capacity, terms.Count));
        fresh.Decompose(g);
        if (terms.Count > 0 && fresh.FailedAt >= 0) return BassRefusal.Numerical;
        double minAngle = fresh.MinLeaveOneOutPivot();
        if (minAngle < _rules.AngleFloor) return BassRefusal.Numerical;
        if (_rules.MaxVarianceInflation is { } cap && minAngle * cap < 1.0) return BassRefusal.VarianceInflation;
        var b = new double[columns.Count];
        for (int j = 0; j < columns.Count; j++) b[j] = CrossResponse(columns[j], _khatriRao ? outputs[j] : null);
        double[,]? priorGram = null;
        PriorPencil? pencil = null;
        if (_priorColumns is not null)
        {
            priorGram = PriorGramOf(priors);
            double[] shifts = [.. priors.Select(c => c.Shift)];
            if (PencilOf(g, priorGram, shifts, b) is not { } admitted) return BassRefusal.Numerical;
            pencil = admitted;
        }

        _factor = fresh;
        _spectrum = null;
        _terms.Clear(); _terms.AddRange(terms);
        _a.Clear(); _a.AddRange(columns);
        _h.Clear(); _h.AddRange(outputs);
        _b.Clear(); _b.AddRange(b);
        if (_priorColumns is not null)
        {
            _priorColumns.Clear(); _priorColumns.AddRange(priors);
            _priorGram = priorGram;
            _pencil = pencil;
        }
        Refit();
        return null;
    }

    /// <summary>
    /// The restart state: the terms in column order, an exact copy of the factor, and the coefficients, residual and its
    /// statistics as held. A strict refactorization would not reproduce the factor's incremental updates bit for bit, so
    /// <see cref="Restore"/> takes them as they are. The columns are functions of the terms: a checkpoint's snapshot leaves
    /// them out and <see cref="Restore"/> re-evaluates them, while a hand-over between replicas
    /// (<paramref name="shareColumns"/>) carries this design's own, which the receiver shares.
    /// </summary>
    internal BassDesignSnapshot Snapshot(bool shareColumns = false)
        => new([.. _terms], _factor.Copy(), (double[])_beta.Clone(), (double[])_residual.Clone(), ResidualSquares, FittedSquares,
               shareColumns ? new BassDesignColumns([.. _a], _khatriRao ? [.. _h] : null, [.. _b], _priorColumns is null ? null : [.. _priorColumns]) : null,
               _spectrum, _pencil);

    /// <summary>Replaces the structure with <paramref name="snapshot"/>'s, every statistic as the design it was taken from
    /// held it, so the two continue bit for bit; its columns shared when it carries them, else re-evaluated.</summary>
    internal void Restore(BassDesignSnapshot snapshot)
    {
        ArgumentNullException.ThrowIfNull(snapshot);
        if (snapshot.Factor.Dimension != snapshot.Terms.Length || snapshot.Beta.Length != snapshot.Terms.Length || snapshot.Residual.Length != _rows)
            throw new ArgumentException("The snapshot is not a design's over these data.", nameof(snapshot));
        if (snapshot.Columns is { } shared
            && (shared.Cases.Length != snapshot.Terms.Length || shared.CrossResponses.Length != snapshot.Terms.Length
                || (shared.Outputs is null) == _khatriRao || shared.Outputs is { } o && o.Length != snapshot.Terms.Length
                || (shared.Priors is null) == (_priorColumns is not null) || shared.Priors is { } q && q.Length != snapshot.Terms.Length
                || Array.Exists(shared.Cases, c => c.Length != (_khatriRao ? _cases : _rows))
                || shared.Outputs is { } h && Array.Exists(h, c => c.Length != _locations)))
            throw new ArgumentException("The snapshot's columns are not its terms' over these data.", nameof(snapshot));
        _pending = null;
        _terms.Clear();
        _a.Clear();
        _h.Clear();
        _b.Clear();
        _priorColumns?.Clear();
        if (snapshot.Columns is { } columns)
        {
            _terms.AddRange(snapshot.Terms);
            _a.AddRange(columns.Cases);
            if (_khatriRao) _h.AddRange(columns.Outputs!);
            _b.AddRange(columns.CrossResponses);
            if (columns.Priors is { } held) _priorColumns!.AddRange(held);
        }
        else foreach (BassTerm term in snapshot.Terms)
        {
            if (Evaluate(term) is not var (column, output, prior)) throw new ArgumentException($"The snapshot holds {term}, which lacks support here.", nameof(snapshot));
            _terms.Add(term);
            _a.Add(column);
            if (_khatriRao) _h.Add(output!);
            _b.Add(CrossResponse(column, output));
            if (prior is not null) _priorColumns!.Add(prior);
        }
        _factor = snapshot.Factor.Copy();
        _spectrum = _ridge ? snapshot.Spectrum : null;
        if (_priorColumns is not null)
        {
            _priorGram = PriorGramOf(_priorColumns);
            _pencil = snapshot.Pencil;
            _spectrum = null;
        }
        _beta = (double[])snapshot.Beta.Clone();
        snapshot.Residual.CopyTo(_residual, 0);
        ResidualSquares = snapshot.ResidualSquares;
        FittedSquares = snapshot.FittedSquares;
    }

    /// <summary>Refactors the Gram strictly and refits; returns the drift the rebuild removed.</summary>
    public FactorDrift Rebuild()
    {
        _pending = null;
        FactorDrift drift = _factor.Rebuild();
        Refit();
        return drift;
    }

    // `border` is the new column's Gram entries against the candidate's other columns, `full` against every current one.
    private DesignProposal Decide(Move kind, int index, BassTerm term, double[] column, double[]? output, double[] border,
                                  double[] full, double diagonal, double bEntry, FactorProposal p, double rssReduced, double zrReduced,
                                  BassPriorColumn? prior = null)
    {
        double floor = _rules.AngleFloor;
        double[] b = kind == Move.Change ? Without(_b, index) : [.. _b];
        Array.Resize(ref b, b.Length + 1);
        b[^1] = bEntry;

        double minAngle, newest, q;
        double[,]? strict = null;
        if (p.PositiveDefinite && p.MinLeaveOneOutPivot >= 2.0 * floor)
        {
            minAngle = p.MinLeaveOneOutPivot;
            newest = p.RelativePivot;
            q = _factor.ProposedInverseQuadraticForm(b);
        }
        else
        {
            // The bordered factor fails the gate or comes near it: a strict factorization of the candidate decides.
            strict = Candidate(kind == Move.Change ? index : -1, border, diagonal);
            BorderedCholesky s = StrictFactor(strict);
            if (s.FailedAt >= 0) return DesignProposal.Refused(BassRefusal.Numerical);
            minAngle = s.MinLeaveOneOutPivot();
            if (minAngle < floor) return DesignProposal.Refused(BassRefusal.Numerical);
            newest = s.RelativePivot(s.Dimension - 1);
            q = s.InverseQuadraticForm(b);
        }
        if (_rules.MaxVarianceInflation is { } cap && minAngle * cap < 1.0) return DesignProposal.Refused(BassRefusal.VarianceInflation);

        double rss = rssReduced - zrReduced * zrReduced / (newest * diagonal);
        if (!(rss > ExplicitBelow * rssReduced))
            rss = ExplicitResidualSquares(kind == Move.Change ? index : -1, column, output, strict ?? Candidate(kind == Move.Change ? index : -1, border, diagonal), b);
        // The ridge's terms: bordered on the current spectrum, or, where the strict factorization decided, the candidate's
        // own decomposition, which an acceptance keeps.
        RidgeColumns? ridge = null;
        RidgeSpectrum? fresh = null;
        double[,]? priorGram = null;
        PriorPencil? pencil = null;
        if (_priorColumns is not null)
        {
            int skip = kind == Move.Change ? index : -1;
            priorGram = PriorGramWith(skip, prior!);
            double[] shifts = skip >= 0 ? Without(CurrentShifts(), skip) : CurrentShifts();
            Array.Resize(ref shifts, shifts.Length + 1);
            shifts[^1] = prior!.Shift;
            pencil = PencilOf(strict ?? Candidate(skip, border, diagonal), priorGram, shifts, b);
            if (pencil is null) return DesignProposal.Refused(BassRefusal.Numerical);
            ridge = pencil;
        }
        else if (_ridge)
            ridge = strict is not null ? fresh = RidgeSpectrum.Of(strict, b, _rows, Weighting.Mass)
                  : kind == Move.Birth ? Spectrum.Birth(full, diagonal, bEntry) : Spectrum.Change(index, full, diagonal, bEntry);
        var account = new DesignProposal(null, kind == Move.Birth ? _terms.Count + 1 : _terms.Count, Math.Max(rss, 0.0), q, minAngle, newest, ridge);
        _pending = new Pending(kind, index, term, column, output, bEntry, strict, account, fresh, prior, priorGram, pencil);
        return account;
    }

    // The term's columns, or null when it falls short of the support threshold in a part. Scalar: the centred case
    // column. Direct: the case-side and output-side factors as evaluated. Under weights: the whitened column over all N rows,
    // stacked for a direct response, and its prior side when the prior stands apart.
    private (double[] Column, double[]? Output, BassPriorColumn? Prior)? Evaluate(BassTerm term)
    {
        _basis.Validate(term);
        var a = new double[_cases];
        _basis.CaseColumn(term, _basis.Inputs.Cases, a);
        if (!Supported(a, _rules.SupportMinimum, _counted)) return null;
        if (!_functional)
        {
            if (_weights is not null) return Whitened(a);
            double mean = PairwiseSum.Of(a) / _cases;
            for (int i = 0; i < a.Length; i++) a[i] -= mean;
            return (a, null, null);
        }
        var h = new double[_locations];
        _basis.OutputColumn(term, _basis.Inputs.Outputs!, h);
        if (!Supported(h, _rules.SupportMinimumFunctional ?? 0.0, null)) return null;
        if (_weights is null) return (a, h, null);
        var stacked = new double[_rows];
        for (int i = 0; i < _cases; i++)
            for (int l = 0; l < _locations; l++) stacked[i * _locations + l] = a[i] * h[l];
        return Whitened(stacked);
    }

    // W^½(z − z̄_w), in place, and z's prior side when the prior stands apart.
    private (double[] Column, double[]? Output, BassPriorColumn? Prior) Whitened(double[] z)
    {
        double mean = PairwiseSum.Dot(_weights!, z) / Weighting.Mass;
        BassPriorColumn? prior = null;
        if (_priorColumns is not null)
        {
            double own = PairwiseSum.Of(z) / _rows;
            var centred = new double[_rows];
            for (int i = 0; i < _rows; i++) centred[i] = z[i] - own;
            prior = new BassPriorColumn(centred, own - mean);
        }
        for (int i = 0; i < _rows; i++) z[i] = _roots![i] * (z[i] - mean);
        return (z, null, prior);
    }

    private bool Supported(double[] part, double minimum, double[]? counts)
    {
        if (_rules.Support == SupportKind.ColumnMass) return (counts is null ? PairwiseSum.Of(part) : PairwiseSum.Dot(counts, part)) >= minimum;
        double nonzero = 0.0;
        for (int i = 0; i < part.Length; i++)
            if (part[i] != 0.0) nonzero += counts?[i] ?? 1.0;
        return nonzero >= minimum;
    }

    // The current columns' Gram entries, read from the factor.
    private double[,] CurrentGram()
    {
        int n = _terms.Count;
        var g = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) g[i, j] = _factor.Entry(i, j);
        return g;
    }

    private double[] CurrentShifts()
    {
        var shifts = new double[_priorColumns!.Count];
        for (int j = 0; j < shifts.Length; j++) shifts[j] = _priorColumns[j].Shift;
        return shifts;
    }

    // The prior side's Gram of `columns`.
    private static double[,] PriorGramOf(IReadOnlyList<BassPriorColumn> columns)
    {
        int n = columns.Count;
        var g = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++) g[i, j] = g[j, i] = PairwiseSum.Dot(columns[i].Centred, columns[j].Centred);
        return g;
    }

    // The prior side's Gram without column `index`.
    private double[,] PriorGramWithout(int index) => Without(_priorGram!, index);

    // The prior side's Gram with column `skip` out (−1: none) and `added` last.
    private double[,] PriorGramWith(int skip, BassPriorColumn added)
    {
        double[,] kept = skip >= 0 ? Without(_priorGram!, skip) : _priorGram!;
        int n = kept.GetLength(0);
        var g = new double[n + 1, n + 1];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) g[i, j] = kept[i, j];
        for (int j = 0, k = 0; j < _priorColumns!.Count; j++)
        {
            if (j == skip) continue;
            g[n, k] = g[k, n] = PairwiseSum.Dot(added.Centred, _priorColumns[j].Centred);
            k++;
        }
        g[n, n] = PairwiseSum.Dot(added.Centred, added.Centred);
        return g;
    }

    // The pencil of the structure whose weighted Gram is `g`, prior-side Gram `u`, centres' gaps `shifts` and cross products
    // `b`; null when its prior fails the gate. In the coordinates the design holds, (α, β) with the columns centred under W,
    // the likelihood is diag(Σw, G) and the prior D_α(1, δ)(1, δ)ᵀ + diag(0, D_β), since the prior reads α on the columns' own
    // centring, α_u = α + δᵀβ: D_α = N and D_β = G_u under the g-form, D_α = 1 and D_β = diag(G_u)/N under the ridge. A flat
    // intercept has no prior to couple, so it stands apart: G against D_β.
    private PriorPencil? PencilOf(double[,] g, double[,] u, double[] shifts, double[] b)
    {
        int m = g.GetLength(0);
        if (!_ridge && m > 0)
        {
            BorderedCholesky gate = StrictFactor(u);
            if (gate.FailedAt >= 0 || gate.MinLeaveOneOutPivot() < _rules.AngleFloor) return null;
        }
        bool inside = _intercept != GPriorCentre.Flat;
        int offset = inside ? 1 : 0, p = m + offset;
        var likelihood = new double[p, p];
        var prior = new double[p, p];
        var cross = new double[p];
        double priorMass = _ridge ? 1.0 : _rows;
        for (int i = 0; i < m; i++)
        {
            cross[i + offset] = b[i];
            for (int j = 0; j < m; j++)
            {
                likelihood[i + offset, j + offset] = g[i, j];
                prior[i + offset, j + offset] = !_ridge ? u[i, j] : i == j ? u[i, i] / _rows : 0.0;
                if (inside) prior[i + offset, j + offset] += priorMass * shifts[i] * shifts[j];
            }
            if (inside) prior[0, i + offset] = prior[i + offset, 0] = priorMass * shifts[i];
        }
        if (inside)
        {
            likelihood[0, 0] = Weighting.Mass;
            prior[0, 0] = priorMass;
            cross[0] = Weighting.Mass * (ResponseMean - PriorLevel);
        }
        return PriorPencil.Of(likelihood, prior, cross, m, inside);
    }

    // The new column's Gram entries against the current columns (all of them, or all but `skip`) and its diagonal.
    private double[] Border(double[] column, double[]? output, int skip, out double diagonal)
    {
        int count = _terms.Count - (skip >= 0 ? 1 : 0);
        var u = new double[count];
        if (!_khatriRao)
        {
            for (int j = 0, k = 0; j < _terms.Count; j++)
                if (j != skip) u[k++] = PairwiseSum.Dot(column, _a[j]);
            diagonal = PairwiseSum.Dot(column, column);
            return u;
        }
        diagonal = CentredSelf(column, output!);
        if (count == 0) return u;
        var a = new double[count][];
        var h = new double[count][];
        for (int j = 0, k = 0; j < _terms.Count; j++)
            if (j != skip) { a[k] = _a[j]; h[k++] = _h[j]; }
        return KhatriRao.CentredBorder(a, h, column, output!, out diagonal);
    }

    // ‖a ⊗ h − āh̄‖² by the three-term form: ‖α‖²‖η‖² + n·ā²‖η‖² + m·h̄²‖α‖², no cancellation.
    private double CentredSelf(double[] a, double[] h)
    {
        double aBar = PairwiseSum.Of(a) / a.Length, hBar = PairwiseSum.Of(h) / h.Length;
        var alpha = new double[a.Length];
        var eta = new double[h.Length];
        for (int i = 0; i < a.Length; i++) alpha[i] = a[i] - aBar;
        for (int r = 0; r < h.Length; r++) eta[r] = h[r] - hBar;
        double aa = PairwiseSum.Dot(alpha, alpha), hh = _output.Dot(eta, eta);
        return aa * hh + a.Length * aBar * aBar * hh + h.Length * hBar * hBar * aa;
    }

    // z_cᵀv for the new column and a vector v with zero mean over all N rows (ỹ, or the residual): scalar z·v; direct
    // aᵀVh, the centring dropping out against v's zero mean. Under weights both are whitened, so their dot is the W-form.
    private double Cross(double[] column, double[]? output, double[] v)
    {
        if (!_khatriRao) return PairwiseSum.Dot(column, v);
        var vh = new double[_cases];
        for (int i = 0; i < _cases; i++) vh[i] = _output.Dot(v.AsSpan(i * _locations, _locations), output);
        return PairwiseSum.Dot(column, vh);
    }

    private double CrossResponse(double[] column, double[]? output) => Cross(column, output, _y);

    private double CrossResidual(double[] column, double[]? output) => Cross(column, output, _residual);

    // The candidate's Gram: the current entries (without `skip`) and the new column last.
    private double[,] Candidate(int skip, double[] border, double diagonal)
    {
        int count = border.Length + 1;
        var g = new double[count, count];
        for (int i = 0, gi = 0; i < _terms.Count; i++)
        {
            if (i == skip) continue;
            for (int j = 0, gj = 0; j < _terms.Count; j++)
            {
                if (j == skip) continue;
                g[gi, gj++] = _factor.Entry(i, j);
            }
            gi++;
        }
        for (int k = 0; k < border.Length; k++) g[count - 1, k] = g[k, count - 1] = border[k];
        g[count - 1, count - 1] = diagonal;
        return g;
    }

    private double[,] Gram(List<double[]> columns, List<double[]> outputs)
    {
        if (columns.Count == 0) return new double[0, 0];
        if (_khatriRao) return KhatriRao.CentredGram([.. columns], [.. outputs]);
        var g = new double[columns.Count, columns.Count];
        for (int i = 0; i < columns.Count; i++)
            for (int j = 0; j <= i; j++) g[i, j] = g[j, i] = PairwiseSum.Dot(columns[i], columns[j]);
        return g;
    }

    private BorderedCholesky StrictFactor(double[,] g)
    {
        if (_strict is null || _strict.Capacity < g.GetLength(0)) _strict = new BorderedCholesky(Math.Max(_factor.Capacity, g.GetLength(0)), 0);
        _strict.Decompose(g);
        return _strict;
    }

    // The candidate's residual sum of squares from its explicit residual, for a move that explains nearly all of it.
    private double ExplicitResidualSquares(int skip, double[] column, double[]? output, double[,] g, double[] b)
    {
        BorderedCholesky s = StrictFactor(g);
        if (s.FailedAt >= 0) return 0.0;
        double[] beta = s.Solve(b);
        var a = new List<double[]>();
        var h = new List<double[]>();
        for (int j = 0; j < _terms.Count; j++)
            if (j != skip) { a.Add(_a[j]); if (_khatriRao) h.Add(_h[j]); }
        a.Add(column);
        if (_khatriRao) h.Add(output!);
        var r = new double[_rows];
        Fit(a, h, beta, r);
        return _output.Dot(r, r);
    }

    // r = ỹ − Z_cβ, into r.
    private void Fit(List<double[]> a, List<double[]> h, double[] beta, double[] r)
    {
        Array.Copy(_y, r, _rows);
        for (int j = 0; j < beta.Length; j++)
        {
            double bj = beta[j];
            if (!_khatriRao)
            {
                double[] z = a[j];
                for (int i = 0; i < _rows; i++) r[i] -= bj * z[i];
                continue;
            }
            double[] aj = a[j], hj = h[j];
            double level = PairwiseSum.Of(aj) / _cases * (PairwiseSum.Of(hj) / _locations);
            for (int i = 0; i < _cases; i++)
            {
                double ai = aj[i];
                int row = i * _locations;
                for (int l = 0; l < _locations; l++) r[row + l] -= bj * (ai * hj[l] - level);
            }
        }
    }

    private void Refit()
    {
        double[] b = [.. _b];
        _beta = _terms.Count == 0 ? [] : _factor.Solve(b);
        Fit(_a, _h, _beta, _residual);
        ResidualSquares = _output.Dot(_residual, _residual);
        FittedSquares = _terms.Count == 0 ? 0.0 : PairwiseSum.Dot(b, _beta);
    }

    private void EnsureCapacity(int order)
    {
        if (order <= _factor.Capacity) return;
        int capacity = Math.Min(Math.Max(_rules.MaxBasis, order), Math.Max(2 * _factor.Capacity, order));
        var grown = new BorderedCholesky(capacity, _factor.RebuildInterval);
        int n = _terms.Count;
        var g = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) g[i, j] = _factor.Entry(i, j);
        grown.Decompose(g);
        _factor = grown;
    }

    private static double[] Without(List<double> values, int index)
    {
        var w = new double[values.Count - 1];
        for (int j = 0, k = 0; j < values.Count; j++)
            if (j != index) w[k++] = values[j];
        return w;
    }

    private static double[] Without(double[] values, int index)
    {
        var w = new double[values.Length - 1];
        for (int j = 0, k = 0; j < values.Length; j++)
            if (j != index) w[k++] = values[j];
        return w;
    }

    private static double[,] Without(double[,] values, int index)
    {
        int n = values.GetLength(0) - 1;
        var w = new double[n, n];
        for (int i = 0, gi = 0; i <= n; i++)
        {
            if (i == index) continue;
            for (int j = 0, gj = 0; j <= n; j++)
            {
                if (j == index) continue;
                w[gi, gj++] = values[i, j];
            }
            gi++;
        }
        return w;
    }
}
