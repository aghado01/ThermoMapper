using System;
using System.Collections.Generic;
using Maths.LinAlg;

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
/// (null in a checkpoint, which restores them from the terms).</summary>
internal sealed record BassDesignSnapshot(BassTerm[] Terms, BorderedCholesky Factor, double[] Beta, double[] Residual,
                                          double ResidualSquares, double FittedSquares, BassDesignColumns? Columns = null);

/// <summary>A design's evaluated columns in column order: each term's case side, its output side for a direct response,
/// and its cross product with the response. A design never writes into a column once evaluated, so designs over the
/// same data can share them.</summary>
internal sealed record BassDesignColumns(double[][] Cases, double[][]? Outputs, double[] CrossResponses);

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
/// residual and fitted sums of squares in the centred basis, its smallest leave-one-out angle (D06's gate) and the new
/// column's angle to the span of the others (NaN for a death).
/// </summary>
public readonly record struct DesignProposal(
    BassRefusal? Refusal, int Columns, double ResidualSquares, double FittedSquares, double MinAngle, double NewestAngle)
{
    /// <summary>Whether the structure lies in F.</summary>
    public bool Admitted => Refusal is null;

    internal static DesignProposal Refused(BassRefusal reason)
        => new(reason, -1, double.NaN, double.NaN, double.NaN, double.NaN);
}

/// <summary>
/// One chain's design: its terms' columns in the centred basis, the factor of their Gram, and the least-squares fit of
/// the centred response, maintained through births, deaths and changes (stage 4's S4; D06; the numerics segment's
/// dense factor maintenance).
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
/// Not thread-safe: one per chain.
/// </summary>
public sealed class BassDesign
{
    private const int InitialCapacity = 32;
    private const double ExplicitBelow = 1e-8;   // ‖r′‖²/‖r_reduced‖² under which the candidate's residual is formed

    private enum Move { Birth, Death, Change }

    private sealed record Pending(Move Kind, int Index, BassTerm? Term, double[]? Column, double[]? Output, double BEntry,
                                  double[,]? Strict, DesignProposal Account);

    private readonly BassBasis _basis;
    private readonly BassAdmission _rules;
    private readonly bool _functional;
    private readonly int _cases, _locations, _rows;
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

    /// <summary>The intercept-only design over <paramref name="basis"/>'s inputs for <paramref name="data"/>'s response.</summary>
    public BassDesign(BassBasis basis, BassData data, BassAdmission rules)
    {
        ArgumentNullException.ThrowIfNull(basis);
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(rules);
        _basis = basis;
        _rules = rules;
        _functional = basis.Inputs.Outputs is not null;
        if (_functional != data.Response is FunctionalResponse)
            throw new ArgumentException("The basis and the data disagree on the response's kind.", nameof(data));
        _cases = data.Response.Cases;
        _locations = data.Response.Locations;
        _rows = data.Response.Observations;
        if (basis.Inputs.Cases.Cases != _cases) throw new ArgumentException("The basis was not prepared from these data.", nameof(basis));
        ReadOnlySpan<double> y = data.Response.Values;
        ResponseMean = PairwiseSum.Of(y) / _rows;
        _y = new double[_rows];
        for (int i = 0; i < _rows; i++) _y[i] = y[i] - ResponseMean;
        ResponseSquares = PairwiseSum.Dot(_y, _y);
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

    /// <summary>ȳ, the response's mean over all N observations.</summary>
    public double ResponseMean { get; }

    /// <summary>‖ỹ‖² = ‖y − ȳ‖².</summary>
    public double ResponseSquares { get; }

    /// <summary>‖r‖², the fit's residual sum of squares, from the explicit residual.</summary>
    public double ResidualSquares { get; private set; }

    /// <summary>Q_c = b_cᵀG_c⁻¹b_c, the fitted sum of squares in the centred basis.</summary>
    public double FittedSquares { get; private set; }

    /// <summary>β̂_c, the least-squares coefficients of the centred columns.</summary>
    public ReadOnlySpan<double> Coefficients => _beta;

    /// <summary>The explicit residual r = ỹ − Z_cβ̂, case by case and location by location.</summary>
    public ReadOnlySpan<double> Residual => _residual;

    /// <summary>The factor of G_c, in the terms' order.</summary>
    public ISpdFactor Factor => _factor;

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
        if (Evaluate(term) is not var (column, output)) return DesignProposal.Refused(BassRefusal.Support);
        EnsureCapacity(_terms.Count + 1);
        double[] border = Border(column, output, -1, out double diagonal);
        double bEntry = CrossResponse(column, output);
        FactorProposal p = _factor.ProposeBirth(border, diagonal);
        return Decide(Move.Birth, _terms.Count, term, column, output, border, diagonal, bEntry, p,
                      ResidualSquares, CrossResidual(column, output));
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
        var account = new DesignProposal(null, _terms.Count - 1, rss, q, p.MinLeaveOneOutPivot, double.NaN);
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
            var same = new DesignProposal(null, _terms.Count, ResidualSquares, FittedSquares, MinAngle, _factor.LeaveOneOutPivot(index));
            _pending = new Pending(Move.Change, index, term, null, null, 0.0, null, same);
            return same;
        }
        if (_terms.Contains(term)) return DesignProposal.Refused(BassRefusal.Numerical);
        if (Evaluate(term) is not var (column, output)) return DesignProposal.Refused(BassRefusal.Support);
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
        return Decide(Move.Change, index, term, column, output, border, diagonal, bEntry, p, rssReduced, zrReduced);
    }

    /// <summary>Commits the last proposal, which must have been admitted.</summary>
    public void Accept()
    {
        Pending p = _pending ?? throw new InvalidOperationException("Nothing admitted is proposed.");
        _pending = null;
        if (p.Kind == Move.Change && p.Column is null) return;   // a change to the same term
        if (p.Strict is not null) _factor.Decompose(p.Strict);
        else _factor.Accept();
        if (p.Kind is Move.Death or Move.Change)
        {
            _terms.RemoveAt(p.Index);
            _a.RemoveAt(p.Index);
            if (_functional) _h.RemoveAt(p.Index);
            _b.RemoveAt(p.Index);
        }
        if (p.Kind is Move.Birth or Move.Change)
        {
            _terms.Add(p.Term!);
            _a.Add(p.Column!);
            if (_functional) _h.Add(p.Output!);
            _b.Add(p.BEntry);
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
        foreach (BassTerm term in terms)
        {
            if (!seen.Add(term)) return BassRefusal.Numerical;
            if (Evaluate(term) is not var (column, output)) return BassRefusal.Support;
            columns.Add(column);
            if (_functional) outputs.Add(output!);
        }
        double[,] g = Gram(columns, outputs);
        var fresh = new BorderedCholesky(Math.Max(_factor.Capacity, terms.Count));
        fresh.Decompose(g);
        if (terms.Count > 0 && fresh.FailedAt >= 0) return BassRefusal.Numerical;
        double minAngle = fresh.MinLeaveOneOutPivot();
        if (minAngle < _rules.AngleFloor) return BassRefusal.Numerical;
        if (_rules.MaxVarianceInflation is { } cap && minAngle * cap < 1.0) return BassRefusal.VarianceInflation;

        _factor = fresh;
        _terms.Clear(); _terms.AddRange(terms);
        _a.Clear(); _a.AddRange(columns);
        _h.Clear(); _h.AddRange(outputs);
        _b.Clear();
        for (int j = 0; j < columns.Count; j++) _b.Add(CrossResponse(columns[j], _functional ? outputs[j] : null));
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
               shareColumns ? new BassDesignColumns([.. _a], _functional ? [.. _h] : null, [.. _b]) : null);

    /// <summary>Replaces the structure with <paramref name="snapshot"/>'s, every statistic as the design it was taken from
    /// held it, so the two continue bit for bit; its columns shared when it carries them, else re-evaluated.</summary>
    internal void Restore(BassDesignSnapshot snapshot)
    {
        ArgumentNullException.ThrowIfNull(snapshot);
        if (snapshot.Factor.Dimension != snapshot.Terms.Length || snapshot.Beta.Length != snapshot.Terms.Length || snapshot.Residual.Length != _rows)
            throw new ArgumentException("The snapshot is not a design's over these data.", nameof(snapshot));
        if (snapshot.Columns is { } shared
            && (shared.Cases.Length != snapshot.Terms.Length || shared.CrossResponses.Length != snapshot.Terms.Length
                || (shared.Outputs is null) == _functional || shared.Outputs is { } o && o.Length != snapshot.Terms.Length
                || Array.Exists(shared.Cases, c => c.Length != _cases)
                || shared.Outputs is { } h && Array.Exists(h, c => c.Length != _locations)))
            throw new ArgumentException("The snapshot's columns are not its terms' over these data.", nameof(snapshot));
        _pending = null;
        _terms.Clear();
        _a.Clear();
        _h.Clear();
        _b.Clear();
        if (snapshot.Columns is { } columns)
        {
            _terms.AddRange(snapshot.Terms);
            _a.AddRange(columns.Cases);
            if (_functional) _h.AddRange(columns.Outputs!);
            _b.AddRange(columns.CrossResponses);
        }
        else foreach (BassTerm term in snapshot.Terms)
        {
            if (Evaluate(term) is not var (column, output)) throw new ArgumentException($"The snapshot holds {term}, which lacks support here.", nameof(snapshot));
            _terms.Add(term);
            _a.Add(column);
            if (_functional) _h.Add(output!);
            _b.Add(CrossResponse(column, output));
        }
        _factor = snapshot.Factor.Copy();
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

    private DesignProposal Decide(Move kind, int index, BassTerm term, double[] column, double[]? output, double[] border,
                                  double diagonal, double bEntry, FactorProposal p, double rssReduced, double zrReduced)
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
        var account = new DesignProposal(null, kind == Move.Birth ? _terms.Count + 1 : _terms.Count, Math.Max(rss, 0.0), q, minAngle, newest);
        _pending = new Pending(kind, index, term, column, output, bEntry, strict, account);
        return account;
    }

    // The term's columns, or null when it falls short of the support threshold in a part. Scalar: the centred case
    // column. Direct: the case-side and output-side factors as evaluated.
    private (double[] Column, double[]? Output)? Evaluate(BassTerm term)
    {
        _basis.Validate(term);
        var a = new double[_cases];
        _basis.CaseColumn(term, _basis.Inputs.Cases, a);
        if (!Supported(a, _rules.SupportMinimum)) return null;
        if (!_functional)
        {
            double mean = PairwiseSum.Of(a) / _cases;
            for (int i = 0; i < a.Length; i++) a[i] -= mean;
            return (a, null);
        }
        var h = new double[_locations];
        _basis.OutputColumn(term, _basis.Inputs.Outputs!, h);
        if (!Supported(h, _rules.SupportMinimumFunctional ?? 0.0)) return null;
        return (a, h);
    }

    private bool Supported(double[] part, double minimum)
    {
        if (_rules.Support == SupportKind.ColumnMass) return PairwiseSum.Of(part) >= minimum;
        int nonzero = 0;
        foreach (double v in part)
            if (v != 0.0) nonzero++;
        return nonzero >= minimum;
    }

    // The new column's Gram entries against the current columns (all of them, or all but `skip`) and its diagonal.
    private double[] Border(double[] column, double[]? output, int skip, out double diagonal)
    {
        int count = _terms.Count - (skip >= 0 ? 1 : 0);
        var u = new double[count];
        if (!_functional)
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
        double aa = PairwiseSum.Dot(alpha, alpha), hh = PairwiseSum.Dot(eta, eta);
        return aa * hh + a.Length * aBar * aBar * hh + h.Length * hBar * hBar * aa;
    }

    // z_cᵀv for the new column and a vector v with zero mean over all N rows (ỹ, or the residual): scalar z·v; direct
    // aᵀVh, the centring dropping out against v's zero mean.
    private double Cross(double[] column, double[]? output, double[] v)
    {
        if (!_functional) return PairwiseSum.Dot(column, v);
        var vh = new double[_cases];
        for (int i = 0; i < _cases; i++) vh[i] = PairwiseSum.Dot(v.AsSpan(i * _locations, _locations), output);
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
        if (_functional) return KhatriRao.CentredGram([.. columns], [.. outputs]);
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
            if (j != skip) { a.Add(_a[j]); if (_functional) h.Add(_h[j]); }
        a.Add(column);
        if (_functional) h.Add(output!);
        var r = new double[_rows];
        Fit(a, h, beta, r);
        return PairwiseSum.Dot(r, r);
    }

    // r = ỹ − Z_cβ, into r.
    private void Fit(List<double[]> a, List<double[]> h, double[] beta, double[] r)
    {
        Array.Copy(_y, r, _rows);
        for (int j = 0; j < beta.Length; j++)
        {
            double bj = beta[j];
            if (!_functional)
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
        ResidualSquares = PairwiseSum.Dot(_residual, _residual);
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
}
