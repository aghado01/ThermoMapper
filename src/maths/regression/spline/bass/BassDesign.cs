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
/// prior apart from its Gram its pencil, immutable values the restored design shares; and what it reads of the working response it
/// was on, which a chain may have moved from the model's.</summary>
internal sealed record BassDesignSnapshot(BassTerm[] Terms, BorderedCholesky Factor, double[] Beta, double[] Residual,
                                          double ResidualSquares, double FittedSquares, BassDesignColumns? Columns = null,
                                          RidgeSpectrum? Spectrum = null, PriorPencil? Pencil = null, BassDesignResponse? Response = null,
                                          DesignOutputsSnapshot? Outputs = null);

/// <summary>What a design reads of one working response: the response itself, W^½ (null at unit weights), ỹ (whitened under weights),
/// the level a prior apart from the Gram centres on, and ‖ỹ‖²_W. Immutable, so designs that hand a state over share it.</summary>
internal sealed record BassDesignResponse(WorkingResponse Working, double[] Y, double[]? Roots, double PriorLevel, double Squares);

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
/// its columns' terms at any (τ, b) (<see cref="RidgeColumns"/>; null under the g-form on the design's own Gram). Several outputs (E34)
/// carry their K × K forms in <see cref="Outputs"/>, the sums of squares then NaN.
/// </summary>
public readonly record struct DesignProposal(
    BassRefusal? Refusal, int Columns, double ResidualSquares, double FittedSquares, double MinAngle, double NewestAngle,
    RidgeColumns? Ridge = null, OutputForms? Outputs = null)
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
/// <para><b>The prior's precision</b> comes from the design's builder (<see cref="BassPriorPrecision"/>, E21), one per coefficient
/// form. KHAOS's two stand apart from the Gram at any weights. Its raw ridge (<see cref="RawRidgeForm"/>) reads the coefficients on
/// the columns as evaluated, uncentred: D = I on (β₀, β), with β₀ = α − z̄_wᵀβ in the design's coordinates, so its pencil has
/// D_α = 1, δ = −z̄_w and D_β = I; it reads no prior-side columns or Gram, and adds nothing to admission, so F is the g-form's. Its
/// modified g-prior (<see cref="ModifiedGForm"/>) is the g-form's pencil on the prior-side columns scaled by 1/g_m, whose Gram
/// clears the angle gate as the g-form's does; scaling a column leaves its angles alone, so F is the g-form's there too.</para>
/// <para><b>When the working response moves</b> (<see cref="ProposeFollow"/>) the design refits its current structure on the new one,
/// its own record of the response replaced whole. F does not read the values, so an offset alone keeps the columns, the factor and
/// the gate. Weights whiten the columns again and factor their Gram strictly, and the gate reads them, so a move can put the
/// current structure outside F; that proposal is refused and the design stays as it was.</para>
/// Not thread-safe: one per chain.
/// </summary>
public sealed class BassDesign
{
    private const int InitialCapacity = 32;
    private const double ExplicitBelow = 1e-8;   // ‖r′‖²/‖r_reduced‖² under which the candidate's residual is formed

    private enum Move { Birth, Death, Change, Follow }

    private sealed record Pending(Move Kind, int Index, BassTerm? Term, double[]? Column, double[]? Output, double BEntry,
                                  double[,]? Strict, DesignProposal Account, RidgeSpectrum? Fresh = null, BassPriorColumn? Prior = null,
                                  double[,]? PriorGram = null, PriorPencil? Pencil = null, Followed? Response = null, double[]? OutputCross = null);

    // A move to another working response, pending: what it reads, and the current structure's fit on it. The columns and the factor
    // are null when only the offset moved, which keeps them.
    private sealed record Followed(BassDesignResponse Response, List<double[]>? Columns, List<BassPriorColumn>? Priors, BorderedCholesky? Factor, double[] B,
                                   double[] Beta, double[] Residual, double ResidualSquares, double FittedSquares, RidgeSpectrum? Spectrum, PriorPencil? Pencil,
                                   DesignOutputs? Outputs = null);

    private readonly BassBasis _basis;
    private readonly BassAdmission _rules;
    private readonly bool _functional, _ridge;
    private readonly BassPriorPrecision _precision; // how the prior's precision reads the columns (E21)
    // Whether the reading the design is on carries a random-effect block's whitening, by which each column is whitened (E21).
    private bool IsWhitened() => _response.Working.Whitening is not null;
    private readonly bool _khatriRao;               // a direct response held as case and output factors, its rows unweighted
    private readonly int _cases, _locations, _rows;
    private BassDesignResponse _response;           // what the design reads of the working response it is on
    private readonly double[]? _counted;            // a scalar response's weights, when they count rows for the support threshold
    private readonly GPriorCentre _intercept;       // what a prior apart from the Gram does with the intercept
    private readonly List<BassPriorColumn>? _priorColumns;   // the prior's side, when it stands apart from the Gram
    private readonly double[]? _priorWeights, _priorRoots;   // W₀ and W₀^½ when the prior's Gram is frozen on them; null otherwise
    private readonly double _priorMass;             // Σw₀, N with no W₀: the prior side's centring and the ridge's normalization
    private double[,]? _priorGram;                  // G_c on the prior's side
    private PriorPencil? _pencil;                   // the current structure's; null until read
    private readonly IInnerProduct _output;         // Q on the output index; the identity for independent noise
    private double[] _residual;
    private readonly List<BassTerm> _terms = [];
    private readonly List<double[]> _a = [];        // scalar: centred columns; direct: case-side factors
    private readonly List<double[]> _h = [];        // direct: output-side factors
    private readonly List<double> _b = [];          // Z_cᵀỹ
    private double[] _beta = [];
    private BorderedCholesky _factor;
    private BorderedCholesky? _strict;
    private Pending? _pending;
    private RidgeSpectrum? _spectrum;   // the current structure's, under the ridge; null until read
    private DesignOutputs? _outputs;   // several outputs' statistics, in place of the scalar ones (E34); null for one response; re-read when the reading moves
    private OutputPencil? _outputPencil;

    // Lever 1 of the efficiency pass: the cross products of a one-way hinge with the held columns, from prefix tables. A table is a pure
    // cache, per (column, variable, sign), built the first time a proposal needs it and dropped with its column; which proposals read tables
    // is decided once, by the design's size, never by what the cache holds, so a restored design, a handed-over one and a live one compute
    // the same numbers. Tables serve the leading variables whose tables, at the cap, fit TableBytesAllowed.
    private const int MaxTableVariables = 16;   // past this a (column, variable, sign) is reached too seldom for a table to pay (measured)
    private const int MinTableCases = 500;      // below this a dot is cheaper than the scan's setup (measured)
    private const long DefaultTableBytes = 64L << 20;
    private const double NearGate = 1e4;   // a table proposal whose bordered minimum angle is under this many floors takes the dots
    private readonly List<ColumnTables?> _tables = [];   // aligned with _a
    private HingeOrder?[]? _orders;                       // per continuous variable, built when first read
    private bool[]? _ordersTried;
    private int _tableVariables = -1;

    private sealed class ColumnTables(int variables)
    {
        public readonly double[]?[] Above = new double[]?[variables], Below = new double[]?[variables];
    }

    // BASS_TABLES = off or on, for a measurement to set for a whole process; unset, a plain chain asks for tables (EnableTables) and the
    // design's size decides.
    private static readonly string? TablesSetting = Environment.GetEnvironmentVariable("BASS_TABLES");

    /// <summary>Whether a one-way hinge's cross products come from prefix tables. Off for a design that is not asked: a table costs about
    /// a dot's worth of work ten times over to build, so it pays only to a design that lives long enough to read it back, which a chain's does and a
    /// particle's step, a handed-over replica's or a draw's does not.</summary>
    internal bool TablesOn { get; set; } = TablesSetting == "on";

    /// <summary>Asks for tables on a chain's design, which they serve when the fit is large enough and not too wide
    /// (<see cref="MaxTableVariables"/> continuous predictors at most, <see cref="MinTableCases"/> cases at least); <c>BASS_TABLES=off</c> refuses.</summary>
    internal void EnableTables()
    {
        if (TablesSetting == "off" || _outputs is not null) return;   // several outputs read exact borders (E34's slice 1)
        TablesOn = TablesSetting == "on" || (_basis.Inputs.Cases.ContinuousCount <= MaxTableVariables && _cases >= MinTableCases);
    }

    /// <summary>Builds every table the design could be asked for: each column's, for every variable the budget serves and both signs. A
    /// chain reaches this state in its first few hundred sweeps; a measurement of a restored chain begins from it.</summary>
    internal void WarmTables()
    {
        int degree = _basis.Degree;
        if (!TablesOn || _functional || _khatriRao || _priorColumns is not null || _counted is not null || IsWhitened() || degree is not (1 or 2)) return;
        int variables = _basis.Inputs.Cases.ContinuousCount, served = Math.Min(variables, TableVariables());
        double[]? roots = _response.Roots;
        for (int v = 0; v < served; v++)
        {
            if (OrderOf(v) is not { } order) continue;
            for (int j = 0; j < _a.Count; j++)
            {
                ColumnTables tables = _tables[j] ??= new ColumnTables(variables);
                if (tables.Above[v] is null) { tables.Above[v] = HingeScan.Above(_a[j], roots, order, degree); TableBuilds++; }
                if (tables.Below[v] is null) { tables.Below[v] = HingeScan.Below(_a[j], roots, order, degree); TableBuilds++; }
            }
        }
    }

    /// <summary>The bytes the tables may hold at the cap; the leading variables that fit are served.</summary>
    internal long TableBytesAllowed
    {
        get => _tableBytesAllowed;
        set
        {
            _tableBytesAllowed = value;
            _tableVariables = -1;
        }
    }

    private long _tableBytesAllowed = DefaultTableBytes;

    /// <summary>The border of <paramref name="term"/>, if a one-way hinge that tables serve, from the tables and from the dots; null
    /// where tables do not apply.</summary>
    internal (double[] Tabled, double[] Dotted)? BorderComparison(BassTerm term)
    {
        if (Evaluate(term, _response) is not var (column, output, _)) return null;
        double[]? tabled = TableBorder(term, column, out _, out bool approximate);
        return tabled is null || !approximate ? null : (tabled, Border(column, output, -1, out _));
    }

    /// <summary>Table entries read, dot entries beside them (a variable no table serves), tables built, and proposals the dots re-decided
    /// near the gate.</summary>
    internal long TableLookups, TableBuilds, TableFallbacks;

    /// <summary>The intercept-only design over <paramref name="basis"/>'s inputs for <paramref name="working"/>, the
    /// observation model's working response for <paramref name="data"/>'s (the response itself when null);
    /// <paramref name="ridge"/> keeps the normalized ridge's spectrum and proposal terms beside the factor; and
    /// <paramref name="intercept"/> places the intercept for a prior apart from the Gram, which alone reads it; and
    /// <paramref name="precision"/> is the prior-precision builder: the g-form's or the normalized ridge's, as
    /// <paramref name="ridge"/> says, when null, and one that agrees with <paramref name="ridge"/> when given. A random-effect block's
    /// whitening comes with the working response (<see cref="WorkingResponse.Whitening"/>).</summary>
    public BassDesign(BassBasis basis, BassData data, BassAdmission rules, bool ridge = false, WorkingResponse? working = null,
                      GPriorCentre intercept = GPriorCentre.Mean, BassPriorPrecision? precision = null)
    {
        ArgumentNullException.ThrowIfNull(basis);
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(rules);
        _basis = basis;
        _rules = rules;
        _ridge = ridge;
        _precision = precision ?? (ridge ? BassPriorPrecision.NormalizedRidge : BassPriorPrecision.Gram);
        if (_precision.IsRidge != ridge) throw new ArgumentException("The ridge flag and the prior-precision builder disagree.", nameof(precision));
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
        _khatriRao = _functional && working.WeightArray is null;
        if (working.WeightArray is { } weights)
        {
            if (!_functional && working.Meaning == WeightMeaning.Replication) _counted = weights;
            if (working.PriorStandsApart)
            {
                if (intercept is not (GPriorCentre.Mean or GPriorCentre.Zero or GPriorCentre.Flat)) throw new ArgumentOutOfRangeException(nameof(intercept));
                _priorColumns = [];
                _priorGram = new double[0, 0];
                _priorMass = _rows;
                if (working.Reference?.GramWeightArray is { } frozen)
                {
                    _priorWeights = frozen;
                    _priorMass = PairwiseSum.Of(frozen);
                    _priorRoots = new double[_rows];
                    for (int i = 0; i < _rows; i++) _priorRoots[i] = Math.Sqrt(frozen[i]);
                }
            }
        }
        if (_precision.AlwaysApart)
        {
            if (_functional) throw new ArgumentException("KHAOS's priors read a scalar response's columns as evaluated.", nameof(precision));
            if (intercept is not (GPriorCentre.Mean or GPriorCentre.Zero or GPriorCentre.Flat)) throw new ArgumentOutOfRangeException(nameof(intercept));
            if (!_precision.Unit && working.WeightArray is not null && !working.PriorStandsApart)
                throw new NotSupportedException("KHAOS's modified g-prior reads the unweighted or a frozen prior Gram; the weighted one is not built for it.");
            if (_priorColumns is null)
            {
                _priorColumns = [];
                _priorMass = _rows;
            }
            _priorGram = new double[0, 0];
        }
        if (data.Response is MultiOutputResponse)
        {
            // Several outputs (E34) read each output as the scalar design reads its response; under a random-effect block (slice 4) its
            // whitening and the g-form, under weights the g-form on the weighted Gram.
            if (_functional) throw new NotSupportedException("Several outputs need a case design.");
            if (working.Weighted is not null)
                throw new NotSupportedException("Several outputs under a noise's weights with random effects arrive with E34's slice 4, its last stage.");
            if ((working.Whitening is not null || working.WeightArray is not null) && (_ridge || _priorColumns is not null || !_precision.Unit && !ReferenceEquals(_precision, BassPriorPrecision.Gram)))
                throw new NotSupportedException("Several outputs take the g-form on the reading's own Gram under a random-effect block or weights; the "
                                                + "normalized ridge, KHAOS's priors and a prior apart from the weights arrive with the noise by row (E34's slice 4).");
            if (_ridge)
            {
                _priorColumns = [];
                _priorMass = _rows;
                _priorGram = new double[0, 0];
            }
            double[]? roots = working.WeightArray is { } w ? [.. w.Select(Math.Sqrt)] : null;
            _outputs = new DesignOutputs(ReadOutputs(working, roots));
            _response = new BassDesignResponse(working, new double[_rows], roots, 0.0, 0.0);
        }
        else _response = Read(working);
        _residual = (double[])_response.Y.Clone();
        ResidualSquares = _response.Squares;
        _factor = new BorderedCholesky(Math.Min(InitialCapacity, Math.Max(rules.MaxBasis, 1)));
    }

    // What the design reads of `working`, W^½ taken from `sameWeights` when that response's weights are these.
    private BassDesignResponse Read(WorkingResponse working, BassDesignResponse? sameWeights = null)
    {
        if (_outputs is not null)
        {
            // Several outputs (E34): the K-output side reads the outputs (ReadOutputs); this carries the reading and W^½.
            double[]? several = sameWeights?.Roots ?? (working.WeightArray is { } w ? [.. w.Select(Math.Sqrt)] : null);
            return new BassDesignResponse(working, new double[_rows], several, 0, 0);
        }
        ReadOnlySpan<double> values = working.Values;
        double mean = working.ResponseMean;
        if (working.Weighted is { } weighted)
        {
            // A random-effect block under weights (E21's slice 4): L(y − offset), centred against c = L·1 at the GLS mean.
            double[] t = weighted.Apply(values);
            ReadOnlySpan<double> c = weighted.Direction;
            for (int i = 0; i < _rows; i++) t[i] -= c[i] * mean;
            double[] rootsOf = sameWeights?.Roots ?? [.. working.WeightArray!.Select(Math.Sqrt)];
            double centre = _priorColumns is null || _intercept == GPriorCentre.Zero ? 0.0 : working.Reference?.Level ?? PairwiseSum.Of(values) / _rows;
            return new BassDesignResponse(working, t, rootsOf, centre, _output.Dot(t, t));
        }
        var y = new double[_rows];
        for (int i = 0; i < _rows; i++) y[i] = values[i] - mean;
        double[]? roots = null;
        if (working.WeightArray is { } weights)
        {
            roots = sameWeights?.Roots;
            if (roots is null)
            {
                roots = new double[_rows];
                for (int i = 0; i < _rows; i++) roots[i] = Math.Sqrt(weights[i]);
            }
            for (int i = 0; i < _rows; i++) y[i] *= roots[i];
        }
        double level = _priorColumns is null || _intercept == GPriorCentre.Zero ? 0.0 : working.Reference?.Level ?? PairwiseSum.Of(values) / _rows;
        return new BassDesignResponse(working, y, roots, level, _output.Dot(y, y));
    }

    // Several outputs (E34): each output read as Read reads a scalar response, ỹ_k = W^½(y_k − ȳ_k) on the reading's values (the block's
    // whitening in its offset) and means, at unit weights y_k − ȳ_k.
    private double[][] ReadOutputs(WorkingResponse working, double[]? roots)
    {
        int k = ((MultiOutputResponse)working.Response).Outputs;
        ReadOnlySpan<double> values = working.Values, means = working.OutputMeans;
        var reads = new double[k][];
        for (int l = 0; l < k; l++)
        {
            var y = new double[_rows];
            double mean = means[l];
            for (int i = 0; i < _rows; i++) y[i] = values[i * k + l] - mean;
            if (roots is not null)
                for (int i = 0; i < _rows; i++) y[i] *= roots[i];
            reads[l] = y;
        }
        return reads;
    }

    /// <summary>M, the terms.</summary>
    public int Count => _terms.Count;

    /// <summary>The terms, in the factor's column order.</summary>
    public IReadOnlyList<BassTerm> Terms => _terms;

    /// <summary>N, the observations.</summary>
    public int Observations => _rows;

    /// <summary>The working response the design is on: its model's, or the one a chain has moved it to.</summary>
    public WorkingResponse Working => _response.Working;

    /// <summary>The count, the mass and Σ ln wᵢ the marginal reads beyond these statistics; (N, N, 0) at unit weights.</summary>
    public GPriorWeights Weighting => _response.Working.Weighting;

    /// <summary>ȳ, the response's mean over all N observations under W.</summary>
    public double ResponseMean => _response.Working.ResponseMean;

    /// <summary>The level a prior apart from the Gram centres the intercept on, free of the weights; 0 otherwise.</summary>
    public double PriorLevel => _response.PriorLevel;

    /// <summary>‖ỹ‖²_W = ‖y − ȳ‖²_W.</summary>
    public double ResponseSquares => _response.Squares;

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

    /// <summary>The current structure's K × K forms for several outputs (E34); null for one response, whose sums of squares are
    /// <see cref="ResidualSquares"/> and <see cref="FittedSquares"/>.</summary>
    public OutputForms? OutputForms
    {
        get
        {
            if (_outputs is null) return null;
            var forms = _outputs.Forms;
            if (_priorColumns is null) return forms;
            _outputPencil ??= OutputPencil.Of(Pencil!, _outputs.CandidateCross(-1, null), OutputMeans, _intercept, _rows);
            return forms with { Prior = _outputPencil };
        }
    }

    /// <summary>β̂ for several outputs, one row per term in the factor's order, one column per output; null for one response.</summary>
    public double[,]? OutputCoefficients => _outputs?.Coefficients;

    /// <summary>The K outputs' means, which the centred basis reads them about; empty for one response.</summary>
    public ReadOnlySpan<double> OutputMeans => _response.Working.OutputMeans;

    /// <summary>Whether the coefficient prior's form is the normalized ridge.</summary>
    public bool IsRidge => _ridge;

    /// <summary>The current structure's terms for a prior the g-form's closed form does not serve: the normalized ridge's
    /// <see cref="Spectrum"/>, or the <see cref="Pencil"/> of a prior apart from the Gram; null otherwise.</summary>
    public RidgeColumns? PriorTerms => _outputs is not null ? null : _priorColumns is not null ? Pencil : _ridge ? Spectrum : null;

    /// <summary>The current structure's pencil, under weights with the prior free of them; null otherwise.</summary>
    public PriorPencil? Pencil
    {
        get
        {
            if (_priorColumns is null) return null;
            if (_pencil is { } cached) return cached;
            return _pencil = PencilOf(CurrentGram(), _priorGram!, CurrentShifts(), [.. _b], _response)
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
        if (Evaluate(term, _response) is not var (column, output, prior)) return DesignProposal.Refused(BassRefusal.Support);
        EnsureCapacity(_terms.Count + 1);
        double[]? tabled = TableBorder(term, column, out double diagonal, out bool approximate);
        double[] border = tabled ?? Border(column, output, -1, out diagonal);
        double bEntry = _outputs is null ? CrossResponse(column, output) : 0.0;   // several outputs' crosses are Decide's
        double crossResidual = _outputs is null ? CrossResidual(column, output) : 0.0;
        FactorProposal p = _factor.ProposeBirth(border, diagonal);
        DesignProposal? account = Decide(Move.Birth, _terms.Count, term, column, output, border, border, diagonal, bEntry, p,
                                         ResidualSquares, crossResidual, prior, approximate);
        if (account is null)   // near the gate or the explicit residual: the dots decide
        {
            TableFallbacks++;
            border = Border(column, output, -1, out diagonal);
            p = _factor.ProposeBirth(border, diagonal);
            account = Decide(Move.Birth, _terms.Count, term, column, output, border, border, diagonal, bEntry, p,
                             ResidualSquares, crossResidual, prior);
        }
        return account.Value;
    }

    /// <summary>Proposes removing term <paramref name="index"/>. A death never leaves F, since the gate is downward closed.</summary>
    public DesignProposal ProposeDeath(int index)
    {
        if ((uint)index >= (uint)_terms.Count) throw new ArgumentOutOfRangeException(nameof(index));
        _pending = null;
        FactorProposal p = _factor.ProposeDeath(index);
        if (_outputs is { } outputs)
        {
            // Several outputs: R′ = R + β̂_m β̂_mᵀ/(G⁻¹)_mm, and F′ from the reduced factor's whitened crosses.
            double inverseDiagonal = _factor.InverseDiagonal(index);
            double[] beta = outputs.CoefficientsOf(index);
            int k = outputs.Count;
            var residual = (double[,])outputs.Residual.Clone();
            for (int a = 0; a < k; a++)
                for (int d = 0; d < k; d++) residual[a, d] += beta[a] * beta[d] / inverseDiagonal;
            double[][] crosses = outputs.CandidateCross(index, null);
            var forms = new OutputForms(residual, ProposedFitted(crosses, null));
            double[,]? ownGram = null;
            PriorPencil? ownPencil = null;
            if (_priorColumns is not null)
            {
                ownGram = PriorGramWithout(index);
                ownPencil = PencilOf(Without(CurrentGram(), index), ownGram, Without(CurrentShifts(), index), new double[_terms.Count - 1], _response)!;
                forms = forms with { Prior = OutputPencil.Of(ownPencil, crosses, OutputMeans, _intercept, _rows) };
            }
            var several = new DesignProposal(null, _terms.Count - 1, double.NaN, double.NaN, p.MinLeaveOneOutPivot, double.NaN, null, forms);
            _pending = new Pending(Move.Death, index, null, null, null, 0.0, null, several, PriorGram: ownGram, Pencil: ownPencil);
            return several;
        }
        double[] b = Without(_b, index);
        double q = _factor.ProposedInverseQuadraticForm(b);
        double rss = ResidualSquares + _beta[index] * _beta[index] / _factor.InverseDiagonal(index);
        if (_priorColumns is not null)
        {
            // Downward closed: a principal block of a positive definite prior Gram stays so.
            double[,] priorGram = PriorGramWithout(index);
            PriorPencil pencil = PencilOf(Without(CurrentGram(), index), priorGram, Without(CurrentShifts(), index), b, _response)!;
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
            var same = _outputs is null
                ? new DesignProposal(null, _terms.Count, ResidualSquares, FittedSquares, MinAngle, _factor.LeaveOneOutPivot(index), PriorTerms)
                : new DesignProposal(null, _terms.Count, double.NaN, double.NaN, MinAngle, _factor.LeaveOneOutPivot(index), null, OutputForms);
            _pending = new Pending(Move.Change, index, term, null, null, 0.0, null, same);
            return same;
        }
        if (_terms.Contains(term)) return DesignProposal.Refused(BassRefusal.Numerical);
        if (Evaluate(term, _response) is not var (column, output, prior)) return DesignProposal.Refused(BassRefusal.Support);
        double[]? tabled = TableBorder(term, column, out double diagonal, out bool approximate);
        double[] full = tabled ?? Border(column, output, -1, out diagonal);
        double bEntry = _outputs is null ? CrossResponse(column, output) : 0.0;   // several outputs' crosses are Decide's
        double crossResidual = _outputs is null ? CrossResidual(column, output) : 0.0;

        // The reduced model, term index removed: ‖r‖² gains β̂_m²/(G⁻¹)_mm, and zᵀr gains β̂_m·uᵀG⁻¹e_m/(G⁻¹)_mm.
        var unit = new double[_terms.Count];
        unit[index] = 1.0;
        _factor.SolveInPlace(unit);
        double inverseDiagonal = unit[index], beta = _outputs is null ? _beta[index] : 0.0;
        double rssReduced = ResidualSquares + beta * beta / inverseDiagonal;

        DesignProposal? account = ChangeAccount(index, term, column, output, prior, full, diagonal, bEntry, rssReduced, crossResidual, beta, inverseDiagonal, unit, approximate);
        if (account is null)   // near the gate or the explicit residual: the dots decide
        {
            TableFallbacks++;
            full = Border(column, output, -1, out diagonal);
            account = ChangeAccount(index, term, column, output, prior, full, diagonal, bEntry, rssReduced, crossResidual, beta, inverseDiagonal, unit, false);
        }
        return account.Value;
    }

    // A change's account from the candidate's border against every current column, `full`; null when the border came from tables
    // and the proposal needs the dots' exactness.
    private DesignProposal? ChangeAccount(int index, BassTerm term, double[] column, double[]? output, BassPriorColumn? prior, double[] full,
                                          double diagonal, double bEntry, double rssReduced, double crossResidual, double beta, double inverseDiagonal,
                                          double[] unit, bool approximate)
    {
        double zrReduced = crossResidual + beta * PairwiseSum.Dot(full, unit) / inverseDiagonal;
        double[] border = Without(full, index);
        FactorProposal p = _factor.ProposeChange(index, border, diagonal);
        return Decide(Move.Change, index, term, column, output, border, full, diagonal, bEntry, p, rssReduced, zrReduced, prior, approximate, unit, inverseDiagonal);
    }

    /// <summary>
    /// Proposes refitting the current structure on <paramref name="working"/>, its model's reading at another latent state
    /// (<see cref="WorkingResponse.Moved"/>); refused when the structure lies outside F under it. An offset alone keeps the columns,
    /// the factor and the gate, at O(NM) with O(M²) for the factor's solve and a prior's terms. Weights, or a random-effect block's
    /// whitening at another scale (<see cref="WorkingResponse.Whitening"/>), whiten every column again and factor their Gram strictly, at
    /// O(NM²), by the arithmetic of a design built on them and loaded with these terms. Weights that
    /// count rows do not move: the support threshold and τ₀ read their sum at resolution.
    /// </summary>
    public DesignProposal ProposeFollow(WorkingResponse working)
    {
        ArgumentNullException.ThrowIfNull(working);
        _pending = null;
        bool reweighted = Reweights(working);
        BassDesignResponse response = Read(working, reweighted ? null : _response);
        int m = _terms.Count;
        List<double[]> columns = _a;
        List<BassPriorColumn>? priors = _priorColumns;
        BorderedCholesky factor = _factor;
        double[,]? gram = null;
        if (reweighted)
        {
            columns = new List<double[]>(m);
            priors = _priorColumns is null ? null : new List<BassPriorColumn>(m);
            foreach (BassTerm term in _terms)
            {
                (double[] column, _, BassPriorColumn? prior) = Evaluate(term, response)!.Value;
                columns.Add(column);
                if (prior is not null) priors!.Add(prior);
            }
            gram = Gram(columns, []);
            factor = new BorderedCholesky(Math.Max(_factor.Capacity, m));
            factor.Decompose(gram);
            if (m > 0 && factor.FailedAt >= 0) return DesignProposal.Refused(BassRefusal.Numerical);
        }
        double minAngle = factor.MinLeaveOneOutPivot();
        if (minAngle < _rules.AngleFloor) return DesignProposal.Refused(BassRefusal.Numerical);
        if (_rules.MaxVarianceInflation is { } cap && minAngle * cap < 1.0) return DesignProposal.Refused(BassRefusal.VarianceInflation);

        var b = new double[m];
        for (int j = 0; j < m; j++) b[j] = Cross(columns[j], _khatriRao ? _h[j] : null, response.Y);
        RidgeSpectrum? spectrum = null;
        PriorPencil? pencil = null;
        if (_priorColumns is not null)
        {
            pencil = reweighted ? PencilOf(gram!, _priorGram!, [.. priors!.Select(c => c.Shift)], b, response) : Pencil!.WithCross(PencilCross(b, response));
            if (pencil is null) return DesignProposal.Refused(BassRefusal.Numerical);
        }
        else if (_ridge) spectrum = reweighted ? RidgeSpectrum.Of(gram!, b, _rows, working.Weighting.Mass) : Spectrum.WithCross(b);

        double[] beta = m == 0 ? [] : factor.Solve(b);
        var residual = new double[_rows];
        Fit(columns, _h, beta, residual, response.Y);
        double rss = _output.Dot(residual, residual), q = m == 0 ? 0.0 : PairwiseSum.Dot(b, beta);
        // Several outputs (E34's slice 4): every output read again on the new reading and refitted on these columns and factor.
        DesignOutputs? outputs = null;
        if (_outputs is not null)
        {
            outputs = new DesignOutputs(ReadOutputs(working, response.Roots));
            outputs.Load(columns);
            outputs.Refit(factor, columns);
        }
        var account = new DesignProposal(null, m, rss, q, minAngle, double.NaN, (RidgeColumns?)pencil ?? spectrum, outputs?.Forms);
        var followed = new Followed(response, reweighted ? columns : null, reweighted ? priors : null, reweighted ? factor : null, b, beta, residual, rss, q, spectrum, pencil,
                                    outputs);
        _pending = new Pending(Move.Follow, -1, null, null, null, 0.0, null, account, Response: followed);
        return account;
    }

    // Whether `working` moves the weights, refused when it is not this design's reading of the response.
    private bool Reweights(WorkingResponse working)
    {
        WorkingResponse current = _response.Working;
        if (!ReferenceEquals(working.Response, current.Response)) throw new ArgumentException("The working response is of another response.", nameof(working));
        if (working.Meaning != current.Meaning || working.PriorGram != current.PriorGram || !ReferenceEquals(working.Reference, current.Reference)
            || !ReferenceEquals(working.Output, current.Output))
            throw new ArgumentException("What a weight states, the prior's Gram, the reference and the inner product are the model's, which WorkingResponse.Moved keeps.", nameof(working));
        if ((working.Whitening is null) != (current.Whitening is null))
        {
            // Under weights every column is evaluated afresh on a reweighting, so a noise's random effects may pass between their integral at
            // b = 1 and their reading as an offset (E21's slice 4); at unit weights whitening is the design's layout.
            if (current.WeightArray is null || working.WeightArray is null)
                throw new ArgumentException("A design's layout depends on whether its columns are whitened, so a model with a random-effect block declares it from the start.", nameof(working));
            return true;
        }
        if (!ReferenceEquals(working.Whitening, current.Whitening)) return true;   // the block at another scale: every column whitened afresh
        double[]? now = current.WeightArray, next = working.WeightArray;
        if ((now is null) != (next is null)) throw new ArgumentException("A design's layout depends on whether it holds weights, so a model whose weights move declares them from the start.", nameof(working));
        if (ReferenceEquals(now, next) || now.AsSpan().SequenceEqual(next)) return false;
        if (current.Meaning == WeightMeaning.Replication) throw new NotSupportedException("Weights that count rows do not move: the support threshold and τ₀ read their sum at resolution.");
        return true;
    }

    /// <summary>Commits the last proposal, which must have been admitted.</summary>
    public void Accept()
    {
        Pending p = _pending ?? throw new InvalidOperationException("Nothing admitted is proposed.");
        _pending = null;
        if (p.Response is { } followed)
        {
            Commit(followed);
            return;
        }
        if (p.Kind == Move.Change && p.Column is null) return;   // a change to the same term
        if (p.Strict is not null) _factor.Decompose(p.Strict);
        else _factor.Accept();
        _spectrum = p.Fresh;
        if (p.Kind is Move.Death or Move.Change)
        {
            _terms.RemoveAt(p.Index);
            _a.RemoveAt(p.Index);
            DropTables(p.Index);
            if (_khatriRao) _h.RemoveAt(p.Index);
            _b.RemoveAt(p.Index);
            _priorColumns?.RemoveAt(p.Index);
        }
        if (p.Kind is Move.Birth or Move.Change)
        {
            _terms.Add(p.Term!);
            _a.Add(p.Column!);
            _tables.Add(null);
            if (_khatriRao) _h.Add(p.Output!);
            _b.Add(p.BEntry);
            if (p.Prior is { } prior) _priorColumns!.Add(prior);
        }
        _outputs?.Accept(p.Kind is Move.Death or Move.Change ? p.Index : -1, p.Kind is Move.Birth or Move.Change ? p.OutputCross : null);
        if (_priorColumns is not null)
        {
            _priorGram = p.PriorGram;
            _pencil = p.Pencil;
        }
        Refit();
        _outputPencil = p.Account.Outputs?.Prior;
    }

    private void Commit(Followed f)
    {
        _response = f.Response;
        if (f.Columns is { } columns)
        {
            _a.Clear(); _a.AddRange(columns);
            ResetTables();
            _factor = f.Factor!;
            if (_priorColumns is not null)
            {
                _priorColumns.Clear(); _priorColumns.AddRange(f.Priors!);
            }
        }
        _b.Clear(); _b.AddRange(f.B);
        _spectrum = f.Spectrum;
        _pencil = f.Pencil;
        _beta = f.Beta;
        _residual = f.Residual;
        ResidualSquares = f.ResidualSquares;
        FittedSquares = f.FittedSquares;
        if (f.Outputs is { } outputs)
        {
            _outputs = outputs;
            _outputPencil = null;
        }
    }

    /// <summary>
    /// Replaces the structure with <paramref name="terms"/>, built and factored strictly from scratch, when it lies in
    /// F; otherwise returns the refusal and leaves the design as it was. Starts and restores use it. With
    /// <paramref name="working"/> it is built on that working response, as <see cref="ProposeFollow"/> reads it; null keeps the design's.
    /// </summary>
    public BassRefusal? Load(IReadOnlyList<BassTerm> terms, WorkingResponse? working = null)
    {
        ArgumentNullException.ThrowIfNull(terms);
        _pending = null;
        BassDesignResponse response = working is null || ReferenceEquals(working, _response.Working) ? _response : Read(working, Reweights(working) ? null : _response);
        bool reread = _outputs is not null && !ReferenceEquals(response, _response);   // several outputs read on a new reading (E34's slice 4)
        if (terms.Count > _rules.MaxBasis) return BassRefusal.Cap;
        var seen = new HashSet<BassTerm>();
        var columns = new List<double[]>(terms.Count);
        var outputs = new List<double[]>(terms.Count);
        var priors = new List<BassPriorColumn>(terms.Count);
        foreach (BassTerm term in terms)
        {
            if (!seen.Add(term)) return BassRefusal.Numerical;
            if (Evaluate(term, response) is not var (column, output, prior)) return BassRefusal.Support;
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
        for (int j = 0; j < columns.Count; j++) b[j] = Cross(columns[j], _khatriRao ? outputs[j] : null, response.Y);
        double[,]? priorGram = null;
        PriorPencil? pencil = null;
        if (_priorColumns is not null)
        {
            priorGram = PriorGramOf(priors);
            double[] shifts = [.. priors.Select(c => c.Shift)];
            if (PencilOf(g, priorGram, shifts, b, response) is not { } admitted) return BassRefusal.Numerical;
            pencil = admitted;
        }

        _response = response;
        _factor = fresh;
        _spectrum = null;
        _terms.Clear(); _terms.AddRange(terms);
        _a.Clear(); _a.AddRange(columns);
        ResetTables();
        _h.Clear(); _h.AddRange(outputs);
        _b.Clear(); _b.AddRange(b);
        if (_priorColumns is not null)
        {
            _priorColumns.Clear(); _priorColumns.AddRange(priors);
            _priorGram = priorGram;
            _pencil = pencil;
        }
        if (reread) _outputs = new DesignOutputs(ReadOutputs(response.Working, response.Roots));
        _outputs?.Load(columns);
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
               _spectrum, _pencil, _response, _outputs?.Snapshot());

    /// <summary>Replaces the structure with <paramref name="snapshot"/>'s, every statistic as the design it was taken from
    /// held it, so the two continue bit for bit; its columns shared when it carries them, else re-evaluated.</summary>
    internal void Restore(BassDesignSnapshot snapshot)
    {
        _outputPencil = null;
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
        _response = snapshot.Response ?? _response;
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
            if (Evaluate(term, _response) is not var (column, output, prior)) throw new ArgumentException($"The snapshot holds {term}, which lacks support here.", nameof(snapshot));
            _terms.Add(term);
            _a.Add(column);
            if (_khatriRao) _h.Add(output!);
            _b.Add(CrossResponse(column, output));
            if (prior is not null) _priorColumns!.Add(prior);
        }
        ResetTables();
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
        if (_outputs is not null) _outputs = new DesignOutputs(ReadOutputs(_response.Working, _response.Roots));   // the restored reading's outputs (E34's slice 4)
        _outputs?.Restore(snapshot.Outputs ?? throw new ArgumentException("The snapshot holds one response's statistics, and the design several outputs'.", nameof(snapshot)));
    }

    /// <summary>Refactors the Gram strictly and refits; returns the drift the rebuild removed.</summary>
    public FactorDrift Rebuild()
    {
        _pending = null;
        FactorDrift drift = _factor.Rebuild();
        Refit();
        return drift;
    }

    // `border` is the new column's Gram entries against the candidate's other columns, `full` against every current one. With
    // `approximate`, the border holds table entries, within rounding of the dots': null is returned, and the dots decide, for a proposal
    // that is not clear of the gate (bordered minimum angle under NearGate floors, or not positive definite) or whose residual
    // explanation passes the explicit threshold.
    private DesignProposal? Decide(Move kind, int index, BassTerm term, double[] column, double[]? output, double[] border,
                                   double[] full, double diagonal, double bEntry, FactorProposal p, double rssReduced, double zrReduced,
                                   BassPriorColumn? prior = null, bool approximate = false, double[]? unit = null, double inverseDiagonal = double.NaN)
    {
        double floor = _rules.AngleFloor;
        if (approximate && (_outputs is not null || !(p.PositiveDefinite && p.MinLeaveOneOutPivot >= NearGate * floor))) return null;
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

        if (_outputs is not null)
        {
            (OutputForms forms, double[] cross) = OutputsAccount(kind, index, column, full, newest * diagonal, strict, border, diagonal, unit, inverseDiagonal);
            double[,]? ownGram = null;
            PriorPencil? ownPencil = null;
            if (_priorColumns is not null)
            {
                int skip = kind == Move.Change ? index : -1;
                ownGram = PriorGramWith(skip, prior!);
                double[] shifts = skip >= 0 ? Without(CurrentShifts(), skip) : CurrentShifts();
                Array.Resize(ref shifts, shifts.Length + 1);
                shifts[^1] = prior!.Shift;
                int columns = kind == Move.Birth ? _terms.Count + 1 : _terms.Count;
                ownPencil = PencilOf(strict ?? Candidate(skip, border, diagonal), ownGram, shifts, new double[columns], _response);
                if (ownPencil is null) return DesignProposal.Refused(BassRefusal.Numerical);
                forms = forms with { Prior = OutputPencil.Of(ownPencil, _outputs.CandidateCross(skip, cross), OutputMeans, _intercept, _rows) };
            }
            var several = new DesignProposal(null, kind == Move.Birth ? _terms.Count + 1 : _terms.Count, double.NaN, double.NaN, minAngle, newest, null, forms);
            _pending = new Pending(kind, index, term, column, output, 0.0, strict, several,
                                   Prior: prior, PriorGram: ownGram, Pencil: ownPencil, OutputCross: cross);
            return several;
        }

        double rss = rssReduced - zrReduced * zrReduced / (newest * diagonal);
        if (!(rss > ExplicitBelow * rssReduced))
        {
            if (approximate) return null;
            rss = ExplicitResidualSquares(kind == Move.Change ? index : -1, column, output, strict ?? Candidate(kind == Move.Change ? index : -1, border, diagonal), b);
        }
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
            pencil = PencilOf(strict ?? Candidate(skip, border, diagonal), priorGram, shifts, b, _response);
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

    // Several outputs' account of a birth or a change (E34): the new column's crosses with the centred outputs and with the residuals, then
    // R′ by the scalar bordered update carried to K × K, R′_kl = R̃_kl − z̃_k z̃_l/(newest·δ) on the reduced model's R̃ and z̃ (a change's gain
    // β̂_mβ̂_mᵀ/(G⁻¹)_mm and β̂_m·uᵀG⁻¹e_m/(G⁻¹)_mm), formed explicitly where a move explains nearly all of an output's residual; and F′ from
    // the candidate factor's whitened crosses.
    private (OutputForms Forms, double[] Cross) OutputsAccount(Move kind, int index, double[] column, double[] full, double pivot, double[,]? strict,
                                                               double[] border, double diagonal, double[]? unit, double inverseDiagonal)
    {
        DesignOutputs outputs = _outputs!;
        int k = outputs.Count;
        double[] cross = DesignOutputs.Dot(column, outputs.Centred);
        double[] zr = DesignOutputs.Dot(column, outputs.Residuals);
        var reduced = (double[,])outputs.Residual.Clone();
        if (kind == Move.Change)
        {
            double[] beta = outputs.CoefficientsOf(index);
            double along = PairwiseSum.Dot(full, unit!);
            for (int a = 0; a < k; a++)
            {
                zr[a] += beta[a] * along / inverseDiagonal;
                for (int d = 0; d < k; d++) reduced[a, d] += beta[a] * beta[d] / inverseDiagonal;
            }
        }
        var residual = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int d = 0; d < k; d++) residual[a, d] = reduced[a, d] - zr[a] * zr[d] / pivot;
        bool explicitly = false;
        for (int a = 0; a < k && !explicitly; a++) explicitly = !(residual[a, a] > ExplicitBelow * reduced[a, a]);
        int skip = kind == Move.Change ? index : -1;
        double[][] b = outputs.CandidateCross(skip, cross);
        if (explicitly) residual = ExplicitResidualForm(skip, column, strict ?? Candidate(skip, border, diagonal), b);
        return (new OutputForms(residual, ProposedFitted(b, strict)), cross);
    }

    // F′_kl = b′_kᵀG′⁻¹b′_l from the candidate factor's whitened crosses: the bordered proposal's, or the strict factor's where it decided.
    private double[,] ProposedFitted(double[][] b, double[,]? strict)
    {
        BorderedCholesky? s = strict is null ? null : StrictFactor(strict);
        int k = b.Length;
        var coordinates = new double[k][];
        for (int l = 0; l < k; l++)
        {
            var c = new double[b[l].Length];
            if (s is null) _factor.ProposedWhiten(b[l], c);
            else
            {
                b[l].CopyTo(c, 0);
                if (c.Length > 0) s.ForwardSolveInPlace(c);
            }
            coordinates[l] = c;
        }
        return DesignOutputs.Gram(coordinates);
    }

    // Several outputs' residual form from the candidate's explicit residuals, for a move that explains nearly all of one; zero where the
    // candidate's strict factor fails, as the scalar's is.
    private double[,] ExplicitResidualForm(int skip, double[] column, double[,] g, double[][] b)
    {
        BorderedCholesky s = StrictFactor(g);
        int k = b.Length;
        var columns = new List<double[]>();
        for (int j = 0; j < _terms.Count; j++)
            if (j != skip) columns.Add(_a[j]);
        columns.Add(column);
        var residuals = new double[k][];
        for (int l = 0; l < k; l++)
        {
            var r = new double[_rows];
            if (s.FailedAt < 0) Fit(columns, _h, s.Solve(b[l]), r, _outputs!.Centred[l]);
            residuals[l] = r;
        }
        return DesignOutputs.Gram(residuals);
    }

    // The term's columns, or null when it falls short of the support threshold in a part. Scalar: the centred case
    // column. Direct: the case-side and output-side factors as evaluated. Under weights: the whitened column over all N rows,
    // stacked for a direct response, and its prior side when the prior stands apart.
    private (double[] Column, double[]? Output, BassPriorColumn? Prior)? Evaluate(BassTerm term, BassDesignResponse response)
    {
        _basis.Validate(term);
        var a = new double[_cases];
        _basis.CaseColumn(term, _basis.Inputs.Cases, a);
        if (!Supported(a, _rules.SupportMinimum, _counted)) return null;
        // A random-effect block under weights reads the raw column: L for the fit, the column itself for the prior's side (E21's slice 4).
        if (!_functional && response.Working.Weighted is { } weighted) return Whitened(term, a, response, weighted);
        response.Working.Whitening?.Whiten(a);   // the support reads the column as evaluated, the fit its whitened image
        if (!_functional)
        {
            if (response.Roots is not null) return Whitened(term, a, response);
            double mean = PairwiseSum.Of(a) / _cases;
            ColumnOps.Subtract(a, mean);
            if (_priorColumns is null) return (a, null, null);
            if (_precision.Unit) return (a, null, new BassPriorColumn([], -mean));   // the raw coefficients' centre is 0
            // Unit weights on both sides: the prior's centre is the column's, so its side is the centred column at the builder's scale.
            double scale = _precision.Scale(term);
            var centred = (double[])a.Clone();
            for (int i = 0; i < centred.Length; i++) centred[i] *= scale;
            return (a, null, new BassPriorColumn(centred, scale * mean - mean));
        }
        var h = new double[_locations];
        _basis.OutputColumn(term, _basis.Inputs.Outputs!, h);
        if (!Supported(h, _rules.SupportMinimumFunctional ?? 0.0, null)) return null;
        if (response.Roots is null) return (a, h, null);
        var stacked = new double[_rows];
        for (int i = 0; i < _cases; i++) ColumnOps.Scale(stacked.AsSpan(i * _locations, _locations), a[i], h);
        return Whitened(term, stacked, response);
    }

    // W^½(z − z̄_w), in place, and z's prior side when the prior stands apart, at the builder's scale on `term`.
    private (double[] Column, double[]? Output, BassPriorColumn? Prior) Whitened(BassTerm term, double[] z, BassDesignResponse response)
    {
        double mean = PairwiseSum.Dot(response.Working.WeightArray!, z) / response.Working.Weighting.Mass;
        BassPriorColumn? prior = null;
        if (_precision.Unit) prior = new BassPriorColumn([], -mean);   // the raw coefficients' centre is 0
        else if (_priorColumns is not null)
        {
            double own = (_priorWeights is null ? PairwiseSum.Of(z) : PairwiseSum.Dot(_priorWeights, z)) / _priorMass;
            var centred = new double[_rows];
            ColumnOps.CentreInto(centred, z, own, _priorRoots);
            double scale = _precision.Scale(term);
            if (scale == 1.0) prior = new BassPriorColumn(centred, own - mean);
            else
            {
                for (int i = 0; i < centred.Length; i++) centred[i] *= scale;
                prior = new BassPriorColumn(centred, scale * own - mean);
            }
        }
        ColumnOps.CentreAndScale(z, mean, response.Roots!);
        return (z, null, prior);
    }

    // A random-effect block under weights (E21's slice 4): L·z − c·m, m = cᵀLz/cᵀc the GLS centring, and the prior's side the raw column's
    // against that same centring. Under a noise u is a latent and the coefficient prior belongs to the model given it, so the reading that
    // integrates u at b = 1 carries the prior the reading given u does, and a chain that follows between them keeps one prior Gram.
    private (double[] Column, double[]? Output, BassPriorColumn? Prior) Whitened(BassTerm term, double[] z, BassDesignResponse response, WeightedWhitening weighted)
    {
        double[] column = weighted.Apply(z);
        double mean = weighted.CentringOf(column);
        ReadOnlySpan<double> c = weighted.Direction;
        for (int i = 0; i < column.Length; i++) column[i] -= c[i] * mean;
        BassPriorColumn? prior = null;
        if (_precision.Unit) prior = new BassPriorColumn([], -mean);   // the raw coefficients' centre is 0
        else if (_priorColumns is not null)
        {
            double own = (_priorWeights is null ? PairwiseSum.Of(z) : PairwiseSum.Dot(_priorWeights, z)) / _priorMass;
            var centred = new double[_rows];
            ColumnOps.CentreInto(centred, z, own, _priorRoots);
            double scale = _precision.Scale(term);
            if (scale == 1.0) prior = new BassPriorColumn(centred, own - mean);
            else
            {
                for (int i = 0; i < centred.Length; i++) centred[i] *= scale;
                prior = new BassPriorColumn(centred, scale * own - mean);
            }
        }
        return (column, null, prior);
    }

    private bool Supported(double[] part, double minimum, double[]? counts)
    {
        if (_rules.Support == SupportKind.ColumnMass) return (counts is null ? PairwiseSum.Of(part) : PairwiseSum.Dot(counts, part)) >= minimum;
        if (counts is not null)
        {
            double counted = 0.0;
            for (int i = 0; i < part.Length; i++)
                if (part[i] != 0.0) counted += counts[i];
            return counted >= minimum;
        }
        return ColumnOps.NonzeroCount(part) >= minimum;
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
    // centring, α_u = α + δᵀβ: D_α = N and D_β = G_u under the g-form, D_α = 1 and D_β = diag(G_u)/N under the ridge, with Σw₀ for
    // N on a frozen W₀; D_α = 1, D_β = I and δ = −z̄_w under KHAOS's raw ridge, whose coefficients sit on the uncentred columns; the
    // g-form's on columns scaled by 1/g_m under its modified g-prior (the builder, BassPriorPrecision). A flat intercept has no prior
    // to couple, so it stands apart: G against D_β.
    private PriorPencil? PencilOf(double[,] g, double[,] u, double[] shifts, double[] b, BassDesignResponse response)
    {
        int m = g.GetLength(0);
        if (_precision.Gated && m > 0)
        {
            BorderedCholesky gate = StrictFactor(u);
            if (gate.FailedAt >= 0 || gate.MinLeaveOneOutPivot() < _rules.AngleFloor) return null;
        }
        bool inside = _intercept != GPriorCentre.Flat;
        int offset = inside ? 1 : 0, p = m + offset;
        var likelihood = new double[p, p];
        var prior = new double[p, p];
        double priorMass = _precision.Mass(_priorMass);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                likelihood[i + offset, j + offset] = g[i, j];
                prior[i + offset, j + offset] = _precision.Entry(u, i, j, _priorMass);
                if (inside) prior[i + offset, j + offset] += priorMass * shifts[i] * shifts[j];
            }
            if (inside) prior[0, i + offset] = prior[i + offset, 0] = priorMass * shifts[i];
        }
        if (inside)
        {
            likelihood[0, 0] = response.Working.Weighting.Mass;
            prior[0, 0] = priorMass;
        }
        return PriorPencil.Of(likelihood, prior, PencilCross(b, response), m, inside);
    }

    // The pencil's cross products: the intercept's Σw(ȳ_w − c) first when it sits inside the prior, then the columns' with ỹ.
    private double[] PencilCross(double[] b, BassDesignResponse response)
    {
        if (_outputs is not null) return new double[b.Length + (_intercept == GPriorCentre.Flat ? 0 : 1)];
        if (_intercept == GPriorCentre.Flat) return b;
        var cross = new double[b.Length + 1];
        cross[0] = response.Working.Weighting.Mass * (response.Working.ResponseMean - response.PriorLevel);
        b.CopyTo(cross, 1);
        return cross;
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

    // The border of a one-way hinge against the current columns, each entry the column's prefix table at the hinge's knot (built here the
    // first time it is needed) over c, the hinge's peak. The diagonal is the column's own dot. Null where tables do not serve this
    // proposal: a design too small or with too many variables for them to pay, a direct design, a prior apart from the Gram, weights that
    // count rows, degree three or more, a term of order two or more or with a categorical or functional part, or a variable past the budget.
    // `approximate` is true with a column to read from: the entries are the dots' within rounding, not in bits. Under weights the table is over
    // roots∘z, since Σ rootᵢ zᵢ = 0 drops the new column's centring.
    private double[]? TableBorder(BassTerm term, double[] column, out double diagonal, out bool approximate)
    {
        diagonal = 0.0;
        approximate = false;
        int degree = _basis.Degree;
        // Under a block the held columns are whitened, V^(−½)xⱼ, where the tables would sum them against the raw hinge: xⱼᵀV^(−½)h for xⱼᵀV⁻¹h.
        if (!TablesOn || _functional || _khatriRao || _priorColumns is not null || _counted is not null || IsWhitened() || degree is not (1 or 2)
            || term.Continuous.Count != 1 || term.Categorical.Count != 0 || term.Functional.Count != 0)
            return null;
        HingeFactor f = term.Continuous[0];
        if (f.Variable >= TableVariables() || OrderOf(f.Variable) is not { } order) return null;
        double sign = f.Sign, c = Hinge.Power(degree, (sign + 1.0) / 2.0 - sign * order.Knots[f.Knot]);
        int variables = _basis.Inputs.Cases.ContinuousCount;
        double[]? roots = _response.Roots;
        var u = new double[_terms.Count];
        for (int j = 0; j < u.Length; j++)
        {
            ColumnTables tables = _tables[j] ??= new ColumnTables(variables);
            double[]? table = f.Sign > 0 ? tables.Above[f.Variable] : tables.Below[f.Variable];
            if (table is null)
            {
                if (f.Sign > 0) tables.Above[f.Variable] = table = HingeScan.Above(_a[j], roots, order, degree);
                else tables.Below[f.Variable] = table = HingeScan.Below(_a[j], roots, order, degree);
                TableBuilds++;
            }
            u[j] = c == 0.0 ? table[f.Knot] : table[f.Knot] / c;
            TableLookups++;
            approximate = true;
        }
        diagonal = PairwiseSum.Dot(column, column);
        return u;
    }

    // The leading variables whose tables, for every column the cap allows and both signs, fit TableBytesAllowed: a function of the inputs
    // and the cap, not of what the design holds.
    private int TableVariables()
    {
        if (_tableVariables >= 0) return _tableVariables;
        long used = 0;
        int served = 0;
        for (; served < _basis.Inputs.Cases.ContinuousCount; served++)
        {
            long bytes = 16L * Math.Max(_rules.MaxBasis, 1) * _basis.Inputs.ContinuousKnots(served).Length;
            if (used + bytes > _tableBytesAllowed) break;
            used += bytes;
        }
        return _tableVariables = served;
    }

    // The cases of continuous variable `variable` in order, built when first read; null where its knots are not its values' distinct entries.
    private HingeOrder? OrderOf(int variable)
    {
        CaseInputs cases = _basis.Inputs.Cases;
        _orders ??= new HingeOrder?[cases.ContinuousCount];
        _ordersTried ??= new bool[cases.ContinuousCount];
        if (!_ordersTried[variable])
        {
            _ordersTried[variable] = true;
            _orders[variable] = HingeOrder.Of(cases.Continuous(variable), _basis.Inputs.ContinuousKnots(variable));
        }
        return _orders[variable];
    }

    // The tables follow their columns: every column starts without one.
    private void ResetTables()
    {
        _tables.Clear();
        for (int j = 0; j < _a.Count; j++) _tables.Add(null);
    }

    private void DropTables(int index) => _tables.RemoveAt(index);

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

    private double CrossResponse(double[] column, double[]? output) => Cross(column, output, _response.Y);

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
        Fit(a, h, beta, r, _response.Y);
        return _output.Dot(r, r);
    }

    // r = y − Z_cβ, into r, for y the design's ỹ on one working response.
    // Element-wise, so the passes are ColumnOps': a scalar design's columns subtract in blocks of four, each entry meeting its columns
    // in column order, and a direct design's inner loop runs over the case's locations.
    private void Fit(List<double[]> a, List<double[]> h, double[] beta, double[] r, double[] y)
    {
        Array.Copy(y, r, _rows);
        if (!_khatriRao)
        {
            ColumnOps.SubtractScaledColumns(r, System.Runtime.InteropServices.CollectionsMarshal.AsSpan(a)[..beta.Length], beta);
            return;
        }
        for (int j = 0; j < beta.Length; j++)
        {
            double bj = beta[j];
            double[] aj = a[j], hj = h[j];
            double level = PairwiseSum.Of(aj) / _cases * (PairwiseSum.Of(hj) / _locations);
            for (int i = 0; i < _cases; i++) ColumnOps.SubtractScaledProduct(r.AsSpan(i * _locations, _locations), bj, aj[i], hj, level);
        }
    }

    private void Refit()
    {
        _outputPencil = null;
        if (_outputs is not null)
        {
            // Several outputs refit their own side; the scalar statistics stay the empty response's, at the structure's length.
            _outputs.Refit(_factor, _a);
            _beta = new double[_terms.Count];
            return;
        }
        double[] b = [.. _b];
        _beta = _terms.Count == 0 ? [] : _factor.Solve(b);
        Fit(_a, _h, _beta, _residual, _response.Y);
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
