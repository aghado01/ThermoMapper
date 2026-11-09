using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// E4's residual birth (the proposal-kernels brief's slice 2), one component of the birth mixture (<see cref="BassBirthMixture"/>):
/// a parent, then one new factor chosen by its fit to the structure's residual.
/// <list type="bullet">
/// <item><b>The residual</b> r(S) is the least-squares residual of the raw response on [1, B_S], the structure's raw columns in a
/// canonical order (<see cref="Residual"/>). It is a function of the terms alone, computed the same way for a birth's structure and a
/// death's reduced one, so a birth's forward probability and its reverse read one number.</item>
/// <item><b>The parent</b> p is drawn uniformly over P(S): the intercept and the structure's terms that can take one more factor of
/// the main part (the continuous hinges, or the Legendre factors) within the cap and, on a Legendre basis, the degree.</item>
/// <item><b>The new factor</b> f is drawn ∝ exp(½·b·Δℓ(f)), Δℓ(f) = −(n/2)·ln(1 − ΔQ(f)/‖r‖²), ΔQ(f) = (c̃ᵀr)²/‖c̃‖², c̃ the centred
/// column B_p ∘ φ_f: every variable p leaves, with every sign and knot of the prior's (a hinge, scored for every knot at once by
/// cumulative sums along the variable's sort order) or every degree within d_max − d_p (a Legendre factor). The exponent ½ is the
/// locally balanced √ weighting (Zanella 2020), and b the move's inverse temperature.</item>
/// <item><b>The forward probability,</b> exactly: q_R(t | S) = Σ_{p ∈ sub(t) ∩ P(S)} (1/|P(S)|)·w_p(t ∖ p)/W_p, sub(t) the terms
/// that t less one main-part factor leaves (the intercept for a one-factor term).</item>
/// </list>
/// An exact fit's residual is rounding: at or below 10⁻²⁴ of the response's own spread it weighs every candidate alike. Categorical and
/// functional factors are the base birth's to propose. A scalar response only: a direct functional response's or several outputs'
/// residual is a matrix, refused (<see cref="Refusal"/>).
/// <para>Its caches are pure, so they change no draw: the residual by structure and each parent's ln W_p by structure, parent and b,
/// with the last full table for draws, per instance; the raw columns per model, shared. A probability then needs one variable's pass.
/// One instance serves one chain or particle, not several at once.</para>
/// </summary>
internal sealed class BassResidualBirth
{
    private const int ResidualCacheLimit = 8, TotalCacheLimit = 4096, ColumnCacheLimit = 4096;

    // The model's tables, built once and shared by every chain and particle of it.
    private static readonly ConditionalWeakTable<BassModel, Tables> Built = [];

    private readonly BassModel _model;
    private readonly bool _legendre;
    private readonly int _vars, _cap, _maxDegree, _hingeDegree, _n;
    private readonly double[] _y;
    private readonly double _exactFit;
    private readonly int[][] _order;               // per variable, the cases by ascending value (hinge)
    private readonly double[][] _sorted;           // per variable, the sorted values (hinge)
    private readonly double[][] _knots;            // per variable, the prior's knots (hinge)
    private readonly int[] _offset;                // per variable, its first candidate slot
    private readonly int _slots;
    private readonly Dictionary<StructureKey, double[]> _residuals = [];
    private readonly Dictionary<(StructureKey, BassTerm?, double), double> _totals = [];   // ln W_p by structure, parent and b
    private readonly ConcurrentDictionary<BassTerm, double[]> _columns;                     // the model's, shared
    private (StructureKey Key, BassTerm? Parent, double B, double[] LogWeights)? _last;     // the last full table, for draws
    private double[] _scratch = [];

    private sealed record Table(double[] LogWeights, double LogTotal);

    // The immutable part: the raw response, each hinge variable's sort order, sorted values and knots, and the slots.
    private sealed class Tables
    {
        public Tables(BassModel model)
        {
            BassStructurePrior prior = model.Prior;
            Legendre = prior.Legendre is not null;
            BassPart part = Legendre ? BassPart.Legendre : BassPart.Continuous;
            Variables = prior.Variables(part);
            Cap = prior.Cap(part);
            MaxDegree = prior.Legendre?.MaxDegree ?? 0;
            Y = model.Data.Response.Values.ToArray();
            int n = Y.Length;
            double mean = 0.0, spread = 0.0;
            foreach (double v in Y) mean += v;
            mean /= n;
            foreach (double v in Y) spread += (v - mean) * (v - mean);
            ExactFit = 1e-24 * spread;
            Offset = new int[Variables + 1];
            Order = new int[Legendre ? 0 : Variables][];
            Sorted = new double[Legendre ? 0 : Variables][];
            Knots = new double[Legendre ? 0 : Variables][];
            CaseInputs cases = model.Inputs.Cases;
            for (int u = 0; u < Variables; u++)
            {
                if (Legendre)
                {
                    Offset[u + 1] = Offset[u] + MaxDegree;
                    continue;
                }
                double[] x = cases.Continuous(u).ToArray();
                int[] order = new int[n];
                for (int i = 0; i < n; i++) order[i] = i;
                Array.Sort((double[])x.Clone(), order);
                Order[u] = order;
                Sorted[u] = new double[n];
                for (int i = 0; i < n; i++) Sorted[u][i] = x[order[i]];
                Knots[u] = model.Inputs.ContinuousKnots(u).ToArray();
                Offset[u + 1] = Offset[u] + 2 * Knots[u].Length;
            }
        }

        public bool Legendre { get; }
        public int Variables { get; }
        public int Cap { get; }
        public int MaxDegree { get; }
        public double[] Y { get; }
        public ConcurrentDictionary<BassTerm, double[]> Columns { get; } = new();   // the raw columns, a pure cache every chain shares
        public double ExactFit { get; }   // ‖r‖² at or below which the residual is rounding: 10⁻²⁴ of the response's own spread
        public int[] Offset { get; }
        public int[][] Order { get; }
        public double[][] Sorted { get; }
        public double[][] Knots { get; }
    }

    /// <summary>The residual birth over <paramref name="model"/>'s main part; <see cref="Refusal"/> must be null.</summary>
    public BassResidualBirth(BassModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (Refusal(model) is { } refusal) throw new ArgumentException(refusal, nameof(model));
        _model = model;
        Tables tables = Built.GetValue(model, m => new Tables(m));
        _legendre = tables.Legendre;
        _vars = tables.Variables;
        _cap = tables.Cap;
        _maxDegree = tables.MaxDegree;
        _hingeDegree = model.Basis.Degree;
        _y = tables.Y;
        _exactFit = tables.ExactFit;
        _n = _y.Length;
        _offset = tables.Offset;
        _columns = tables.Columns;
        _order = tables.Order;
        _sorted = tables.Sorted;
        _knots = tables.Knots;
        _slots = _offset[_vars];
    }

    /// <summary>Why <paramref name="model"/> cannot run the residual birth, or null.</summary>
    public static string? Refusal(BassModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (model.Data.Response is not ScalarResponse)
            return "E4's residual birth reads a scalar response's residual; a direct functional response's, or several outputs', is a matrix, which no score here reads";
        if (!(model.Prior.Has(BassPart.Continuous) || model.Prior.Has(BassPart.Legendre)))
            return "E4's residual birth adds continuous factors, and the fit has no continuous part";
        return null;
    }

    // ---- the residual ------------------------------------------------------------------------------------------------------------

    /// <summary>r(S): the raw response's least-squares residual on [1, B_S], the columns in a canonical order, so it is a function of
    /// the terms alone. A tiny ridge, set by the Gram's own diagonal, keeps a collinear structure's solve defined.</summary>
    public double[] Residual(IReadOnlyList<BassTerm> terms)
    {
        ArgumentNullException.ThrowIfNull(terms);
        var key = new StructureKey(terms);
        if (_residuals.TryGetValue(key, out double[]? held)) return held;
        if (_residuals.Count >= ResidualCacheLimit) _residuals.Clear();
        return _residuals[key] = Fit(key.Sorted);
    }

    private double[] Fit(BassTerm[] sorted)
    {
        int n = _n, m = sorted.Length;
        double mean = 0.0;
        for (int i = 0; i < n; i++) mean += _y[i];
        mean /= n;
        var r = new double[n];
        for (int i = 0; i < n; i++) r[i] = _y[i] - mean;
        if (m == 0) return r;
        var z = new double[m][];
        for (int j = 0; j < m; j++)
        {
            double[] c = Column(sorted[j]);
            double cm = 0.0;
            for (int i = 0; i < n; i++) cm += c[i];
            cm /= n;
            var centred = new double[n];
            for (int i = 0; i < n; i++) centred[i] = c[i] - cm;
            z[j] = centred;
        }
        var g = new double[m, m];
        var b = new double[m];
        double top = 0.0;
        for (int j = 0; j < m; j++)
        {
            for (int k = 0; k <= j; k++)
            {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += z[j][i] * z[k][i];
                g[j, k] = g[k, j] = s;
            }
            double t = 0.0;
            for (int i = 0; i < n; i++) t += z[j][i] * r[i];
            b[j] = t;
            top = Math.Max(top, g[j, j]);
        }
        double[] beta = Solve(g, b, top);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++) r[i] -= beta[j] * z[j][i];
        return r;
    }

    // (G + λI)β = b by Cholesky, λ = 10⁻¹² of the largest diagonal, raised tenfold until the factor holds.
    private static double[] Solve(double[,] g, double[] b, double top)
    {
        int m = b.Length;
        double ridge = Math.Max(top, 1e-300) * 1e-12;
        var l = new double[m, m];
        for (int attempt = 0; ; attempt++)
        {
            bool held = true;
            for (int j = 0; j < m && held; j++)
            {
                double d = g[j, j] + ridge;
                for (int k = 0; k < j; k++) d -= l[j, k] * l[j, k];
                if (!(d > 0.0)) { held = false; break; }
                l[j, j] = Math.Sqrt(d);
                for (int i = j + 1; i < m; i++)
                {
                    double s = g[i, j];
                    for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                    l[i, j] = s / l[j, j];
                }
            }
            if (held) break;
            if (attempt > 12) throw new ArithmeticException("The residual's Gram did not factor.");
            ridge *= 10.0;
        }
        var w = new double[m];
        for (int i = 0; i < m; i++)
        {
            double s = b[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * w[k];
            w[i] = s / l[i, i];
        }
        var beta = new double[m];
        for (int i = m - 1; i >= 0; i--)
        {
            double s = w[i];
            for (int k = i + 1; k < m; k++) s -= l[k, i] * beta[k];
            beta[i] = s / l[i, i];
        }
        return beta;
    }

    private double[] Column(BassTerm term)
    {
        if (_columns.TryGetValue(term, out double[]? held)) return held;
        if (_columns.Count >= ColumnCacheLimit) _columns.Clear();
        var column = new double[_n];
        _model.Basis.CaseColumn(term, _model.Inputs.Cases, column);
        return _columns.GetOrAdd(term, column);
    }

    // ---- parents ------------------------------------------------------------------------------------------------------------------

    private int MainOrder(BassTerm term) => _legendre ? term.Legendre.Count : term.Continuous.Count;

    private bool Uses(BassTerm term, int variable)
    {
        if (_legendre)
        {
            foreach (LegendreFactor f in term.Legendre)
                if (f.Variable == variable) return true;
            return false;
        }
        foreach (HingeFactor h in term.Continuous)
            if (h.Variable == variable) return true;
        return false;
    }

    // Whether `term` (null, the intercept) can take one more main-part factor.
    private bool Eligible(BassTerm? term)
    {
        if (term is null) return _cap >= 1 && _vars >= 1 && (!_legendre || _maxDegree >= 1);
        int order = MainOrder(term);
        if (order >= _cap || order >= _vars) return false;
        return !_legendre || term.LegendreDegree < _maxDegree;
    }

    // P(S): the intercept, when eligible, then the structure's eligible terms in its order.
    private List<BassTerm?> Parents(IReadOnlyList<BassTerm> terms)
    {
        var parents = new List<BassTerm?>(terms.Count + 1);
        if (Eligible(null)) parents.Add(null);
        foreach (BassTerm t in terms)
            if (Eligible(t)) parents.Add(t);
        return parents;
    }

    // ---- the scores ----------------------------------------------------------------------------------------------------------------

    // The candidates' log weights ½·b·Δℓ from `parent` (null, the intercept) on the structure's residual, and their log total: the
    // last full table is kept for draws, and each total by its structure, parent and b, so a probability needs one variable's pass.
    private Table TableFor(StructureKey key, IReadOnlyList<BassTerm> terms, BassTerm? parent, double b)
    {
        if (_last is { } last && last.Key.Equals(key) && Equals(last.Parent, parent) && last.B == b)
            return new Table(last.LogWeights, _totals[(key, parent, b)]);
        var logw = new double[_slots];
        Array.Fill(logw, double.NegativeInfinity);
        for (int u = 0; u < _vars; u++) Score(u, terms, parent, b, logw);
        double max = double.NegativeInfinity;
        foreach (double v in logw) max = Math.Max(max, v);
        double total = 0.0;
        foreach (double v in logw)
            if (!double.IsNegativeInfinity(v)) total += Math.Exp(v - max);
        double logTotal = max + Math.Log(total);
        if (_totals.Count >= TotalCacheLimit) _totals.Clear();
        _totals[(key, parent, b)] = logTotal;
        _last = (key, parent, b, logw);
        return new Table(logw, logTotal);
    }

    // ln W_p, from the cache or a full table.
    private double LogTotal(StructureKey key, IReadOnlyList<BassTerm> terms, BassTerm? parent, double b)
        => _totals.TryGetValue((key, parent, b), out double held) ? held : TableFor(key, terms, parent, b).LogTotal;

    // One slot's log weight: its variable's pass alone, the same arithmetic as the full table's.
    private double LogWeightAt(IReadOnlyList<BassTerm> terms, BassTerm? parent, double b, int slot)
    {
        int u = 0;
        while (_offset[u + 1] <= slot) u++;
        if (_scratch.Length != _slots) _scratch = new double[_slots];
        Array.Fill(_scratch, double.NegativeInfinity, _offset[u], _offset[u + 1] - _offset[u]);
        Score(u, terms, parent, b, _scratch);
        return _scratch[slot];
    }

    // Variable u's slots from `parent` on the structure's residual; nothing for a variable the parent uses.
    private void Score(int u, IReadOnlyList<BassTerm> terms, BassTerm? parent, double b, double[] logw)
    {
        if (parent is not null && Uses(parent, u)) return;
        double[] r = Residual(terms);
        double rss = 0.0, sumR = 0.0;
        foreach (double v in r)
        {
            rss += v * v;
            sumR += v;
        }
        double[] parentColumn = parent is null ? Ones() : Column(parent);
        if (_legendre) ScoreLegendre(u, _maxDegree - (parent?.LegendreDegree ?? 0), parentColumn, r, rss, sumR, b, logw);
        else ScoreHinge(u, parentColumn, r, rss, sumR, b, logw);
    }

    private double[] _ones = [];
    private double[] _sums = [];

    private double[] Ones()
    {
        if (_ones.Length != _n)
        {
            _ones = new double[_n];
            Array.Fill(_ones, 1.0);
        }
        return _ones;
    }

    // ½·b·Δℓ from the candidate's cross product with r, its sum and its sum of squares.
    private double LogWeight(double cross, double sum, double squares, double rss, double sumR, double b)
    {
        if (!(rss > _exactFit)) return 0.0;   // an exact fit: the residual is rounding, and every candidate weighs alike
        double centredCross = cross - sum / _n * sumR;
        double centredSquares = squares - sum * sum / _n;
        if (!(centredSquares > 1e-13 * Math.Max(squares, 1e-300))) return 0.0;
        double share = centredCross * centredCross / centredSquares / rss;
        share = Math.Min(share, 1.0 - 1e-12);
        double gain = -0.5 * _n * Math.Log(1.0 - share);
        return 0.5 * b * gain;
    }

    // Every degree d ≤ budget of variable u: c = B_p ∘ ψ_d(x_u).
    private void ScoreLegendre(int u, int budget, double[] parent, double[] r, double rss, double sumR, double b, double[] logw)
    {
        if (budget < 1) return;
        ReadOnlySpan<double> x = _model.Inputs.Cases.Continuous(u);
        var cross = new double[budget + 1];
        var sum = new double[budget + 1];
        var squares = new double[budget + 1];
        Span<double> psi = stackalloc double[budget + 1];
        for (int i = 0; i < _n; i++)
        {
            if (parent[i] == 0.0) continue;
            ShiftedLegendre.Values(x[i], psi);
            for (int d = 1; d <= budget; d++)
            {
                double c = parent[i] * psi[d];
                cross[d] += c * r[i];
                sum[d] += c;
                squares[d] += c * c;
            }
        }
        for (int d = 1; d <= budget; d++) logw[_offset[u] + d - 1] = LogWeight(cross[d], sum[d], squares[d], rss, sumR, b);
    }

    // Every sign and knot of variable u's hinge: c = B_p ∘ (s(x_u − t))₊^q, its cross product, sum and sum of squares from the
    // sums of B·r·x^j, B·x^j (j ≤ q) and B²·x^j (j ≤ 2q) beyond the knot (s = +1) or below it (s = −1), by binomial expansion.
    private void ScoreHinge(int u, double[] parent, double[] r, double rss, double sumR, double b, double[] logw)
    {
        int q = _hingeDegree, n = _n;
        int[] order = _order[u];
        double[] xs = _sorted[u], knots = _knots[u];
        int k1 = q + 1, k2 = 2 * q + 1;
        // Suffix sums from position m (cases at sorted positions ≥ m) and prefix sums below m, each power; the buffers are the
        // instance's, every row written before it is read, the ends set to zero.
        if (_sums.Length != 6 * (n + 1) * k2) _sums = new double[6 * (n + 1) * k2];
        int block = (n + 1) * k2;
        Span<double> sufA = _sums.AsSpan(0, block), sufB = _sums.AsSpan(block, block), sufC = _sums.AsSpan(2 * block, block);
        Span<double> preA = _sums.AsSpan(3 * block, block), preB = _sums.AsSpan(4 * block, block), preC = _sums.AsSpan(5 * block, block);
        sufA.Slice(n * k1, k1).Clear();
        sufB.Slice(n * k1, k1).Clear();
        sufC.Slice(n * k2, k2).Clear();
        preA[..k1].Clear();
        preB[..k1].Clear();
        preC[..k2].Clear();
        Span<double> pw = stackalloc double[k2];
        for (int m = n - 1; m >= 0; m--)
        {
            int i = order[m];
            double bp = parent[i], br = bp * r[i], bb = bp * bp, x = xs[m];
            pw[0] = 1.0;
            for (int j = 1; j < k2; j++) pw[j] = pw[j - 1] * x;
            for (int j = 0; j < k1; j++)
            {
                sufA[m * k1 + j] = sufA[(m + 1) * k1 + j] + br * pw[j];
                sufB[m * k1 + j] = sufB[(m + 1) * k1 + j] + bp * pw[j];
            }
            for (int j = 0; j < k2; j++) sufC[m * k2 + j] = sufC[(m + 1) * k2 + j] + bb * pw[j];
        }
        for (int m = 0; m < n; m++)
        {
            int i = order[m];
            double bp = parent[i], br = bp * r[i], bb = bp * bp, x = xs[m];
            pw[0] = 1.0;
            for (int j = 1; j < k2; j++) pw[j] = pw[j - 1] * x;
            for (int j = 0; j < k1; j++)
            {
                preA[(m + 1) * k1 + j] = preA[m * k1 + j] + br * pw[j];
                preB[(m + 1) * k1 + j] = preB[m * k1 + j] + bp * pw[j];
            }
            for (int j = 0; j < k2; j++) preC[(m + 1) * k2 + j] = preC[m * k2 + j] + bb * pw[j];
        }
        int kn = knots.Length, above = 0, below = 0;
        for (int k = 0; k < kn; k++)
        {
            double t = knots[k];
            while (above < n && xs[above] <= t) above++;   // the first position with x > t
            while (below < n && xs[below] < t) below++;    // the first position with x ≥ t
            // s = +1: (x − t)^q = Σ_j C(q, j) x^j (−t)^(q−j), over positions ≥ above.
            double crossP = Expand(sufA, above * k1, q, -t), sumP = Expand(sufB, above * k1, q, -t), sqP = Expand(sufC, above * k2, 2 * q, -t);
            // s = −1: (t − x)^q = Σ_j C(q, j) (−x)^j t^(q−j), over positions < below.
            double crossM = ExpandMinus(preA, below * k1, q, t), sumM = ExpandMinus(preB, below * k1, q, t), sqM = ExpandMinus(preC, below * k2, 2 * q, t);
            logw[_offset[u] + k] = LogWeight(crossM, sumM, sqM, rss, sumR, b);            // sign −1
            logw[_offset[u] + kn + k] = LogWeight(crossP, sumP, sqP, rss, sumR, b);       // sign +1
        }
    }

    // Σ_j C(q, j)·S_j·a^(q−j): Σ of (x + a)^q weighted, from power sums S_j at `at`.
    private static double Expand(ReadOnlySpan<double> sums, int at, int q, double a)
    {
        double total = 0.0, binomial = 1.0;
        for (int j = 0; j <= q; j++)
        {
            total += binomial * sums[at + j] * Math.Pow(a, q - j);
            binomial = binomial * (q - j) / (j + 1);
        }
        return total;
    }

    // Σ_j C(q, j)·(−1)^j·S_j·t^(q−j): Σ of (t − x)^q weighted, from power sums S_j at `at`.
    private static double ExpandMinus(ReadOnlySpan<double> sums, int at, int q, double t)
    {
        double total = 0.0, binomial = 1.0;
        for (int j = 0; j <= q; j++)
        {
            total += binomial * (j % 2 == 0 ? 1.0 : -1.0) * sums[at + j] * Math.Pow(t, q - j);
            binomial = binomial * (q - j) / (j + 1);
        }
        return total;
    }

    // ---- the slots -----------------------------------------------------------------------------------------------------------------

    // The candidate slot of `factor`'s index into the main part: a hinge's variable, sign and knot, or a Legendre factor's variable
    // and degree.
    private int Slot(HingeFactor h) => _offset[h.Variable] + (h.Sign > 0 ? _knots[h.Variable].Length : 0) + h.Knot;

    private int Slot(LegendreFactor f) => _offset[f.Variable] + f.Degree - 1;

    // The term the slot adds to `parent` (null, the intercept).
    private BassTerm Child(BassTerm? parent, int slot)
    {
        int u = 0;
        while (_offset[u + 1] <= slot) u++;   // a variable with no slots shares its offset, and is passed over
        int within = slot - _offset[u];
        if (_legendre)
        {
            var factor = new LegendreFactor(u, within + 1);
            return new BassTerm(legendre: parent is null ? [factor] : [.. parent.Legendre, factor]);
        }
        int kn = _knots[u].Length;
        var hinge = within >= kn ? new HingeFactor(u, 1, within - kn) : new HingeFactor(u, -1, within);
        return parent is null ? new BassTerm(continuous: [hinge]) : new BassTerm([.. parent.Continuous, hinge], parent.Categorical, parent.Functional);
    }

    // ---- for the gates ------------------------------------------------------------------------------------------------------------

    /// <summary>The candidate slots: every variable's signs and knots, or degrees.</summary>
    internal int Slots => _slots;

    /// <summary>The term slot <paramref name="slot"/> adds to <paramref name="parent"/> (null, the intercept).</summary>
    internal BassTerm ChildOf(BassTerm? parent, int slot) => Child(parent, slot);

    /// <summary>P(S) for the structure <paramref name="terms"/>: the intercept (null) when eligible, then its eligible terms.</summary>
    internal IReadOnlyList<BassTerm?> ParentsOf(IReadOnlyList<BassTerm> terms) => Parents(terms);

    /// <summary>A copy of the candidates' log weights ½·b·Δℓ from <paramref name="parent"/> on <paramref name="terms"/>' residual,
    /// −∞ where the parent leaves no such candidate.</summary>
    internal double[] LogWeights(IReadOnlyList<BassTerm> terms, BassTerm? parent, double b)
        => (double[])TableFor(new StructureKey(terms), terms, parent, b).LogWeights.Clone();

    // ---- the draw and its probability ----------------------------------------------------------------------------------------------

    /// <summary>A birth into <paramref name="at"/>'s structure: a parent uniformly over P(S), then a new factor by its weight; null
    /// when no parent can take one.</summary>
    public BassTerm? Draw(Xoshiro256PlusPlus rng, BassBirthContext at)
    {
        ArgumentNullException.ThrowIfNull(rng);
        ArgumentNullException.ThrowIfNull(at);
        List<BassTerm?> parents = Parents(at.Terms);
        if (parents.Count == 0) return null;
        BassTerm? parent = parents[BassDraws.Below(rng, parents.Count)];
        var key = new StructureKey(at.Terms);
        Table table = TableFor(key, at.Terms, parent, at.InverseTemperature);
        double target = rng.NextDouble(), cumulative = 0.0;
        int last = -1;
        for (int s = 0; s < table.LogWeights.Length; s++)
        {
            double lw = table.LogWeights[s];
            if (double.IsNegativeInfinity(lw)) continue;
            last = s;
            cumulative += Math.Exp(lw - table.LogTotal);
            if (target < cumulative) return Child(parent, s);
        }
        return Child(parent, last);
    }

    /// <summary>ln q_R(t | S) for the structure <paramref name="at"/> holds: −∞ when no parent of P(S) is a sub-term of t.</summary>
    public double LogProbability(BassTerm term, BassBirthContext at)
    {
        ArgumentNullException.ThrowIfNull(term);
        ArgumentNullException.ThrowIfNull(at);
        List<BassTerm?> parents = Parents(at.Terms);
        if (parents.Count == 0) return double.NegativeInfinity;
        var held = new HashSet<BassTerm>(at.Terms);
        var key = new StructureKey(at.Terms);
        double logP = -Math.Log(parents.Count);
        double top = double.NegativeInfinity;
        Span<double> parts = stackalloc double[Math.Max(1, MainOrder(term))];
        int count = 0;
        if (_legendre)
        {
            for (int k = 0; k < term.Legendre.Count; k++)
            {
                LegendreFactor f = term.Legendre[k];
                BassTerm? parent = term.Legendre.Count == 1 ? null : new BassTerm(legendre: Without(term.Legendre, k));
                if (Contribution(parent, Slot(f), held, key, at) is double v) { parts[count++] = v; top = Math.Max(top, v); }
            }
        }
        else
        {
            int factors = term.Continuous.Count + term.Categorical.Count + term.Functional.Count;
            for (int k = 0; k < term.Continuous.Count; k++)
            {
                HingeFactor h = term.Continuous[k];
                BassTerm? parent = factors == 1 ? null : new BassTerm(Without(term.Continuous, k), term.Categorical, term.Functional);
                if (Contribution(parent, Slot(h), held, key, at) is double v) { parts[count++] = v; top = Math.Max(top, v); }
            }
        }
        if (count == 0) return double.NegativeInfinity;
        double sum = 0.0;
        for (int i = 0; i < count; i++) sum += Math.Exp(parts[i] - top);
        return logP + top + Math.Log(sum);
    }

    // ln[w_p(slot)/W_p] when `parent` is in P(S); null otherwise.
    private double? Contribution(BassTerm? parent, int slot, HashSet<BassTerm> held, StructureKey key, BassBirthContext at)
    {
        if (parent is not null && !held.Contains(parent)) return null;
        if (!Eligible(parent)) return null;
        double logTotal = LogTotal(key, at.Terms, parent, at.InverseTemperature);
        double lw = LogWeightAt(at.Terms, parent, at.InverseTemperature, slot);
        return double.IsNegativeInfinity(lw) ? null : lw - logTotal;
    }

    private static T[] Without<T>(IReadOnlyList<T> factors, int index)
    {
        var rest = new T[factors.Count - 1];
        for (int i = 0, j = 0; i < factors.Count; i++)
            if (i != index) rest[j++] = factors[i];
        return rest;
    }

    // ---- the structure's key -------------------------------------------------------------------------------------------------------

    // A structure as a set: its terms in a canonical order, equal as sequences.
    private sealed class StructureKey : IEquatable<StructureKey>
    {
        private readonly int _hash;

        public StructureKey(IReadOnlyList<BassTerm> terms)
        {
            Sorted = new BassTerm[terms.Count];
            for (int i = 0; i < Sorted.Length; i++) Sorted[i] = terms[i];
            Array.Sort(Sorted, TermOrder.Compare);
            var h = new HashCode();
            foreach (BassTerm t in Sorted) h.Add(t);
            _hash = h.ToHashCode();
        }

        public BassTerm[] Sorted { get; }

        public bool Equals(StructureKey? other)
        {
            if (other is null || other._hash != _hash || other.Sorted.Length != Sorted.Length) return false;
            for (int i = 0; i < Sorted.Length; i++)
                if (!Sorted[i].Equals(other.Sorted[i])) return false;
            return true;
        }

        public override bool Equals(object? obj) => Equals(obj as StructureKey);

        public override int GetHashCode() => _hash;
    }

    // A total order on terms by their factors: the Legendre, the continuous, the categorical, then the functional ones.
    private static class TermOrder
    {
        public static int Compare(BassTerm? a, BassTerm? b)
        {
            if (ReferenceEquals(a, b)) return 0;
            if (a is null) return -1;
            if (b is null) return 1;
            int c = CompareLists(a.Legendre, b.Legendre, (x, y) => x.Variable != y.Variable ? x.Variable.CompareTo(y.Variable) : x.Degree.CompareTo(y.Degree));
            if (c != 0) return c;
            c = CompareLists(a.Continuous, b.Continuous, CompareHinge);
            if (c != 0) return c;
            c = CompareLists(a.Categorical, b.Categorical, (x, y) =>
            {
                int v = x.Variable.CompareTo(y.Variable);
                return v != 0 ? v : CompareLists(x.Levels, y.Levels, (l, m) => l.CompareTo(m));
            });
            return c != 0 ? c : CompareLists(a.Functional, b.Functional, CompareHinge);
        }

        private static int CompareHinge(HingeFactor x, HingeFactor y)
            => x.Variable != y.Variable ? x.Variable.CompareTo(y.Variable) : x.Sign != y.Sign ? x.Sign.CompareTo(y.Sign) : x.Knot.CompareTo(y.Knot);

        private static int CompareLists<T>(IReadOnlyList<T> a, IReadOnlyList<T> b, Comparison<T> compare)
        {
            int c = a.Count.CompareTo(b.Count);
            if (c != 0) return c;
            for (int i = 0; i < a.Count; i++)
            {
                c = compare(a[i], b[i]);
                if (c != 0) return c;
            }
            return 0;
        }
    }
}
