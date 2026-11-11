using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Conformal;

/// <summary>
/// One stratum's training units against batches of candidates under Hoff's leave-one-out score (E16 slice 12a; Hoff 2023), the full
/// step's opt-in (<see cref="FullStepScore.LeaveOneOut"/>). Unit i's score is cᵢ = p(zᵢ | every other unit, the candidate included) and
/// the candidate's is A = p(z_new | the n units), read from the n-unit draws (<c>hoffScore_none</c>, <c>hoffScore_some</c>):
/// <list type="bullet">
/// <item><b>The candidate</b> scores ln A = ln Σ_s W̄_s p(z_new | θ_s), the draws' plain predictive.</item>
/// <item><b>A data unit</b> scores ln cᵢ = ln A − ln Bᵢ with Bᵢ = Σ_s W̄_s p(z_new | θ_s)/p(zᵢ | θ_s), the posterior given every unit but i
/// being the n-unit posterior times that ratio (<c>postUExcept_some</c>). With the add-one-in weights w̃_s ∝ W̄_s p(z_new | θ_s) it is
/// ln cᵢ = −ln Σ_s w̃_s/p(zᵢ | θ_s), one (n × S)·(S × G) product of the units' inverse likelihoods, each scaled by its largest, against
/// the weights.</item>
/// <item><b>The rule</b> needs only the Bᵢ: the candidate is accepted iff at least k units have Bᵢ ≥ 1 (<c>hoff_accepts_iff</c>), ties
/// counted for it (CB-1, <c>hoff_coverage</c>). A unit within 10⁻⁹ of the candidate is read again as Bᵢ = Σ_s W_s e^(ℓ_new,s − ℓᵢ,s)/Σ_s W_s
/// in draw order, so a unit whose column equals the candidate's has Bᵢ = 1 exactly and ties.</item>
/// <item><b>The error</b> (<c>coverage_close_shift</c>, <c>coverage_inflated_shift</c>). Â is common to every unit of an ordering, so the
/// ranks ignore it: a data unit's bound is ln B̂ᵢ's Monte Carlo error, from the relative variance Σ_s (w̃_s/(pᵢₛHᵢ) − W̄_s)² over the
/// candidate's relative efficiency along the chains, and the candidate's position bound is 0. Two more products give it.</item>
/// <item><b>Diagnostics</b> per grid value: the least Kish ESS of the data units' weights W̄_s p(z_new | θ_s)/p(zᵢ | θ_s); and, asked for,
/// where the inflated rule accepts, the largest k̂ among the units whose side of the candidate their own bound could change, |ln cᵢ − ln A|
/// within it, −∞ where no unit's is. Those are the units whose weights decide the set. A k̂ for every unit costs n PSIS fits per grid value,
/// fifteen times the in-sample step on BASS fits of 99 units; over the deciding units it is still most of the step where the weights
/// degenerate, so it is an opt-in.</item>
/// </list>
/// A unit with zero density under a draw of positive weight is refused: its weight is unbounded, and the n-unit posterior gives that
/// draw no mass. Not thread-safe: one per stratum and call.
/// </summary>
internal sealed class LeaveOneOutScorer
{
    private const double Floor = 1e-100;   // a product entry below it is summed again in logs, as the in-sample path's errors are
    private const double TieWindow = 1e-9;
    private readonly double[,] _ll;
    private readonly int[] _rows, _members;
    private readonly double[]? _logWeights;
    private readonly IReadOnlyList<int>? _chains;
    private readonly double[] _inverse, _inverseSquared, _bottom, _drawWeight;   // Q, Q² (units × draws), each unit's smallest log term, W̄
    private readonly double _weightSquares;   // Σ_s W̄_s²
    private readonly int _s, _level;
    private readonly InflationBounds _inflation;
    private double[] _w = [], _w2 = [], _wW = [], _g1 = [], _g2 = [], _g3 = [], _p = [];
    private PackedColumns? _right, _rightSquared, _rightWeighted;

    public LeaveOneOutScorer(double[,] ll, int[] rows, double[]? logWeights, int[] members, IReadOnlyList<int>? chains, int level, double multiple,
                             InflationBounds inflation = InflationBounds.PerPosition)
    {
        _ll = ll;
        _rows = rows;
        _logWeights = logWeights;
        _members = members;
        _chains = chains;
        _level = level;
        _inflation = inflation;
        Multiple = multiple;
        _s = rows.Length;
        int n = members.Length, s = _s;
        _drawWeight = new double[s];
        double squares = 0.0;
        for (int d = 0; d < s; d++)
        {
            _drawWeight[d] = logWeights is null ? 1.0 / s : Math.Exp(logWeights[d]);
            squares += _drawWeight[d] * _drawWeight[d];
        }
        _weightSquares = squares;
        _bottom = new double[n];
        Array.Fill(_bottom, double.PositiveInfinity);
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++)
            {
                double v = ll[rows[d], members[i]];
                if (double.IsNaN(v) || double.IsPositiveInfinity(v)) throw new ArgumentException($"Unit {members[i]}'s log term is NaN or +∞ under draw {rows[d]}.", "training");
                if (double.IsNegativeInfinity(v))
                    throw new ArgumentException($"Unit {members[i]} has zero density under draw {rows[d]}: its leave-one-out weight p(z_new | θ)/p(zᵢ | θ) is "
                                                + "unbounded there, and the n-unit posterior gives that draw no mass, so the leave-one-out score cannot be read "
                                                + "from these draws.", "training");
                if (v < _bottom[i]) _bottom[i] = v;
            }
        _inverse = new double[(long)n * s];
        _inverseSquared = new double[(long)n * s];
        for (int i = 0; i < n; i++)
        {
            int u = members[i];
            for (int d = 0; d < s; d++)
            {
                double q = Math.Exp(_bottom[i] - ll[rows[d], u]);
                _inverse[(long)i * s + d] = q;
                _inverseSquared[(long)i * s + d] = q * q;
            }
        }
    }

    /// <summary>z: each bound is this many Monte Carlo standard errors.</summary>
    public double Multiple { get; }

    /// <summary>Whether the readings take k̂ over the deciding units (<see cref="ConformalOptions.LeaveOneOutParetoK"/>); NaN, not read,
    /// where off.</summary>
    internal bool ReadParetoK { get; init; }

    /// <summary>The last batch's unit scores, units × candidates, ln cᵢ; NaN where the candidate is outside every draw's support.</summary>
    public double[] UnitScores { get; private set; } = [];

    /// <summary>The last batch's unit errors, units × candidates: one Monte Carlo standard error of ln B̂ᵢ, on the log scale.</summary>
    public double[] UnitErrors { get; private set; } = [];

    /// <summary>The readings of <paramref name="count"/> candidates from column <paramref name="first"/> into <paramref name="into"/> from
    /// its start: the counts, ranks, bounds and diagnostics.</summary>
    public void Rank(double[,] candidates, int first, int count, ConformalBayes.Reading into)
    {
        int n = _members.Length, s = _s, at = into.Start, size = s * count;
        if (_w.Length < size)
        {
            _w = new double[size];
            _w2 = new double[size];
            _wW = new double[size];
            _p = new double[size];
        }
        Span<double> own = into.Own.AsSpan(at, count), efficiency = into.Efficiency.AsSpan(at, count);
        var top = new double[count];
        Array.Fill(top, double.NegativeInfinity);
        for (int d = 0; d < s; d++)
            for (int b = 0; b < count; b++) top[b] = Math.Max(top[b], candidates[_rows[d], first + b]);
        // p = exp(ln p − its largest); W̄·p is the add-one-in weights' numerator, and its sum A·e^(−top).
        var totals = new double[count];
        for (int d = 0; d < s; d++)
            for (int b = 0; b < count; b++)
            {
                double p = double.IsNegativeInfinity(top[b]) ? 0.0 : Math.Exp(candidates[_rows[d], first + b] - top[b]);
                _p[d * count + b] = p;
                totals[b] += _drawWeight[d] * p;
            }
        for (int b = 0; b < count; b++)
        {
            bool outside = double.IsNegativeInfinity(top[b]) || !(totals[b] > 0.0);
            own[b] = outside ? double.NegativeInfinity : top[b] + Math.Log(totals[b]);
            into.LogPredictive[at + b] = own[b];
            efficiency[b] = 1.0;
        }
        for (int d = 0; d < s; d++)
            for (int b = 0; b < count; b++)
            {
                int k = d * count + b;
                double w = totals[b] > 0.0 ? _drawWeight[d] * _p[k] / totals[b] : 0.0;
                _w[k] = w;
                _w2[k] = w * w;
                _wW[k] = w * _drawWeight[d];
            }
        if (_chains is not null)
            for (int b = 0; b < count; b++)
            {
                if (double.IsNegativeInfinity(own[b])) continue;
                var values = new double[s];
                for (int d = 0; d < s; d++) values[d] = _p[d * count + b];
                efficiency[b] = PredictiveChecks.RelativeEfficiency(values, _chains);
            }

        if (UnitScores.Length < n * count)
        {
            UnitScores = new double[n * count];
            UnitErrors = new double[n * count];
            _g1 = new double[n * count];
            _g2 = new double[n * count];
            _g3 = new double[n * count];
        }
        if (n > 0)
        {
            _right = Fill(_right, _w, s, count);
            _rightSquared = Fill(_rightSquared, _w2, s, count);
            _rightWeighted = Fill(_rightWeighted, _wW, s, count);
            DenseProduct.Multiply(_inverse, n, _right, _g1);
            DenseProduct.Multiply(_inverseSquared, n, _rightSquared, _g2);
            DenseProduct.Multiply(_inverse, n, _rightWeighted, _g3);
        }

        Span<int> below = into.Below.AsSpan(at, count), tied = into.Tied.AsSpan(at, count), inflated = into.InflatedRanks.AsSpan(at, count);
        Span<double> ess = into.Ess.AsSpan(at, count), uniform = into.UniformErrors.AsSpan(at, count), ownError = into.CandidateErrors.AsSpan(at, count);
        below.Clear();
        tied.Clear();
        for (int b = 0; b < count; b++)
        {
            inflated[b] = 1;   // the candidate itself
            ownError[b] = 0.0;   // Â is common to every unit, so the candidate's position carries no bound
            uniform[b] = 0.0;
            ess[b] = double.IsNegativeInfinity(own[b]) ? 0.0 : double.PositiveInfinity;
        }
        for (int i = 0; i < n; i++)
            for (int b = 0; b < count; b++)
            {
                int k = i * count + b;
                if (double.IsNegativeInfinity(own[b]))
                {
                    UnitScores[k] = UnitErrors[k] = double.NaN;
                    continue;
                }
                double g1 = _g1[k], score, relative, unitEss;
                if (g1 >= Floor)
                {
                    score = _bottom[i] - Math.Log(g1);
                    relative = _g2[k] / (g1 * g1) - 2.0 * _g3[k] / g1 + _weightSquares;
                    unitEss = g1 * g1 / _g2[k];
                }
                else (score, relative, unitEss) = InLogs(candidates, first + b, i, own[b]);
                if (Math.Abs(score - own[b]) <= TieWindow * Math.Max(1.0, Math.Abs(own[b])))
                    score = own[b] - ExactLogB(candidates, first + b, _members[i]);
                double error = double.IsNaN(relative) ? double.PositiveInfinity
                             : PredictiveChecks.LogScaleError(Math.Max(0.0, relative), efficiency[b]);
                UnitScores[k] = score;
                UnitErrors[k] = error;
                double bound = Multiple * error;
                if (score < own[b]) below[b]++;
                else if (score == own[b]) tied[b]++;
                if (score - bound <= own[b]) inflated[b]++;
                if (bound > uniform[b]) uniform[b] = bound;
                if (unitEss < ess[b]) ess[b] = unitEss;
            }
        for (int b = 0; b < count; b++)
        {
            if (double.IsNegativeInfinity(own[b]))
            {
                into.ParetoK[at + b] = double.NaN;
                continue;
            }
            if (n == 0) ess[b] = 0.0;
            if (_inflation == InflationBounds.Uniform)
            {
                // One bound for every position, the largest.
                int rank = 1;
                for (int i = 0; i < n; i++)
                    if (UnitScores[i * count + b] - uniform[b] <= own[b] + uniform[b]) rank++;
                inflated[b] = rank;
            }
            into.ParetoK[at + b] = inflated[b] > _level && ReadParetoK ? DecidingParetoK(candidates, first + b, b, count, own[b], efficiency[b]) : double.NaN;
        }
    }

    private static PackedColumns Fill(PackedColumns? packed, double[] values, int s, int count)
    {
        if (packed is null || packed.Columns != count) packed = new PackedColumns(s, count);
        packed.Fill(values.AsSpan(0, s * count));
        return packed;
    }

    // Unit i's score, its relative variance and its weights' ESS, summed in logs where the product's entry underflows.
    private (double Score, double Relative, double Ess) InLogs(double[,] candidates, int column, int i, double logA)
    {
        int unit = _members[i];
        double LogWeight(int d) => (_logWeights is null ? -Math.Log(_s) : _logWeights[d]) + candidates[_rows[d], column] - logA;   // ln w̃
        double max = double.NegativeInfinity;
        for (int d = 0; d < _s; d++) max = Math.Max(max, LogWeight(d) - _ll[_rows[d], unit]);
        if (double.IsNegativeInfinity(max)) return (double.PositiveInfinity, double.NaN, 0.0);
        double sum = 0.0, squares = 0.0;
        for (int d = 0; d < _s; d++)
        {
            double e = Math.Exp(LogWeight(d) - _ll[_rows[d], unit] - max);
            sum += e;
            squares += e * e;
        }
        double logH = max + Math.Log(sum), relative = 0.0;
        for (int d = 0; d < _s; d++)
        {
            double psi = Math.Exp(LogWeight(d) - _ll[_rows[d], unit] - logH) - _drawWeight[d];
            relative += psi * psi;
        }
        return (-logH, relative, sum * sum / squares);
    }

    // ln Bᵢ = ln Σ_s W_s e^(ℓ_new,s − ℓᵢ,s) − ln Σ_s W_s, summed in draw order about its largest term, with W_s = 1 at equal weights: a
    // column equal to the candidate's gives every ratio 1, so the two sums agree and ln Bᵢ is 0 exactly.
    private double ExactLogB(double[,] candidates, int column, int unit)
    {
        double max = double.NegativeInfinity;
        for (int d = 0; d < _s; d++) max = Math.Max(max, Weight(d) + (candidates[_rows[d], column] - _ll[_rows[d], unit]));
        double numerator = 0.0, denominator = 0.0;
        for (int d = 0; d < _s; d++)
        {
            // The difference first: an equal column gives exactly 0, so each term is the denominator's.
            numerator += Math.Exp(Weight(d) + (candidates[_rows[d], column] - _ll[_rows[d], unit]) - max);
            denominator += Math.Exp(Weight(d) - max);
        }
        return Math.Log(numerator) - Math.Log(denominator);
    }

    private double Weight(int d) => _logWeights is null ? 0.0 : _logWeights[d];

    // The largest k̂ of the weights ln W̄_s + ln p(z_new | θ_s) − ln p(zᵢ | θ_s), PSIS's diagnostic for each unit's ratio, among the units
    // within their own bound of the candidate (an unbounded error counts); −∞ where none is.
    private double DecidingParetoK(double[,] candidates, int column, int b, int count, double own, double efficiency)
    {
        double largest = double.NegativeInfinity;
        double[]? ratios = null;
        for (int i = 0; i < _members.Length; i++)
        {
            int k = i * count + b;
            if (!(Math.Abs(UnitScores[k] - own) <= Multiple * UnitErrors[k])) continue;
            if (_s < 2) return double.PositiveInfinity;
            ratios ??= new double[_s];
            int unit = _members[i];
            for (int d = 0; d < _s; d++) ratios[d] = Weight(d) + candidates[_rows[d], column] - _ll[_rows[d], unit];
            largest = Math.Max(largest, PredictiveChecks.ParetoK(ratios, efficiency));
        }
        return largest;
    }
}
