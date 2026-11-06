using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using Maths.LinAlg;

namespace Maths.Regression.Conformal;

/// <summary>What the full step asks for (E16).</summary>
public sealed record ConformalOptions
{
    /// <summary>
    /// The miscoverage α, in (0, 1). A case's stratum of m units, its training units and the candidate, takes the level
    /// k = ⌊α·m⌋, and its set {y : π(y) > k/m} covers at least 1 − k/m ≥ 1 − α (CB-1). Below m = ⌈1/α⌉, k is 0 and the set is
    /// the whole space.
    /// </summary>
    public double Alpha { get; init; } = 0.1;
}

/// <summary>A run of accepted grid values, from its lowest to its highest; the set's boundary lies in the gaps beside it.</summary>
public readonly record struct GridInterval(double Lower, double Upper);

/// <summary>
/// One test case's conformal set on its grid (<see cref="ConformalBayes.FullStep"/>). The candidate at grid value y_g is ranked
/// among the <see cref="Units"/> − 1 training units of its <see cref="Stratum"/>: <see cref="Below"/> score strictly below it and
/// <see cref="Tied"/> equal to it, so π(y_g) = (Below + Tied + 1)/m, ties counted for the candidate (CB-1), and y_g is accepted
/// when π(y_g) > k/m. π need not be monotone in y, so the set can be a union of intervals (<see cref="Intervals"/>). A stratum
/// too small for the level (k = 0) has the whole space as its set, and nothing is computed: <see cref="WholeSpace"/>, with the
/// per-value arrays null.
/// </summary>
/// <param name="Case">The test case, in the candidates' order.</param>
/// <param name="Stratum">Its stratum.</param>
/// <param name="Units">m: the stratum's training units and the candidate.</param>
/// <param name="Level">k = ⌊α·m⌋.</param>
/// <param name="Grid">The case's candidate values, ascending.</param>
/// <param name="PValues">π(y_g) for each grid value; null for the whole space.</param>
/// <param name="Below">The training units of the stratum scoring strictly below the candidate at y_g.</param>
/// <param name="Tied">Those scoring equal to it; the strict rule (accept when Below ≥ k) and random tie-breaking read
/// these with <see cref="Below"/>.</param>
/// <param name="CandidateScores">The candidate's own log conformity score at y_g, ln p(y_g | all m units); −∞ outside every
/// draw's support.</param>
/// <param name="EffectiveSampleSizes">Kish's ESS of the add-one-in weights at y_g, 1/Σ w̃²; 0 outside every draw's
/// support.</param>
public sealed record ConformalSet(int Case, int Stratum, int Units, int Level, double[] Grid, double[]? PValues, int[]? Below,
                                  int[]? Tied, double[]? CandidateScores, double[]? EffectiveSampleSizes)
{
    /// <summary>Whether the stratum is too small for the level (k = 0), so the set is the whole space.</summary>
    public bool WholeSpace => PValues is null;

    /// <summary>The guaranteed coverage 1 − k/m: at least this, exactly this with distinct scores.</summary>
    public double Coverage => 1.0 - (double)Level / Units;

    /// <summary>Whether grid value <paramref name="g"/> is in the set.</summary>
    public bool Accepts(int g) => WholeSpace || Below![g] + Tied![g] + 1 > Level;

    /// <summary>The set's maximal runs of accepted grid values; empty for the whole space and for an empty set.</summary>
    public IReadOnlyList<GridInterval> Intervals
    {
        get
        {
            var runs = new List<GridInterval>();
            if (WholeSpace) return runs;
            for (int g = 0; g < Grid.Length; g++)
            {
                if (!Accepts(g)) continue;
                int first = g;
                while (g + 1 < Grid.Length && Accepts(g + 1)) g++;
                runs.Add(new GridInterval(Grid[first], Grid[g]));
            }
            return runs;
        }
    }

    /// <summary>The grid's largest gap: each boundary of the set is known to within it.</summary>
    public double Resolution
    {
        get
        {
            double gap = 0.0;
            for (int g = 1; g < Grid.Length; g++) gap = Math.Max(gap, Grid[g] - Grid[g - 1]);
            return gap;
        }
    }

    /// <summary>Whether the lowest grid value is accepted, so the set may extend below the grid.</summary>
    public bool OpenBelow => !WholeSpace && Accepts(0);

    /// <summary>Whether the highest grid value is accepted, so the set may extend above the grid.</summary>
    public bool OpenAbove => !WholeSpace && Accepts(Grid.Length - 1);

    /// <summary>
    /// The set's length read on the grid: each accepted value stands for its cell, the half-gaps to its neighbours, the
    /// grid's end values the half-gap on their inner side twice; within <see cref="Resolution"/> of the true length per
    /// boundary. +∞ for the whole space.
    /// </summary>
    public double Width
    {
        get
        {
            if (WholeSpace) return double.PositiveInfinity;
            int n = Grid.Length;
            if (n == 1) return 0.0;
            double width = 0.0;
            for (int g = 0; g < n; g++)
            {
                if (!Accepts(g)) continue;
                double left = g > 0 ? Grid[g] - Grid[g - 1] : Grid[1] - Grid[0];
                double right = g < n - 1 ? Grid[g + 1] - Grid[g] : Grid[n - 1] - Grid[n - 2];
                width += 0.5 * (left + right);
            }
            return width;
        }
    }
}

/// <summary>The full step's sets, one per test case in the candidates' order, at miscoverage <paramref name="Alpha"/>.</summary>
public sealed record ConformalPrediction(IReadOnlyList<ConformalSet> Sets, double Alpha);

/// <summary>One candidate's conformity scores (<see cref="ConformalBayes.Scores"/>): the log scores of its stratum's training
/// units, in <paramref name="Units"/>' order, the candidate's own, and the add-one-in weights' Kish ESS. A candidate outside
/// every draw's support scores −∞, with ESS 0, and the units' scores are NaN: its weights do not exist.</summary>
public sealed record ConformalScores(int[] Units, double[] Scores, double Candidate, double EffectiveSampleSize);

/// <summary>
/// Conformal Bayes, the full step (Fong &amp; Holmes 2021, <i>Conformal Bayesian Computation</i>; E16): prediction sets whose
/// coverage holds in finite samples under exchangeability, whatever the model, from retained draws with no refit. It reads a
/// training <see cref="PointwiseLogDensity"/> and a <see cref="CandidateLogDensity"/> over the same draws and knows no engine:
/// it accumulates the add-one-in weights over the draws and reduces each unit's density to a conformity score, ranked within
/// the candidate's stratum.
/// <list type="bullet">
/// <item><b>The weights.</b> For a candidate y at case t, w̃_s ∝ W_s·p(y | x_t, θ_s), the draw's weight
/// (<see cref="PointwiseLogDensity.Weights"/>, equal for chains) times the candidate's density: the draws reweighted to the
/// posterior given all n + 1 units.</item>
/// <item><b>The scores.</b> Each unit's conformity is Σ_s w̃_s p(z_i | θ_s), the self-normalized estimate of its in-sample
/// predictive p(z_i | all n + 1 units), the candidate's own included: the retained parameters' content, Fong and Holmes's.
/// For G candidates and n units they are one (n × S)·(S × G) product (<see cref="DenseProduct"/>), each unit's likelihoods
/// scaled by their largest and each entry a sum of fused multiply-adds in draw order. The candidate's own score is the same
/// sum, so a candidate whose column equals a unit's scores exactly as that unit does, and ties are exact. An entry below
/// 10⁻²⁰⁰ of its unit's largest likelihood is summed again in logs.</item>
/// <item><b>The rank</b> (CB-1, <c>coverage_orderings</c>, <c>coverage_stratum</c>): within the candidate's stratum, ties
/// counted for the candidate, π(y) = (#{units scoring ≤ the candidate} + 1)/m, accepted when π(y) &gt; k/m with k = ⌊α·m⌋. It
/// covers at least 1 − k/m whatever the ties, exactly 1 − k/m with distinct scores; the strict rule is never used.</item>
/// <item><b>Refused</b> (CB-4): leave-unit-out terms (<see cref="PointwiseLogDensity.LeavesUnitsOut"/>, the collapsed
/// predictive), whose rest excludes the candidate; a latent of one observation's term conditioned on or mixed, which reads
/// the data units at their retained latents and the candidate otherwise and undercovers (<c>mixed_counterexample</c>); a
/// candidate whose latents are treated otherwise than the training units'; draws that do not align.</item>
/// </list>
/// The guarantee is exact for exact posterior weights; from Monte Carlo draws it holds up to the estimates' error (CB-3).
/// A candidate whose density is zero under every draw scores −∞ and takes π = 1/m. Thread-safe: it reads its inputs and
/// writes only its own buffers.
/// </summary>
public static class ConformalBayes
{
    // Candidates scored per product: a multiple of the product's panel, so only a case's last batch has a partial panel.
    private const int Batch = 256;

    // A product entry below this share of its unit's largest likelihood is summed again in logs: above it, the terms that
    // underflow in the product weigh less than S·2⁻¹⁰²²/10⁻²⁰⁰ of the entry.
    private const double Floor = 1e-200;

    /// <summary>
    /// The full step at every test case of <paramref name="candidates"/>: π(y) on each case's grid within its stratum, the set,
    /// and each candidate's ESS (see <see cref="ConformalBayes"/>).
    /// </summary>
    public static ConformalPrediction FullStep(PointwiseLogDensity training, CandidateLogDensity candidates, ConformalOptions? options = null)
    {
        options ??= new ConformalOptions();
        double alpha = options.Alpha;
        if (!(alpha > 0.0 && alpha < 1.0)) throw new ArgumentOutOfRangeException(nameof(options), "α must lie in (0, 1).");
        Check(training, candidates);
        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(training);
        Dictionary<int, List<int>> units = UnitsByStratum(training.Units);

        var sets = new ConformalSet[candidates.Cases];
        var casesOf = new SortedDictionary<int, List<int>>();
        for (int t = 0; t < candidates.Cases; t++)
        {
            int h = candidates.StratumOf[t];
            if (!casesOf.TryGetValue(h, out List<int>? list)) casesOf.Add(h, list = []);
            list.Add(t);
        }
        foreach ((int h, List<int> cases) in casesOf)
        {
            int[] members = units.TryGetValue(h, out List<int>? found) ? [.. found] : [];
            int m = members.Length + 1, k = (int)Math.Floor(alpha * m);
            if (k == 0)
            {
                foreach (int t in cases) sets[t] = new ConformalSet(t, h, m, k, candidates.Grids[t], null, null, null, null, null);
                continue;
            }
            var scorer = new StratumScorer(training.Values, rows, logWeights, members);
            foreach (int t in cases)
            {
                int g = candidates.Grids[t].Length, offset = candidates.Offset(t);
                var below = new int[g];
                var tied = new int[g];
                var own = new double[g];
                var ess = new double[g];
                for (int b0 = 0; b0 < g; b0 += Batch)
                    scorer.Rank(candidates.Values, offset + b0, Math.Min(Batch, g - b0), below.AsSpan(b0), tied.AsSpan(b0), own.AsSpan(b0), ess.AsSpan(b0));
                var p = new double[g];
                for (int j = 0; j < g; j++) p[j] = (below[j] + tied[j] + 1) / (double)m;
                sets[t] = new ConformalSet(t, h, m, k, candidates.Grids[t], p, below, tied, own, ess);
            }
        }
        return new ConformalPrediction(sets, alpha);
    }

    /// <summary>
    /// The conformity scores at candidate column <paramref name="candidate"/>: every training unit of its case's stratum, the
    /// candidate's own, and the add-one-in weights' ESS, computed as <see cref="FullStep"/> computes them, bit for bit.
    /// </summary>
    public static ConformalScores Scores(PointwiseLogDensity training, CandidateLogDensity candidates, int candidate)
    {
        Check(training, candidates);
        if ((uint)candidate >= (uint)candidates.Values.GetLength(1)) throw new ArgumentOutOfRangeException(nameof(candidate));
        int t = 0;
        while (t + 1 < candidates.Cases && candidates.Offset(t + 1) <= candidate) t++;
        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(training);
        int[] members = UnitsByStratum(training.Units).TryGetValue(candidates.StratumOf[t], out List<int>? found) ? [.. found] : [];
        var scores = new double[members.Length];
        double own = 0.0, ess = 0.0;
        new StratumScorer(training.Values, rows, logWeights, members).Score(candidates.Values, candidate, scores, ref own, ref ess);
        return new ConformalScores(members, scores, own, ess);
    }

    private static Dictionary<int, List<int>> UnitsByStratum(UnitPartition partition)
    {
        var units = new Dictionary<int, List<int>>();
        for (int u = 0; u < partition.Units; u++)
        {
            int h = partition.StratumOf[u];
            if (!units.TryGetValue(h, out List<int>? list)) units.Add(h, list = []);
            list.Add(u);
        }
        return units;
    }

    // CB-4's refusals and the draws' alignment.
    private static void Check(PointwiseLogDensity training, CandidateLogDensity candidates)
    {
        ArgumentNullException.ThrowIfNull(training);
        ArgumentNullException.ThrowIfNull(candidates);
        if (training.LeavesUnitsOut)
            throw new ArgumentException(
                "The training terms are leave-unit-out predictives (the collapsed content): each leaves its unit out of a rest that "
                + "excludes the candidate, which is not the full step's in-sample score. Take the retained parameters' content.", nameof(training));
        foreach (LatentReport latent in training.Latents)
            if (latent.Scope == LatentScope.Observation && latent.Treatment is LatentTreatment.Conditioned or LatentTreatment.Mixed)
                throw new ArgumentException(
                    $"The latent '{latent.Name}' of each observation's term is {(latent.Treatment == LatentTreatment.Mixed ? "mixed" : "conditioned on")}: "
                    + "the scores would read the data units at their retained latents and the candidate by another density, which undercovers "
                    + "(CB-4). Take the density with it integrated out.", nameof(training));
        if (!SameLatents(training.Latents, candidates.Latents))
            throw new ArgumentException("The candidates' latents are not treated as the training units' are (CB-4): every unit, the candidate included, is read by one density.", nameof(candidates));
        if (candidates.Draws != training.Draws)
            throw new ArgumentException($"The candidates hold {candidates.Draws} draws, the training density {training.Draws}: the rows must be the same draws.", nameof(candidates));
        if (!SameSequence(training.ChainLengths, candidates.ChainLengths))
            throw new ArgumentException("The candidates' chain lengths are not the training density's.", nameof(candidates));
        if (!SameSequence(training.Weights, candidates.Weights))
            throw new ArgumentException("The candidates' draw weights are not the training density's.", nameof(candidates));
    }

    private static bool SameLatents(IReadOnlyList<LatentReport> a, IReadOnlyList<LatentReport> b)
    {
        if (a.Count != b.Count) return false;
        var counts = new Dictionary<LatentReport, int>();
        foreach (LatentReport l in a) counts[l] = counts.GetValueOrDefault(l) + 1;
        foreach (LatentReport l in b)
        {
            if (!counts.TryGetValue(l, out int c) || c == 0) return false;
            counts[l] = c - 1;
        }
        return true;
    }

    private static bool SameSequence<T>(IReadOnlyList<T>? a, IReadOnlyList<T>? b) where T : IEquatable<T>
    {
        if (a is null || b is null) return a is null && b is null;
        if (a.Count != b.Count) return false;
        for (int i = 0; i < a.Count; i++)
            if (!a[i].Equals(b[i])) return false;
        return true;
    }

    /// <summary>
    /// One stratum's training units against batches of candidates. It holds each unit's likelihoods over the draws of positive
    /// weight scaled by their largest, P_is = p(z_i | θ_s)/max_s, row by row, the product's left operand, and for each batch the
    /// candidates' normalized add-one-in weights draw by draw, its right one. Not thread-safe: one per stratum and call.
    /// </summary>
    internal sealed class StratumScorer
    {
        private readonly double[,] _ll;
        private readonly int[] _rows, _members;
        private readonly double[]? _logWeights;
        private readonly double[] _left, _top;   // P, units × draws; each unit's largest log term
        private readonly int _s;
        private double[] _weights = [], _likelihoods = [], _logNorm = [], _candidateTop = [], _product = [];
        private PackedColumns? _right;

        public StratumScorer(double[,] ll, int[] rows, double[]? logWeights, int[] members)
        {
            _ll = ll;
            _rows = rows;
            _logWeights = logWeights;
            _members = members;
            _s = rows.Length;
            int n = members.Length, s = _s;
            _top = new double[n];
            Array.Fill(_top, double.NegativeInfinity);
            for (int d = 0; d < s; d++)
            {
                int row = rows[d];
                for (int i = 0; i < n; i++)
                {
                    double v = ll[row, members[i]];
                    if (double.IsNaN(v) || double.IsPositiveInfinity(v)) throw new ArgumentException($"Unit {members[i]}'s log term is NaN or +∞ under draw {row}.", "training");
                    if (v > _top[i]) _top[i] = v;
                }
            }
            for (int i = 0; i < n; i++)
                if (double.IsNegativeInfinity(_top[i])) throw new ArgumentException($"Unit {members[i]} has zero density under every draw of positive weight.", "training");
            // P in tiles of 32 draws by 32 units, so the strided reads of the draws × units matrix stay in cache.
            _left = new double[(long)n * s];
            const int tile = 32;
            for (int d0 = 0; d0 < s; d0 += tile)
                for (int i0 = 0; i0 < n; i0 += tile)
                    for (int i = i0; i < Math.Min(n, i0 + tile); i++)
                    {
                        int u = members[i];
                        double top = _top[i];
                        for (int d = d0; d < Math.Min(s, d0 + tile); d++) _left[(long)i * s + d] = Math.Exp(ll[rows[d], u] - top);
                    }
        }

        /// <summary>The ranks of <paramref name="count"/> candidates from column <paramref name="first"/>: below and tied counts,
        /// own log scores and ESS. Its three steps are apart for the cost row: <see cref="Weigh"/>, <see cref="Multiply"/>,
        /// <see cref="Count"/>.</summary>
        public void Rank(double[,] candidates, int first, int count, Span<int> below, Span<int> tied, Span<double> own, Span<double> ess)
        {
            Weigh(candidates, first, count, own, ess);
            Multiply(count);
            Count(candidates, first, count, below, tied, own);
        }

        /// <summary>The training units' scores against the batch the last <see cref="Weigh"/> prepared: P·W̃, units × candidates.</summary>
        public void Multiply(int count)
        {
            int n = _members.Length;
            if (n == 0) return;
            if (_product.Length < n * count) _product = new double[n * count];
            DenseProduct.Multiply(_left, n, _right!, _product);
        }

        /// <summary>Each candidate's below and tied counts from the product, against its own score.</summary>
        public void Count(double[,] candidates, int first, int count, Span<int> below, Span<int> tied, ReadOnlySpan<double> own)
        {
            below[..count].Clear();
            tied[..count].Clear();
            for (int i = 0; i < _members.Length; i++)
            {
                ReadOnlySpan<double> row = _product.AsSpan(i * count, count);
                double top = _top[i];
                for (int b = 0; b < count; b++)
                {
                    double candidate = own[b];
                    if (double.IsNegativeInfinity(candidate)) continue;   // outside every draw's support: π = 1/m
                    double score = row[b] >= Floor ? top + Math.Log(row[b]) : LogScore(candidates, first + b, b, _members[i]);
                    if (score < candidate) below[b]++;
                    else if (score == candidate) tied[b]++;
                }
            }
        }

        /// <summary>One candidate's scores, unit by unit, by the same arithmetic as <see cref="Rank"/>.</summary>
        public void Score(double[,] candidates, int column, Span<double> scores, ref double own, ref double ess)
        {
            Span<double> o = stackalloc double[1], e = stackalloc double[1];
            Weigh(candidates, column, 1, o, e);
            own = o[0];
            ess = e[0];
            int n = _members.Length;
            Multiply(1);
            for (int i = 0; i < n; i++)
                scores[i] = double.IsNegativeInfinity(own) ? double.NaN
                          : _product[i] >= Floor ? _top[i] + Math.Log(_product[i]) : LogScore(candidates, column, 0, _members[i]);
        }

        /// <summary>The batch's normalized add-one-in weights, draws × candidates, packed as the product's right operand; each
        /// candidate's scaled likelihoods, its log normalizer, its own log score and its ESS.</summary>
        public void Weigh(double[,] candidates, int first, int count, Span<double> own, Span<double> ess)
        {
            int s = _s, size = s * count;
            if (_weights.Length < size)
            {
                _weights = new double[size];
                _likelihoods = new double[size];
            }
            if (_logNorm.Length < count)
            {
                _logNorm = new double[count];
                _candidateTop = new double[count];
            }
            Span<double> top = _candidateTop.AsSpan(0, count), norm = _logNorm.AsSpan(0, count);
            Span<double> w = _weights.AsSpan(0, size), p = _likelihoods.AsSpan(0, size);
            top.Fill(double.NegativeInfinity);
            for (int d = 0; d < s; d++)
            {
                int row = _rows[d];
                for (int b = 0; b < count; b++) top[b] = Math.Max(top[b], candidates[row, first + b]);
            }
            // p = exp(ln p − its largest): the candidate's own likelihoods, scaled as a unit's are, so equal columns agree.
            for (int d = 0; d < s; d++)
            {
                int row = _rows[d];
                Span<double> pd = p.Slice(d * count, count);
                for (int b = 0; b < count; b++) pd[b] = double.IsNegativeInfinity(top[b]) ? 0.0 : Math.Exp(candidates[row, first + b] - top[b]);
            }
            if (_logWeights is null)
            {
                // Equal weights: w̃ = p/Σp, and ln Σ p(y | θ_s) = top + ln Σp.
                p.CopyTo(w);
            }
            else
            {
                // w ∝ W_s·p(y | θ_s), taken in logs about its own largest.
                norm.Fill(double.NegativeInfinity);
                for (int d = 0; d < s; d++)
                {
                    int row = _rows[d];
                    for (int b = 0; b < count; b++) norm[b] = Math.Max(norm[b], _logWeights[d] + candidates[row, first + b]);
                }
                for (int d = 0; d < s; d++)
                {
                    int row = _rows[d];
                    Span<double> wd = w.Slice(d * count, count);
                    for (int b = 0; b < count; b++) wd[b] = double.IsNegativeInfinity(norm[b]) ? 0.0 : Math.Exp(_logWeights[d] + candidates[row, first + b] - norm[b]);
                }
            }
            var totals = new double[count];
            for (int d = 0; d < s; d++)
            {
                Span<double> wd = w.Slice(d * count, count);
                for (int b = 0; b < count; b++) totals[b] += wd[b];
            }
            for (int b = 0; b < count; b++)
            {
                double shift = _logWeights is null ? top[b] : norm[b];
                norm[b] = totals[b] > 0.0 ? shift + Math.Log(totals[b]) : double.NegativeInfinity;
                totals[b] = totals[b] > 0.0 ? 1.0 / totals[b] : 0.0;
            }
            ess[..count].Clear();
            for (int d = 0; d < s; d++)
            {
                Span<double> wd = w.Slice(d * count, count);
                for (int b = 0; b < count; b++)
                {
                    wd[b] *= totals[b];
                    ess[b] += wd[b] * wd[b];
                }
            }
            for (int b = 0; b < count; b++) ess[b] = ess[b] > 0.0 ? 1.0 / ess[b] : 0.0;
            OwnSums(w, p, s, count, own);
            for (int b = 0; b < count; b++)
                own[b] = double.IsNegativeInfinity(top[b]) || double.IsNegativeInfinity(norm[b]) ? double.NegativeInfinity
                       : own[b] >= Floor ? top[b] + Math.Log(own[b]) : LogScore(candidates, first + b, b, -1);
            _right ??= new PackedColumns(s, count);
            if (_right.Columns != count) _right = new PackedColumns(s, count);
            _right.Fill(w);
        }

        // ln Σ_s w̃_s p(z | θ_s) summed in logs, for unit `unit`, or the candidate's own column when unit is −1.
        private double LogScore(double[,] candidates, int column, int b, int unit)
        {
            double norm = _logNorm[b];
            double max = double.NegativeInfinity;
            for (int d = 0; d < _s; d++) max = Math.Max(max, Term(candidates, column, norm, unit, d));
            if (double.IsNegativeInfinity(max)) return max;
            double sum = 0.0;
            for (int d = 0; d < _s; d++) sum += Math.Exp(Term(candidates, column, norm, unit, d) - max);
            return max + Math.Log(sum);
        }

        private double Term(double[,] candidates, int column, double norm, int unit, int d)
        {
            int row = _rows[d];
            double logWeight = (_logWeights is null ? 0.0 : _logWeights[d]) + candidates[row, column] - norm;
            return logWeight + (unit < 0 ? candidates[row, column] : _ll[row, unit]);
        }

        // own[b] = Σ_d w[d, b]·p[d, b] by fused multiply-adds in draw order, four candidates to a vector: each lane is the
        // product's entry for a unit whose column equals the candidate's, bit for bit.
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private static void OwnSums(ReadOnlySpan<double> w, ReadOnlySpan<double> p, int s, int count, Span<double> own)
        {
            ref double wr = ref MemoryMarshal.GetReference(w);
            ref double pr = ref MemoryMarshal.GetReference(p);
            int b = 0;
            for (; b + 4 <= count; b += 4)
            {
                Vector256<double> acc = Vector256<double>.Zero;
                for (int d = 0; d < s; d++)
                {
                    nuint at = (nuint)(d * count + b);
                    acc = Vector256.FusedMultiplyAdd(Vector256.LoadUnsafe(ref wr, at), Vector256.LoadUnsafe(ref pr, at), acc);
                }
                acc.CopyTo(own.Slice(b, 4));
            }
            for (; b < count; b++)
            {
                double acc = 0.0;
                for (int d = 0; d < s; d++) acc = Math.FusedMultiplyAdd(w[d * count + b], p[d * count + b], acc);
                own[b] = acc;
            }
        }
    }
}
