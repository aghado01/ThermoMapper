using System;
using System.Collections.Generic;
using Maths.Rng;

namespace Maths.Regression.Conformal;

/// <summary>The full form's conformity score, higher conforming better, recomputed for every unit at every arrival under the draws
/// reweighted by the stream so far (<see cref="FullFormPValues"/>).</summary>
public enum FullFormScore
{
    /// <summary>The in-sample predictive density ln p̂ₙ(yᵢ | xᵢ) = ln Σ_s w_s⁽ⁿ⁾ p(yᵢ | xᵢ, θ_s): the full step's own score.</summary>
    Density,

    /// <summary>The in-sample PIT's centrality min(F̂ₙ, 1 − F̂ₙ) at yᵢ, F̂ₙ(y | x) = Σ_s w_s⁽ⁿ⁾ F(y | x, θ_s): the full step's CDF score,
    /// read by both tails. A unit far in either tail conforms least. It needs each draw's log CDF.</summary>
    Centrality,
}

/// <summary>
/// The smoothed conformal p-values of a stream in the full form (E35; Vovk, Nouretdinov and Gammerman 2003, Theorem 1). At each arrival
/// every unit's score is recomputed under the posterior given the stream so far, and the newest score is ranked among them all:
/// pₙ = (Lₙ + θₙTₙ)/n, as <see cref="SmoothedPValues"/> ranks.
/// <para><b>The posterior given the stream.</b> The fit's draws θ_s come from a training set proper outside the stream. Each arrival
/// multiplies every draw's weight by the newest unit's term, so after n arrivals w_s⁽ⁿ⁾ ∝ W̄_s Π_{i≤n} p(yᵢ | xᵢ, θ_s), the draws'
/// own weights W̄ times the stream's likelihood: the add-one-in weights over the whole stream. Each unit is then scored in sample under
/// w⁽ⁿ⁾ (<see cref="FullFormScore"/>).</para>
/// <para><b>The guarantee</b> (<c>seqRank_uniform</c>, <c>seqRank_law</c>): if the stream's units are exchangeable, p₁, p₂, … are i.i.d.
/// uniform. M-1 needs the scores at each time to be a symmetric function of the units so far (<c>Equivariant</c>), and these are. Given
/// the draws, every score at time n reads the units only through the product over all of them, whatever their order. So the guarantee
/// is exact whatever the number of draws, and it needs no Monte Carlo margin, unlike the full step's sets (CB-3). It does need the draws
/// fixed before the stream, independent of it. A refit to the stream would break it.</para>
/// <para><b>What the draws decide.</b> How sharp the scores are. As the stream grows, the weights concentrate on fewer draws
/// (<see cref="EffectiveSampleSize"/>). At one draw the score is that draw's own density, still valid but blind to the rest of the
/// posterior.</para>
/// <para><b>Cost.</b> Each arrival rescores every unit, O(n·S) for S draws, against the inductive form's O(log n). Memory is two numbers
/// per unit and draw.</para>
/// Not thread-safe: one per stream.
/// </summary>
public sealed class FullFormPValues
{
    private const double FastFloor = 1e-250;   // below this, a unit's weighted sum is recomputed in logs
    private readonly int _draws;
    private readonly IReadOnlyList<double>? _drawWeights;
    private readonly int[] _rows;               // the draws of positive weight
    private readonly double[] _logWeight;       // each such draw's ln W̄_s + Σ_{i≤n} ln p(yᵢ | xᵢ, θ_s)
    private readonly double[] _weight;          // exp(ln w − its largest), at the latest arrival
    private double _scaledTotal;                // their sum
    private readonly List<double[]> _first = [];   // per unit: exp(ln p − its largest) per draw (density), or F per draw (centrality)
    private readonly List<double[]> _second = [];  // per unit: ln p per draw (density), or 1 − F per draw (centrality)
    private readonly List<double> _peak = [];      // per unit: its largest ln p over the draws (density)
    private readonly Xoshiro256PlusPlus _theta;
    private double[] _scores = [];

    /// <param name="draws">The fit's draws S, the rows of every unit's terms.</param>
    /// <param name="seed">The tie-breaking draws' seed; the draws are independent of the scores by construction.</param>
    /// <param name="drawWeights">The draws' own weights W̄ in row order, as <see cref="PointwiseLogDensity.Weights"/>; null for equal weights.</param>
    /// <param name="score">The score recomputed at each arrival.</param>
    public FullFormPValues(int draws, int seed, IReadOnlyList<double>? drawWeights = null, FullFormScore score = FullFormScore.Density)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(draws, 1);
        if (drawWeights is not null)
        {
            if (drawWeights.Count != draws) throw new ArgumentException($"One weight per draw is required: {drawWeights.Count} for {draws} draws.", nameof(drawWeights));
            foreach (double w in drawWeights)
                if (!(w >= 0.0 && double.IsFinite(w))) throw new ArgumentException("Each draw's weight is finite and nonnegative.", nameof(drawWeights));
        }
        if (score is not (FullFormScore.Density or FullFormScore.Centrality)) throw new ArgumentOutOfRangeException(nameof(score));
        (_rows, double[]? logWeights) = PredictiveChecks.Rows(drawWeights, draws);
        if (_rows.Length == 0) throw new ArgumentException("No draw has positive weight.", nameof(drawWeights));
        _draws = draws;
        _drawWeights = drawWeights is null ? null : [.. drawWeights];
        _logWeight = logWeights ?? new double[_rows.Length];
        _weight = new double[_rows.Length];
        Score = score;
        _theta = new Xoshiro256PlusPlus(seed);
        Reweigh();
    }

    /// <summary>The score recomputed at each arrival.</summary>
    public FullFormScore Score { get; }

    /// <summary>The fit's draws S.</summary>
    public int Draws => _draws;

    /// <summary>The arrivals so far.</summary>
    public int Count => _first.Count;

    /// <summary>The effective sample size (Σw)²/Σw² of the draws reweighted by the stream so far; before any arrival, of the draws' own
    /// weights. It falls as the stream grows; at 1 a single draw scores every unit.</summary>
    public double EffectiveSampleSize { get; private set; }

    /// <summary>Every unit's score at the latest arrival, in arrival order.</summary>
    public IReadOnlyList<double> Scores => _scores;

    /// <summary>
    /// Takes the newest unit's terms, reweighs the draws by them, rescores every unit and ranks the newest among them all.
    /// </summary>
    /// <param name="logDensity">Each draw's ln p(yₙ | xₙ, θ_s) at the unit, in row order: finite or −∞ (outside the draw's support), and finite
    /// for some draw the stream so far leaves weight on.</param>
    /// <param name="logCdf">Each draw's ln F(yₙ | xₙ, θ_s), both tails kept as <see cref="PointwiseLogDensity.LogCdf"/> keeps them; read for
    /// <see cref="FullFormScore.Centrality"/>, ignored otherwise.</param>
    public SmoothedPValue Add(ReadOnlySpan<double> logDensity, ReadOnlySpan<double> logCdf = default)
    {
        if (logDensity.Length != _draws) throw new ArgumentException($"One term per draw is required: {logDensity.Length} for {_draws} draws.", nameof(logDensity));
        bool weighed = false;
        for (int k = 0; k < _rows.Length; k++)
        {
            double v = logDensity[_rows[k]];
            if (double.IsNaN(v) || double.IsPositiveInfinity(v)) throw new ArgumentException("A draw's log density is NaN or +∞.", nameof(logDensity));
            if (_logWeight[k] > double.NegativeInfinity && v > double.NegativeInfinity) weighed = true;
        }
        if (!weighed)
            throw new ArgumentException("The unit has zero density under every draw the stream so far leaves weight on, so no posterior holds it with the stream.", nameof(logDensity));
        double[] first = new double[_rows.Length], second = new double[_rows.Length];
        if (Score == FullFormScore.Centrality)
        {
            if (logCdf.Length != _draws) throw new ArgumentException($"The centrality reads each draw's log CDF: {logCdf.Length} for {_draws} draws.", nameof(logCdf));
            for (int k = 0; k < _rows.Length; k++)
            {
                double c = logCdf[_rows[k]];
                if (!(c <= 0.0)) throw new ArgumentException("A draw's log CDF is NaN or above 0.", nameof(logCdf));
                first[k] = Math.Exp(c);
                second[k] = Math.Exp(GPrior.LogFunctions.Log1mExp(c));
            }
        }
        else
        {
            double peak = double.NegativeInfinity;
            for (int k = 0; k < _rows.Length; k++) peak = Math.Max(peak, logDensity[_rows[k]]);
            for (int k = 0; k < _rows.Length; k++)
            {
                second[k] = logDensity[_rows[k]];
                first[k] = Math.Exp(second[k] - peak);
            }
            _peak.Add(peak);
        }
        for (int k = 0; k < _rows.Length; k++) _logWeight[k] += logDensity[_rows[k]];
        _first.Add(first);
        _second.Add(second);
        double logTotal = Reweigh();

        int n = _first.Count;
        var scores = new double[n];
        for (int i = 0; i < n; i++) scores[i] = ScoreOf(i, logTotal);
        _scores = scores;
        double newest = scores[n - 1];
        int below = 0, tied = 0;
        foreach (double s in scores)
        {
            if (s < newest) below++;
            else if (s == newest) tied++;
        }
        double theta = 1.0 - _theta.NextDouble();
        return new SmoothedPValue(n, newest, below, tied, theta, (below + theta * tied) / n);
    }

    /// <summary>Each unit of <paramref name="units"/> in turn (<see cref="Add"/>), through the candidate seam: each case's grid holds its one
    /// observed response, the seam's draws and weights are this stream's, and it serves the full form, since reweighing the draws by its
    /// terms reads each term as that unit's likelihood.</summary>
    public SmoothedPValue[] AddRange(CandidateLogDensity units)
    {
        var steps = new List<SmoothedPValue>();
        foreach ((double[] density, double[]? cdf) in Units(units)) steps.Add(Add(density, cdf));
        return [.. steps];
    }

    // The seam's units in turn, each as its draws' log densities and, for the centrality, log CDFs: checked whole before the first.
    internal IEnumerable<(double[] Density, double[]? Cdf)> Units(CandidateLogDensity units)
    {
        ArgumentNullException.ThrowIfNull(units);
        CheckSeam(units);
        return Columns(units);
    }

    private IEnumerable<(double[] Density, double[]? Cdf)> Columns(CandidateLogDensity units)
    {
        for (int t = 0; t < units.Cases; t++)
        {
            var density = new double[_draws];
            for (int s = 0; s < _draws; s++) density[s] = units.Values[s, t];
            double[]? cdf = null;
            if (Score == FullFormScore.Centrality)
            {
                cdf = new double[_draws];
                for (int s = 0; s < _draws; s++) cdf[s] = units.LogCdf![s, t];
            }
            yield return (density, cdf);
        }
    }

    private void CheckSeam(CandidateLogDensity units)
    {
        units.CheckServes(ConformalForms.Full, nameof(units));
        if (units.Draws != _draws) throw new ArgumentException($"The seam holds {units.Draws} draws; this stream's fit holds {_draws}.", nameof(units));
        if (!SameWeights(units.Weights)) throw new ArgumentException("The seam's draw weights are not this stream's.", nameof(units));
        for (int t = 0; t < units.Cases; t++)
            if (units.Grids[t].Length != 1)
                throw new ArgumentException($"Unit {t}'s grid holds {units.Grids[t].Length} values: each unit is read at its observed response.", nameof(units));
        if (Score == FullFormScore.Centrality && units.LogCdf is null)
            throw new ArgumentException("The units carry no log CDF; the centrality reads the fit's predictive CDF.", nameof(units));
    }

    private bool SameWeights(IReadOnlyList<double>? weights)
    {
        if (weights is null || _drawWeights is null) return weights is null && _drawWeights is null;
        if (weights.Count != _drawWeights.Count) return false;
        for (int s = 0; s < weights.Count; s++)
            if (weights[s] != _drawWeights[s]) return false;
        return true;
    }

    // Sets the scaled weights, their sum and the effective sample size from the log weights; returns ln Σ_s exp(ln w_s).
    private double Reweigh()
    {
        double largest = double.NegativeInfinity;
        foreach (double lw in _logWeight) largest = Math.Max(largest, lw);
        double total = 0.0, squares = 0.0;
        for (int k = 0; k < _weight.Length; k++)
        {
            double w = Math.Exp(_logWeight[k] - largest);
            _weight[k] = w;
            total += w;
            squares += w * w;
        }
        _scaledTotal = total;
        EffectiveSampleSize = total * total / squares;
        return largest + Math.Log(total);
    }

    // Unit i's score under the current weights; logTotal = ln Σ_s exp(ln w_s).
    private double ScoreOf(int i, double logTotal)
    {
        double[] first = _first[i];
        double sum = 0.0;
        if (Score == FullFormScore.Centrality)
        {
            double upper = 0.0;
            double[] second = _second[i];
            for (int k = 0; k < first.Length; k++)
            {
                sum += _weight[k] * first[k];
                upper += _weight[k] * second[k];
            }
            return Math.Min(sum, upper) / _scaledTotal;
        }
        for (int k = 0; k < first.Length; k++) sum += _weight[k] * first[k];
        if (sum >= FastFloor) return _peak[i] + Math.Log(sum / _scaledTotal);
        // The weighted terms underflow in the scaled sum: ln Σ_s exp(ln w_s + ln p(yᵢ | θ_s)) − ln Σ_s w_s, in logs.
        double[] log = _second[i];
        double top = double.NegativeInfinity;
        for (int k = 0; k < log.Length; k++) top = Math.Max(top, _logWeight[k] + log[k]);
        double exact = 0.0;
        for (int k = 0; k < log.Length; k++) exact += Math.Exp(_logWeight[k] + log[k] - top);
        return top + Math.Log(exact) - logTotal;
    }
}
