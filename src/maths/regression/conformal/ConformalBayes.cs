using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using Maths.Distributions;
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

    /// <summary>
    /// ν, in (0, 1): the chance, under the Monte Carlo error's normal approximation, that some unit of a stratum's m strays beyond
    /// its bound. Each score's bound ε is z = Φ⁻¹(1 − ν/(2m)) of its Monte Carlo standard errors, and ν adds to the miscoverage
    /// the bounds give (CB-3), so it is read in <see cref="ConformalSet.CoverageUnderError"/>. Indicative: the bound is not a sure
    /// one. 0.01 by default.
    /// </summary>
    public double ErrorLevel { get; init; } = 0.01;

    /// <summary>Whether the set is the inflated rule's (Ndiaye's stable set, <c>coverage_inflated</c>), which keeps 1 − k/m under
    /// the bounds at the cost of width, rather than the plain rule's. Off by default; the margin is reported either way.</summary>
    public bool Inflate { get; init; }

    /// <summary>
    /// The fallback rule, opt-in: when set, <see cref="ConformalSet.FallBack"/> fires when the inflated set exceeds the plain one
    /// by more than this share of the plain set (<see cref="ConformalSet.Gap"/>), or when a candidate the inflated rule accepts
    /// has add-one-in weights whose Pareto k̂ exceeds PSIS's threshold. Null, the default, fires nothing: the layer reports its
    /// diagnostics (<see cref="ConformalSet.Gap"/>, <see cref="ConformalSet.SplitHalfMass"/>) and switches to nothing on its own,
    /// the user's decision of 2026-11-07, since no fixture with chains like BASS's yet shows a shortfall to set a threshold by.
    /// </summary>
    public double? FallbackShare { get; init; }

    /// <summary>The bounds the margin's band reads (<see cref="Conformal.MarginBounds"/>): the uniform one by default.</summary>
    public MarginBounds MarginBounds { get; init; } = MarginBounds.Uniform;
}

/// <summary>
/// The bounds the margin reads. <c>coverage_close</c> needs bounds carried by the units, holding in every ordering; both readings
/// take them from the orderings the layer holds, so both are indicative.
/// </summary>
public enum MarginBounds
{
    /// <summary>One bound for every unit, the largest of the stratum's units' and the candidate's under this candidate's
    /// weights: the safe reading, which a few units with wide errors widen for all.</summary>
    Uniform,

    /// <summary>
    /// Each unit its own: the larger of its bound under this candidate's weights and its bound as the candidate, under weights
    /// of its own density, since every unit is the candidate in some ordering (the Lean chip's per-unit finding). It costs each
    /// unit's own weighting and, for chains, its relative efficiency, once per stratum.
    /// </summary>
    Carried,
}

/// <summary>A run of accepted grid values, from its lowest to its highest; the set's boundary lies in the gaps beside it.</summary>
public readonly record struct GridInterval(double Lower, double Upper);

/// <summary>
/// CB-3's reading of one case's set (<see cref="ConformalSet.Margin"/>): each score's Monte Carlo error bound, and what the bounds
/// can cost. Scores and bounds are on the log scale. Per grid value y_g, under that candidate's add-one-in weights:
/// <list type="bullet">
/// <item><see cref="CandidateErrors"/>: the candidate's own bound; with each training unit's, the inflated rule's per-position
/// bounds (<c>coverage_inflated</c>).</item>
/// <item><see cref="UniformErrors"/>: the largest bound among the stratum's units and the candidate, the margin's bound carried
/// by every unit (the safe reading of <c>coverage_close</c>'s unit-carried bounds).</item>
/// <item><see cref="Fragile"/>: the estimated band's count (<c>card_fragile_le_estimated</c>) under the uniform bound or each
/// unit's carried one (<see cref="Bounds"/>), the units that estimates within their bounds could refuse; −1 outside every
/// draw's support.</item>
/// <item><see cref="InflatedRanks"/>: the units whose estimate, lowered by its bound, is at most the candidate's raised by its
/// own, the candidate counted; the inflated rule accepts y_g when it exceeds k.</item>
/// <item><see cref="ParetoK"/>: k̂ of the add-one-in weights, read where the inflated rule accepts and NaN elsewhere; above
/// <see cref="ParetoKThreshold"/> the bounds there cannot be trusted.</item>
/// </list>
/// </summary>
/// <param name="Multiple">z: each bound is z Monte Carlo standard errors.</param>
/// <param name="CandidateErrors">The candidate's bound at each grid value; NaN outside every draw's support.</param>
/// <param name="UniformErrors">The uniform bound at each grid value; NaN outside every draw's support.</param>
/// <param name="Fragile">The estimated band's count at each grid value.</param>
/// <param name="InflatedRanks">The inflated rule's rank at each grid value.</param>
/// <param name="ParetoK">k̂ of the add-one-in weights where the inflated rule accepts.</param>
/// <param name="ParetoKThreshold">PSIS's threshold at the draws' count.</param>
/// <param name="Bounds">The bounds <see cref="Fragile"/> was read with.</param>
public sealed record ConformalMargin(double Multiple, double[] CandidateErrors, double[] UniformErrors, int[] Fragile,
                                     int[] InflatedRanks, double[] ParetoK, double ParetoKThreshold, MarginBounds Bounds = MarginBounds.Uniform)
{
    /// <summary>m̂: the largest band over the grid's values inside some draw's support, since any of them may be the truth.</summary>
    public int Count
    {
        get
        {
            int most = 0;
            foreach (int f in Fragile) most = Math.Max(most, f);
            return most;
        }
    }

    /// <summary>The grid values the inflated rule accepts whose weights' k̂ exceeds the threshold.</summary>
    public int Untrusted(int level)
    {
        int count = 0;
        for (int g = 0; g < ParetoK.Length; g++)
            if (InflatedRanks[g] > level && !(ParetoK[g] <= ParetoKThreshold)) count++;
        return count;
    }
}

/// <summary>
/// One test case's conformal set on its grid (<see cref="ConformalBayes.FullStep"/>). The candidate at grid value y_g is ranked
/// among the <see cref="Units"/> − 1 training units of its <see cref="Stratum"/>: <see cref="Below"/> score strictly below it and
/// <see cref="Tied"/> equal to it, so π(y_g) = (Below + Tied + 1)/m, ties counted for the candidate (CB-1), and y_g is accepted
/// when π(y_g) > k/m, or under <see cref="Inflated"/> when the inflated rank exceeds k. π need not be monotone in y, so the set
/// can be a union of intervals (<see cref="Intervals"/>). A stratum too small for the level (k = 0) has the whole space as its
/// set, and nothing is computed: <see cref="WholeSpace"/>, with the per-value arrays and the margin null. The split form's density
/// score (<see cref="SplitConformal.Density"/>) reports its sets here too, ranked among the stratum's calibration units: its
/// scores are the training split's predictive, with no add-one-in weights, so no ESS, no margin and no diagnostics, and the
/// draws' error costs no coverage (CB-2).
/// </summary>
/// <param name="Case">The test case, in the candidates' order.</param>
/// <param name="Stratum">Its stratum.</param>
/// <param name="Units">m: the stratum's training units and the candidate; in the split form, its calibration units and the
/// candidate.</param>
/// <param name="Level">k = ⌊α·m⌋.</param>
/// <param name="Grid">The case's candidate values, ascending.</param>
/// <param name="PValues">π(y_g) for each grid value; null for the whole space.</param>
/// <param name="Below">The training units of the stratum scoring strictly below the candidate at y_g.</param>
/// <param name="Tied">Those scoring equal to it; the strict rule (accept when Below ≥ k) and random tie-breaking read
/// these with <see cref="Below"/>.</param>
/// <param name="CandidateScores">The candidate's own log conformity score at y_g, ln p(y_g | all m units), in the split form
/// ln p(y_g | the training split); −∞ outside every draw's support.</param>
/// <param name="EffectiveSampleSizes">Kish's ESS of the add-one-in weights at y_g, 1/Σ w̃²; 0 outside every draw's
/// support; null in the split form.</param>
/// <param name="Margin">CB-3's reading; null for the whole space and in the split form.</param>
/// <param name="Inflated">Whether the set is the inflated rule's.</param>
/// <param name="ErrorLevel">ν, which the bounds' coverage gives up.</param>
/// <param name="FallBack">The opt-in fallback rule's reading (<see cref="ConformalOptions.FallbackShare"/>): whether the scores'
/// error moves too much of the set, or its bounds are untrusted. False unless the rule is asked for.</param>
/// <param name="SplitHalfMass">The split-half diagnostic: the mass, under the full draws' predictive at the case, of the grid
/// values where the plain sets of the draws' two halves disagree (each chain's first and second halves, or alternate draws without
/// chains), the coverage Monte Carlo error could shift here, read at its typical size; conservative, since each half carries √2
/// the run's error. The main session's construction (2026-11-07), not from a paper. NaN on a grid of one value, where it has no
/// width, for the whole space, when the draws do not split, and in the split form.</param>
/// <param name="SplitHalfDisagreements">The grid values where the two halves' plain sets disagree.</param>
public sealed record ConformalSet(int Case, int Stratum, int Units, int Level, double[] Grid, double[]? PValues, int[]? Below,
                                  int[]? Tied, double[]? CandidateScores, double[]? EffectiveSampleSizes,
                                  ConformalMargin? Margin = null, bool Inflated = false, double ErrorLevel = 0.0, bool FallBack = false,
                                  double SplitHalfMass = double.NaN, int SplitHalfDisagreements = 0)
{
    /// <summary>Whether the stratum is too small for the level (k = 0), so the set is the whole space.</summary>
    public bool WholeSpace => PValues is null;

    /// <summary>The guaranteed coverage 1 − k/m under exact scores: at least this, exactly this with distinct scores.</summary>
    public double Coverage => 1.0 - (double)Level / Units;

    /// <summary>
    /// The coverage the set keeps under its scores' Monte Carlo error, indicative since the bounds are not sure ones:
    /// 1 − (k + m̂)/m − ν for the plain rule (<c>coverage_close</c> with <c>card_fragile_le_estimated</c>), 1 − k/m − ν for the
    /// inflated one (<c>coverage_inflated</c>). 1 for the whole space; 1 − k/m with no margin, the split form's, whose guarantee
    /// holds for any fit (<c>split_coverage</c>).
    /// </summary>
    public double CoverageUnderError => WholeSpace ? 1.0 : Margin is null ? Coverage
        : Inflated ? Coverage - ErrorLevel : 1.0 - (double)(Level + Margin.Count) / Units - ErrorLevel;

    /// <summary>Whether grid value <paramref name="g"/> is in the set: the inflated rule's under <see cref="Inflated"/>, the plain
    /// rule's otherwise.</summary>
    public bool Accepts(int g) => Inflated ? AcceptsInflated(g) : AcceptsPlain(g);

    /// <summary>Whether the plain rule accepts grid value <paramref name="g"/>, π(y_g) > k/m.</summary>
    public bool AcceptsPlain(int g) => WholeSpace || Below![g] + Tied![g] + 1 > Level;

    /// <summary>Whether the inflated rule accepts grid value <paramref name="g"/>. It accepts every value the plain rule does,
    /// since raising the candidate's estimate and lowering the units' only adds units at or below it. With no margin, the split
    /// form's, it is the plain rule.</summary>
    public bool AcceptsInflated(int g) => WholeSpace || (Margin is null ? AcceptsPlain(g) : Margin.InflatedRanks[g] > Level);

    /// <summary>
    /// How much more the inflated set holds than the plain one, as a share of the plain set: (|I| − |A|)/|A|, read on the grid's
    /// cells as <see cref="Width"/> reads them, or by count on a grid of one value. The values only the inflated rule accepts are
    /// those whose decision rests on the scores' Monte Carlo error, so this is how far the error moves the set. 0 for the whole
    /// space and when both sets are empty; +∞ when only the inflated one holds anything.
    /// </summary>
    public double Gap
    {
        get
        {
            if (WholeSpace) return 0.0;
            double plain = Extent(AcceptsPlain), inflated = Extent(AcceptsInflated);
            return inflated <= plain ? 0.0 : plain > 0.0 ? (inflated - plain) / plain : double.PositiveInfinity;
        }
    }

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
    public double Width => WholeSpace ? double.PositiveInfinity : Grid.Length == 1 ? 0.0 : Extent(Accepts);

    // The accepted values' cells summed, or their count on a grid of one value.
    private double Extent(Func<int, bool> accepts)
    {
        if (Grid.Length == 1) return accepts(0) ? 1.0 : 0.0;
        double width = 0.0;
        for (int g = 0; g < Grid.Length; g++)
            if (accepts(g)) width += ConformalBayes.Cell(Grid, g);
        return width;
    }
}

/// <summary>The sets, one per test case in the candidates' order, at miscoverage <paramref name="Alpha"/>, from the full step or the
/// split form's density score (<paramref name="Form"/>).</summary>
public sealed record ConformalPrediction(IReadOnlyList<ConformalSet> Sets, double Alpha, ConformalForm Form = ConformalForm.Full);

/// <summary>One candidate's conformity scores (<see cref="ConformalBayes.Scores"/>): the log scores of its stratum's training
/// units, in <paramref name="Units"/>' order, the candidate's own, the add-one-in weights' Kish ESS, and each score's Monte Carlo
/// standard error on the log scale (one error, not the bound, which is <see cref="ConformalMargin.Multiple"/> of them). A
/// candidate outside every draw's support scores −∞, with ESS 0, and the units' scores and every error are NaN: its weights do
/// not exist.</summary>
public sealed record ConformalScores(int[] Units, double[] Scores, double Candidate, double EffectiveSampleSize,
                                     double[] Errors, double CandidateError);

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
/// <item><b>The scores' error</b> (CB-3): each score's relative variance is Σ_s ψ_s² over its draws' influences
/// ψ_s = w̃_s(p_s/p̂ − 1) (<see cref="PredictiveChecks"/>' influence), over the relative efficiency of the candidate's density
/// along the chains, the add-one-in ratio's; two more products, P·W̃² and P²·W̃², give it for every unit. Its log-scale error
/// times z is the bound (<see cref="ConformalOptions.ErrorLevel"/>). From the bounds: the inflated rule, the margin and the
/// fallback signal (<see cref="ConformalMargin"/>).</item>
/// <item><b>Refused</b> (CB-4): leave-unit-out terms (<see cref="PointwiseLogDensity.LeavesUnitsOut"/>, the collapsed
/// predictive), whose rest excludes the candidate; a latent of one observation's term conditioned on or mixed, which reads
/// the data units at their retained latents and the candidate otherwise and undercovers (<c>mixed_counterexample</c>); a
/// candidate whose latents are treated otherwise than the training units'; draws that do not align.</item>
/// </list>
/// The guarantee is exact for exact posterior weights; from Monte Carlo draws it holds up to the estimates' error, which the
/// margin reads. A candidate whose density is zero under every draw scores −∞ and takes π = 1/m. The full step's predictive
/// distribution (<see cref="Distribution"/>) is an opt-in beside the split form's default. Thread-safe: it reads its inputs and
/// writes only its own buffers.
/// </summary>
public static class ConformalBayes
{
    // Candidates scored per product: a multiple of the product's panel, so only a case's last batch has a partial panel.
    private const int Batch = 256;

    // A product entry below this share of its unit's largest likelihood is summed again in logs: above it, the terms that
    // underflow in the product weigh less than S·2⁻¹⁰²²/10⁻²⁰⁰ of the entry.
    private const double Floor = 1e-200;

    // An entry below this share has its error summed in logs: the error's products hold squared likelihoods, which underflow
    // past 10⁻¹⁵⁴.
    private const double ErrorFloor = 1e-100;

    /// <summary>
    /// The full step at every test case of <paramref name="candidates"/>: π(y) on each case's grid within its stratum, the set,
    /// each candidate's ESS, and CB-3's margin, inflated rule and fallback signal (see <see cref="ConformalBayes"/>).
    /// </summary>
    public static ConformalPrediction FullStep(PointwiseLogDensity training, CandidateLogDensity candidates, ConformalOptions? options = null)
    {
        options ??= new ConformalOptions();
        double alpha = options.Alpha, nu = options.ErrorLevel;
        if (!(alpha > 0.0 && alpha < 1.0)) throw new ArgumentOutOfRangeException(nameof(options), "α must lie in (0, 1).");
        if (!(nu > 0.0 && nu < 1.0)) throw new ArgumentOutOfRangeException(nameof(options), "ν must lie in (0, 1).");
        if (options.FallbackShare is double share && !(share >= 0.0)) throw new ArgumentOutOfRangeException(nameof(options), "The fallback share must be non-negative.");
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
                foreach (int t in cases)
                    sets[t] = new ConformalSet(t, h, m, k, candidates.Grids[t], null, null, null, null, null, null, options.Inflate, nu);
                continue;
            }
            var scorer = new StratumScorer(training.Values, rows, logWeights, members, training.ChainLengths, k, Multiple(nu, m),
                                           options.MarginBounds == MarginBounds.Carried);
            (StratumScorer A, StratumScorer B)? halves = Halves(training.Values, rows, logWeights, members, training.ChainLengths, k);
            foreach (int t in cases)
            {
                int g = candidates.Grids[t].Length, offset = candidates.Offset(t);
                var r = new Reading(g);
                Reading? ra = halves is null ? null : new Reading(g), rb = halves is null ? null : new Reading(g);
                for (int b0 = 0; b0 < g; b0 += Batch)
                {
                    int count = Math.Min(Batch, g - b0);
                    scorer.Rank(candidates.Values, offset + b0, count, r.Slice(b0));
                    if (halves is not { } h2) continue;
                    h2.A.Rank(candidates.Values, offset + b0, count, ra!.Value.Slice(b0));
                    h2.B.Rank(candidates.Values, offset + b0, count, rb!.Value.Slice(b0));
                }
                (double splitMass, int disagreements) = SplitHalf(candidates.Grids[t], r, ra, rb, k);
                var p = new double[g];
                for (int j = 0; j < g; j++) p[j] = (r.Below[j] + r.Tied[j] + 1) / (double)m;
                var margin = new ConformalMargin(scorer.Multiple, r.CandidateErrors, r.UniformErrors, r.Fragile, r.InflatedRanks, r.ParetoK,
                                                 PredictiveChecks.ParetoKThreshold(Math.Max(2, rows.Length)), options.MarginBounds);
                var set = new ConformalSet(t, h, m, k, candidates.Grids[t], p, r.Below, r.Tied, r.Own, r.Ess, margin, options.Inflate, nu,
                                           SplitHalfMass: splitMass, SplitHalfDisagreements: disagreements);
                sets[t] = options.FallbackShare is double asked ? set with { FallBack = set.Gap > asked || margin.Untrusted(k) > 0 } : set;
            }
        }
        return new ConformalPrediction(sets, alpha);
    }

    /// <summary>
    /// The full step's conformal predictive distribution, an opt-in beside the split form's default (the user's decision of
    /// 2026-11-08; <see cref="SplitConformal.Distribution"/>). The score is the in-sample predictive CDF under the add-one-in weights,
    /// each unit's F̃ᵢ = Σ_s w̃_s F(y_i | x_i, θ_s) and the candidate's F̃_c = Σ_s w̃_s F(y | x_t, θ_s), read by both tails from the
    /// seams' log CDF (<see cref="PointwiseLogDensity.LogCdf"/>, <see cref="CandidateLogDensity.LogCdf"/>), two (n × S)·(S × G)
    /// products a batch. It spends no data on calibration, but it is not monotone in y: it falls where the draws' predictives at x
    /// disagree (slice 4's reading). So each case is read along its grid, and one whose counts fall is flagged with how far
    /// (<see cref="ConformalDistribution.Fall"/>, <see cref="ConformalDistribution.ScoreFall"/>): it keeps its p-values and its set,
    /// with no distribution. A conformal p-value needs no monotone score, so for draws from the posterior the PIT holds on every case,
    /// flagged or not, up to their Monte Carlo error as the full step's sets do. Draws that are no posterior, a wrong fit's or two
    /// chains stuck apart, break the score's symmetry, and the split form is the guarantee there. No margin is read for this score.
    /// CB-4's refusals and alignment are the full step's.
    /// </summary>
    public static ConformalDistributions Distribution(PointwiseLogDensity training, CandidateLogDensity candidates)
    {
        Check(training, candidates);
        double[,] unitLogCdf = training.LogCdf
            ?? throw new ArgumentException("The training density carries no log CDF; the predictive distribution's score is the in-sample predictive CDF.", nameof(training));
        double[,] candidateLogCdf = candidates.LogCdf
            ?? throw new ArgumentException("The candidates carry no log CDF; the predictive distribution's score is the in-sample predictive CDF.", nameof(candidates));
        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(training);
        Dictionary<int, List<int>> units = UnitsByStratum(training.Units);
        var casesOf = new SortedDictionary<int, List<int>>();
        for (int t = 0; t < candidates.Cases; t++)
        {
            int h = candidates.StratumOf[t];
            if (!casesOf.TryGetValue(h, out List<int>? list)) casesOf.Add(h, list = []);
            list.Add(t);
        }
        var distributions = new ConformalDistribution[candidates.Cases];
        foreach ((int h, List<int> cases) in casesOf)
        {
            int[] members = units.TryGetValue(h, out List<int>? found) ? [.. found] : [];
            var scorer = new CdfScorer(unitLogCdf, rows, logWeights, members);
            foreach (int t in cases)
            {
                int g = candidates.Grids[t].Length, offset = candidates.Offset(t);
                var below = new int[g];
                var tied = new int[g];
                var own = new CdfScore[g];
                for (int b0 = 0; b0 < g; b0 += Batch)
                    scorer.Score(candidates.Values, candidateLogCdf, offset + b0, Math.Min(Batch, g - b0), below, tied, own, b0);
                distributions[t] = new ConformalDistribution(t, h, members.Length + 1, candidates.Grids[t], below, tied, ConformalForm.Full,
                                                             CdfScores.Fall(below, tied), CdfScores.ScoreFall(own));
            }
        }
        return new ConformalDistributions(distributions, ConformalForm.Full);
    }

    /// <summary>
    /// The conformity scores at candidate column <paramref name="candidate"/>: every training unit of its case's stratum, the
    /// candidate's own, the add-one-in weights' ESS and each score's Monte Carlo standard error, computed as
    /// <see cref="FullStep"/> computes them, bit for bit.
    /// </summary>
    public static ConformalScores Scores(PointwiseLogDensity training, CandidateLogDensity candidates, int candidate)
    {
        Check(training, candidates);
        if ((uint)candidate >= (uint)candidates.Values.GetLength(1)) throw new ArgumentOutOfRangeException(nameof(candidate));
        int t = 0;
        while (t + 1 < candidates.Cases && candidates.Offset(t + 1) <= candidate) t++;
        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(training);
        int[] members = UnitsByStratum(training.Units).TryGetValue(candidates.StratumOf[t], out List<int>? found) ? [.. found] : [];
        var scorer = new StratumScorer(training.Values, rows, logWeights, members, training.ChainLengths, members.Length + 1, 1.0);
        var r = new Reading(1);
        scorer.Weigh(candidates.Values, candidate, 1, r);
        var scores = new double[members.Length];
        var errors = new double[members.Length];
        scorer.Multiply(1);
        scorer.Unit(candidates.Values, candidate, r, scores, errors);
        return new ConformalScores(members, scores, r.Own[0], r.Ess[0], errors, r.CandidateErrors[0]);
    }

    // The two halves of the draws' rows, each chain's first and second halves where there are chains (an odd middle draw left out),
    // alternate rows otherwise, each scored in the plain mode; null when a half would be empty, or would give some unit zero density
    // under every one of its draws: the diagnostic is then NaN, and the step goes on.
    private static (StratumScorer A, StratumScorer B)? Halves(double[,] ll, int[] rows, double[]? logWeights, int[] members,
                                                             IReadOnlyList<int>? chains, int level)
    {
        var a = new List<int>();
        var b = new List<int>();
        if (chains is null)
            for (int d = 0; d < rows.Length; d++) (d % 2 == 0 ? a : b).Add(d);
        else
        {
            int start = 0;
            foreach (int length in chains)
            {
                int half = length / 2;
                for (int d = 0; d < half; d++) { a.Add(start + d); b.Add(start + half + d); }
                start += length;
            }
        }
        if (a.Count == 0 || b.Count == 0) return null;
        StratumScorer Half(List<int> at)
            => new(ll, [.. at.Select(d => rows[d])], logWeights is null ? null : [.. at.Select(d => logWeights[d])], members, null, level, 1.0, plain: true);
        try
        {
            return (Half(a), Half(b));
        }
        catch (ArgumentException)
        {
            return null;
        }
    }

    // The split-half reading of one case: the full draws' predictive mass on the grid values where the halves' plain sets disagree,
    // and how many there are.
    private static (double Mass, int Disagreements) SplitHalf(double[] grid, Reading full, Reading? a, Reading? b, int level)
    {
        if (a is not { } ra || b is not { } rb) return (double.NaN, 0);
        double mass = 0.0;
        int disagreements = 0;
        for (int g = 0; g < grid.Length; g++)
        {
            if (ra.Below[g] + ra.Tied[g] + 1 > level == rb.Below[g] + rb.Tied[g] + 1 > level) continue;
            disagreements++;
            if (grid.Length > 1) mass += Math.Exp(full.LogPredictive[g]) * Cell(grid, g);
        }
        return (grid.Length > 1 ? mass : double.NaN, disagreements);
    }

    /// <summary>Grid value <paramref name="g"/>'s cell: the half-gaps to its neighbours, the end values the half-gap on their inner
    /// side twice.</summary>
    internal static double Cell(double[] grid, int g)
    {
        int n = grid.Length;
        double left = g > 0 ? grid[g] - grid[g - 1] : grid[1] - grid[0];
        double right = g < n - 1 ? grid[g + 1] - grid[g] : grid[n - 1] - grid[n - 2];
        return 0.5 * (left + right);
    }

    // z = Φ⁻¹(1 − ν/(2m)): the multiple of a score's Monte Carlo standard error that every one of a stratum's m units stays within
    // with chance 1 − ν, by the union of their normal tails.
    private static double Multiple(double nu, int m) => Normal.StandardQuantile(1.0 - nu / (2.0 * m));

    /// <summary>
    /// The count of FragileEstimated units of <paramref name="scores"/> under one bound <paramref name="error"/> for every unit
    /// (<c>card_fragile_le_estimated</c> read with uniform bounds): those the inflated rule accepts, more than k units with
    /// s_j − ε ≤ s_u + ε, with at most k units j = u or s_j + 4ε ≤ s_u. The predicates are the Lean definitions', each monotone in
    /// s_j, so on the sorted scores two searches read them. <paramref name="scores"/> is sorted in place.
    /// </summary>
    internal static int EstimatedBand(Span<double> scores, double error, int k)
    {
        scores.Sort();
        int m = scores.Length, count = 0;
        double shift = 2 * error + 2 * error;
        for (int u = 0; u < m; u++)
        {
            double su = scores[u], raised = su + error;
            int accepted = CountWhere(scores, v => v - error <= raised);
            if (accepted <= k) continue;
            int below = CountWhere(scores, v => v + shift <= su);
            int rank = su + shift <= su ? below : below + 1;
            if (rank <= k) count++;
        }
        return count;
    }

    /// <summary>
    /// The count of FragileEstimated units under per-unit bounds <paramref name="errors"/>: the inflated rule's side read exactly
    /// on the sorted s_j − ε_j, and rankShift's, s_j + (2ε_j + 2ε_u) ≤ s_u, read as s_j + 2ε_j ≤ s_u − 2ε_u on the sorted
    /// s_j + 2ε_j, which agrees with the literal comparison except at rounding distance.
    /// </summary>
    internal static int EstimatedBand(ReadOnlySpan<double> scores, ReadOnlySpan<double> errors, int k)
    {
        int m = scores.Length, count = 0;
        var lowered = new double[m];
        var raised = new double[m];
        for (int j = 0; j < m; j++)
        {
            lowered[j] = scores[j] - errors[j];
            raised[j] = scores[j] + 2 * errors[j];
        }
        Array.Sort(lowered);
        Array.Sort(raised);
        for (int u = 0; u < m; u++)
        {
            double su = scores[u], eu = errors[u], top = su + eu;
            if (CountWhere(lowered, v => v <= top) <= k) continue;
            double floor = su - 2 * eu;
            int below = CountWhere(raised, v => v <= floor);
            int rank = su + 2 * eu <= floor ? below : below + 1;
            if (rank <= k) count++;
        }
        return count;
    }

    // The leading entries of sorted values for which a predicate monotone in the value holds.
    private static int CountWhere(ReadOnlySpan<double> sorted, Func<double, bool> holds)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (holds(sorted[mid])) lo = mid + 1;
            else hi = mid;
        }
        return lo;
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

    internal static bool SameLatents(IReadOnlyList<LatentReport> a, IReadOnlyList<LatentReport> b)
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

    internal static bool SameSequence<T>(IReadOnlyList<T>? a, IReadOnlyList<T>? b) where T : IEquatable<T>
    {
        if (a is null || b is null) return a is null && b is null;
        if (a.Count != b.Count) return false;
        for (int i = 0; i < a.Count; i++)
            if (!a[i].Equals(b[i])) return false;
        return true;
    }

    /// <summary>A case's readings per grid value, or a batch's slice of them.</summary>
    internal readonly struct Reading
    {
        public Reading(int count)
            : this(new int[count], new int[count], new double[count], new double[count], new double[count], new double[count],
                   new int[count], new int[count], new double[count], new double[count], new double[count], 0) { }

        private Reading(int[] below, int[] tied, double[] own, double[] ess, double[] candidateErrors, double[] uniformErrors,
                        int[] fragile, int[] inflatedRanks, double[] paretoK, double[] efficiency, double[] logPredictive, int start)
        {
            Below = below; Tied = tied; Own = own; Ess = ess; CandidateErrors = candidateErrors; UniformErrors = uniformErrors;
            Fragile = fragile; InflatedRanks = inflatedRanks; ParetoK = paretoK; Efficiency = efficiency; LogPredictive = logPredictive;
            Start = start;
        }

        public int[] Below { get; }
        public int[] Tied { get; }
        public double[] Own { get; }
        public double[] Ess { get; }
        public double[] CandidateErrors { get; }
        public double[] UniformErrors { get; }
        public int[] Fragile { get; }
        public int[] InflatedRanks { get; }
        public double[] ParetoK { get; }
        public double[] Efficiency { get; }   // the relative efficiency of each candidate's density along the chains
        public double[] LogPredictive { get; }   // ln of the draws' predictive density at each candidate, ln Σ_s W_s p(y | θ_s)

        /// <summary>The first grid value of this slice.</summary>
        public int Start { get; }

        public Reading Slice(int start) => new(Below, Tied, Own, Ess, CandidateErrors, UniformErrors, Fragile, InflatedRanks, ParetoK, Efficiency, LogPredictive, start);
    }

    /// <summary>
    /// One stratum's training units against batches of candidates. It holds each unit's likelihoods over the draws of positive
    /// weight scaled by their largest, P_is = p(z_i | θ_s)/max_s, and their squares, row by row, the products' left operands,
    /// and for each batch the candidates' normalized add-one-in weights and their squares draw by draw, the right ones. In the
    /// plain mode, which scores a half of the draws for the split-half reading, it takes the scores and their ranks alone: no
    /// squares, errors, bands or k̂. Not thread-safe: one per stratum and call.
    /// </summary>
    internal sealed class StratumScorer
    {
        private readonly double[,] _ll;
        private readonly int[] _rows, _members;
        private readonly double[]? _logWeights;
        private readonly IReadOnlyList<int>? _chains;
        private readonly double[] _left, _leftSquared, _top;   // P and P², units × draws; each unit's largest log term
        private readonly int _s, _level;
        private readonly bool _plain;
        private double[] _weights = [], _weightsSquared = [], _likelihoods = [], _logNorm = [], _candidateTop = [], _squareSums = [],
                         _ownFirst = [], _ownSecond = [], _product = [], _first = [], _second = [], _column = [];
        private PackedColumns? _right, _rightSquared;

        public StratumScorer(double[,] ll, int[] rows, double[]? logWeights, int[] members, IReadOnlyList<int>? chains, int level, double multiple,
                             bool carried = false, bool plain = false)
        {
            _plain = plain;
            _ll = ll;
            _rows = rows;
            _logWeights = logWeights;
            _members = members;
            _chains = chains;
            _level = level;
            Multiple = multiple;
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
            _leftSquared = plain ? [] : new double[(long)n * s];
            const int tile = 32;
            for (int d0 = 0; d0 < s; d0 += tile)
                for (int i0 = 0; i0 < n; i0 += tile)
                    for (int i = i0; i < Math.Min(n, i0 + tile); i++)
                    {
                        int u = members[i];
                        double top = _top[i];
                        for (int d = d0; d < Math.Min(s, d0 + tile); d++)
                        {
                            double p = Math.Exp(ll[rows[d], u] - top);
                            _left[(long)i * s + d] = p;
                            if (!plain) _leftSquared[(long)i * s + d] = p * p;
                        }
                    }
            if (carried) _selfErrors = SelfErrors();
        }

        private readonly double[]? _selfErrors;   // each unit's standard error as the candidate, for the carried bounds
        private double[] _bounds = [];

        // Each unit's Monte Carlo standard error as the candidate: its own score under weights W_s·p(z_i | θ_s), from its row of
        // P, over the relative efficiency of its likelihood along the chains, as PSIS-LOO reads each unit's.
        private double[] SelfErrors()
        {
            int n = _members.Length, s = _s;
            var errors = new double[n];
            var row = new double[s];
            var weight = new double[s];
            for (int d = 0; d < s; d++) weight[d] = _logWeights is null ? 1.0 : Math.Exp(_logWeights[d]);
            for (int i = 0; i < n; i++)
            {
                double total = 0.0, first = 0.0;
                for (int d = 0; d < s; d++)
                {
                    double p = _left[(long)i * s + d], w = weight[d] * p;
                    row[d] = p;
                    total += w;
                    first += w * p;
                }
                double score = first / total, relative = 0.0;
                for (int d = 0; d < s; d++)
                {
                    double t = weight[d] * row[d] / total * (row[d] / score - 1.0);
                    relative += t * t;
                }
                double efficiency = _chains is null ? 1.0 : PredictiveChecks.RelativeEfficiency(row, _chains);
                errors[i] = PredictiveChecks.LogScaleError(relative, efficiency);
            }
            return errors;
        }

        /// <summary>z: each bound is this many Monte Carlo standard errors.</summary>
        public double Multiple { get; }

        /// <summary>The readings of <paramref name="count"/> candidates from column <paramref name="first"/> into
        /// <paramref name="into"/> from its start. Its steps are apart for the cost row: <see cref="Weigh"/>, <see cref="Multiply"/>,
        /// <see cref="Count"/>.</summary>
        public void Rank(double[,] candidates, int first, int count, Reading into)
        {
            Weigh(candidates, first, count, into);
            Multiply(count);
            if (_plain) CountPlain(count, into);
            else Count(candidates, first, count, into);
        }

        /// <summary>The training units' scores and their errors' sums against the batch the last <see cref="Weigh"/> prepared:
        /// P·W̃, P·W̃² and P²·W̃², units × candidates.</summary>
        public void Multiply(int count)
        {
            int n = _members.Length;
            if (n == 0) return;
            if (_product.Length < n * count)
            {
                _product = new double[n * count];
                _first = new double[n * count];
                _second = new double[n * count];
            }
            DenseProduct.Multiply(_left, n, _right!, _product);
            if (_plain) return;
            DenseProduct.Multiply(_left, n, _rightSquared!, _first);
            DenseProduct.Multiply(_leftSquared, n, _rightSquared!, _second);
        }

        // The plain mode's counts: each candidate's below and tied units, as the full step's are taken, scores read from the product
        // with the same floor.
        private void CountPlain(int count, Reading into)
        {
            int at = into.Start;
            Span<int> below = into.Below.AsSpan(at, count), tied = into.Tied.AsSpan(at, count);
            ReadOnlySpan<double> own = into.Own.AsSpan(at, count);
            below.Clear();
            tied.Clear();
            for (int i = 0; i < _members.Length; i++)
            {
                ReadOnlySpan<double> row = _product.AsSpan(i * count, count);
                double top = _top[i];
                for (int b = 0; b < count; b++)
                {
                    if (double.IsNegativeInfinity(own[b])) continue;
                    double acc = row[b];
                    double score = acc >= Floor ? top + Math.Log(acc) : LogScoreAndVariance(_batch!, _batchFirst + b, b, _members[i]).Score;
                    if (score < own[b]) below[b]++;
                    else if (score == own[b]) tied[b]++;
                }
            }
        }

        private double[,]? _batch;   // the last weighed batch's candidates and first column, for the plain mode's log path
        private int _batchFirst;

        /// <summary>Each candidate's below and tied counts, inflated rank, uniform bound, estimated band and, where the inflated
        /// rule accepts, its weights' k̂, from the products, against its own score and bound.</summary>
        public void Count(double[,] candidates, int first, int count, Reading into)
        {
            int n = _members.Length, at = into.Start;
            Span<int> below = into.Below.AsSpan(at, count), tied = into.Tied.AsSpan(at, count), inflated = into.InflatedRanks.AsSpan(at, count);
            Span<double> own = into.Own.AsSpan(at, count), ownError = into.CandidateErrors.AsSpan(at, count), uniform = into.UniformErrors.AsSpan(at, count);
            below.Clear();
            tied.Clear();
            for (int b = 0; b < count; b++)
            {
                inflated[b] = 1;   // the candidate itself
                uniform[b] = ownError[b];
            }
            for (int i = 0; i < n; i++)
            {
                Span<double> row = _product.AsSpan(i * count, count);
                Span<double> a1 = _first.AsSpan(i * count, count);
                ReadOnlySpan<double> a2 = _second.AsSpan(i * count, count);
                double top = _top[i];
                for (int b = 0; b < count; b++)
                {
                    double candidate = own[b];
                    if (double.IsNegativeInfinity(candidate)) continue;   // outside every draw's support: π = 1/m
                    double acc = row[b], score, relative;
                    if (acc >= ErrorFloor)
                    {
                        score = top + Math.Log(acc);
                        relative = a2[b] / (acc * acc) - 2.0 * a1[b] / acc + _squareSums[b];
                    }
                    else
                    {
                        (double exact, relative) = LogScoreAndVariance(candidates, first + b, b, _members[i]);
                        score = acc >= Floor ? top + Math.Log(acc) : exact;
                    }
                    // An error that does not exist (every weighted term of the unit at −∞) is taken as unbounded, so the inflated
                    // rule still holds every unit the plain rule counts.
                    double error = double.IsNaN(relative) ? double.PositiveInfinity
                                 : Multiple * PredictiveChecks.LogScaleError(Math.Max(0.0, relative), into.Efficiency[at + b]);
                    if (score < candidate) below[b]++;
                    else if (score == candidate) tied[b]++;
                    if (score - error <= candidate + ownError[b]) inflated[b]++;
                    if (error > uniform[b]) uniform[b] = error;
                    a1[b] = _selfErrors is null ? error : Math.Max(error, Multiple * _selfErrors[i]);   // the carried bound, in the spent sum's place
                    row[b] = score;
                }
            }
            Span<double> column = _column.Length >= n + 1 ? _column.AsSpan(0, n + 1) : (_column = new double[n + 1]).AsSpan();
            for (int b = 0; b < count; b++)
            {
                if (double.IsNegativeInfinity(own[b]))
                {
                    into.Fragile[at + b] = -1;
                    into.ParetoK[at + b] = double.NaN;
                    continue;
                }
                for (int i = 0; i < n; i++) column[i] = _product[i * count + b];
                column[n] = own[b];
                if (_selfErrors is null) into.Fragile[at + b] = EstimatedBand(column, uniform[b], _level);
                else
                {
                    Span<double> bounds = _bounds.Length >= n + 1 ? _bounds.AsSpan(0, n + 1) : (_bounds = new double[n + 1]).AsSpan();
                    for (int i = 0; i < n; i++) bounds[i] = _first[i * count + b];
                    bounds[n] = ownError[b];
                    into.Fragile[at + b] = EstimatedBand(column, bounds, _level);
                }
                into.ParetoK[at + b] = inflated[b] > _level ? ParetoK(candidates, first + b, into.Efficiency[at + b]) : double.NaN;
            }
        }

        /// <summary>One candidate's unit scores and Monte Carlo standard errors, by the same arithmetic as <see cref="Count"/>.</summary>
        public void Unit(double[,] candidates, int column, Reading from, Span<double> scores, Span<double> errors)
        {
            bool outside = double.IsNegativeInfinity(from.Own[0]);
            for (int i = 0; i < _members.Length; i++)
            {
                if (outside)
                {
                    scores[i] = errors[i] = double.NaN;
                    continue;
                }
                double acc = _product[i], relative;
                if (acc >= ErrorFloor)
                {
                    scores[i] = _top[i] + Math.Log(acc);
                    relative = _second[i] / (acc * acc) - 2.0 * _first[i] / acc + _squareSums[0];
                }
                else
                {
                    (double exact, relative) = LogScoreAndVariance(candidates, column, 0, _members[i]);
                    scores[i] = acc >= Floor ? _top[i] + Math.Log(acc) : exact;
                }
                errors[i] = Multiple * PredictiveChecks.LogScaleError(Math.Max(0.0, relative), from.Efficiency[0]);
            }
        }

        /// <summary>The batch's normalized add-one-in weights and their squares, draws × candidates, packed as the products' right
        /// operands; each candidate's scaled likelihoods, log normalizer, own log score and bound, ESS and relative efficiency.</summary>
        public void Weigh(double[,] candidates, int first, int count, Reading into)
        {
            int s = _s, size = s * count, at = into.Start;
            _batch = candidates;
            _batchFirst = first;
            if (_weights.Length < size)
            {
                _weights = new double[size];
                _weightsSquared = new double[size];
                _likelihoods = new double[size];
            }
            if (_logNorm.Length < count)
            {
                _logNorm = new double[count];
                _candidateTop = new double[count];
                _squareSums = new double[count];
                _ownFirst = new double[count];
                _ownSecond = new double[count];
            }
            Span<double> own = into.Own.AsSpan(at, count), ess = into.Ess.AsSpan(at, count);
            Span<double> top = _candidateTop.AsSpan(0, count), norm = _logNorm.AsSpan(0, count);
            Span<double> w = _weights.AsSpan(0, size), w2 = _weightsSquared.AsSpan(0, size), p = _likelihoods.AsSpan(0, size);
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
            ess.Clear();
            for (int d = 0; d < s; d++)
            {
                Span<double> wd = w.Slice(d * count, count), w2d = w2.Slice(d * count, count);
                for (int b = 0; b < count; b++)
                {
                    wd[b] *= totals[b];
                    w2d[b] = wd[b] * wd[b];
                    ess[b] += w2d[b];
                }
            }
            Span<double> squares = _squareSums.AsSpan(0, count);
            for (int b = 0; b < count; b++)
            {
                squares[b] = ess[b];   // Σ w̃², the relative variance's constant term
                ess[b] = ess[b] > 0.0 ? 1.0 / ess[b] : 0.0;
            }
            // The draws' predictive density at each candidate: ln Σ_s p(y | θ_s) − ln S at equal weights, ln Σ_s W_s p(y | θ_s) with
            // the normalized weights.
            for (int b = 0; b < count; b++) into.LogPredictive[at + b] = _logWeights is null ? norm[b] - Math.Log(s) : norm[b];
            Span<double> ownFirst = _ownFirst.AsSpan(0, count), ownSecond = _ownSecond.AsSpan(0, count);
            OwnSums(w, w2, p, s, count, own, ownFirst, ownSecond);
            Span<double> errors = into.CandidateErrors.AsSpan(at, count), efficiency = into.Efficiency.AsSpan(at, count);
            for (int b = 0; b < count; b++)
            {
                efficiency[b] = _plain ? 1.0 : Efficiency(p, count, b);
                if (double.IsNegativeInfinity(top[b]) || double.IsNegativeInfinity(norm[b]))
                {
                    own[b] = double.NegativeInfinity;
                    errors[b] = double.NaN;
                    continue;
                }
                double acc = own[b], rel;
                if (acc >= ErrorFloor)
                {
                    rel = ownSecond[b] / (acc * acc) - 2.0 * ownFirst[b] / acc + squares[b];
                    own[b] = top[b] + Math.Log(acc);
                }
                else
                {
                    (double exact, rel) = LogScoreAndVariance(candidates, first + b, b, -1);
                    own[b] = acc >= Floor ? top[b] + Math.Log(acc) : exact;
                }
                errors[b] = Multiple * PredictiveChecks.LogScaleError(Math.Max(0.0, rel), efficiency[b]);
            }
            _right ??= new PackedColumns(s, count);
            if (_right.Columns != count) _right = new PackedColumns(s, count);
            _right.Fill(w);
            if (_plain) return;
            _rightSquared ??= new PackedColumns(s, count);
            if (_rightSquared.Columns != count) _rightSquared = new PackedColumns(s, count);
            _rightSquared.Fill(w2);
        }

        // The relative efficiency of candidate b's density along the chains, the add-one-in ratio's (PSIS reads the left-out
        // unit's likelihood, its ratio's, the same way); 1 for independent or weighted draws.
        private double Efficiency(ReadOnlySpan<double> p, int count, int b)
        {
            if (_chains is null) return 1.0;
            var values = new double[_s];
            for (int d = 0; d < _s; d++) values[d] = p[d * count + b];
            return PredictiveChecks.RelativeEfficiency(values, _chains);
        }

        // k̂ of candidate column `column`'s add-one-in log ratios, ln W_s + ln p(y | θ_s).
        private double ParetoK(double[,] candidates, int column, double efficiency)
        {
            if (_s < 2) return double.PositiveInfinity;
            var ratios = new double[_s];
            for (int d = 0; d < _s; d++) ratios[d] = (_logWeights is null ? 0.0 : _logWeights[d]) + candidates[_rows[d], column];
            return PredictiveChecks.ParetoK(ratios, efficiency);
        }

        // ln Σ_s w̃_s p(z | θ_s) summed in logs, and its relative variance Σ_s ψ_s², for unit `unit`, or the candidate's own column
        // when unit is −1.
        private (double Score, double Relative) LogScoreAndVariance(double[,] candidates, int column, int b, int unit)
        {
            double norm = _logNorm[b];
            double max = double.NegativeInfinity;
            for (int d = 0; d < _s; d++) max = Math.Max(max, LogWeight(candidates, column, norm, d) + Term(candidates, column, unit, d));
            if (double.IsNegativeInfinity(max)) return (max, double.NaN);
            double sum = 0.0;
            for (int d = 0; d < _s; d++) sum += Math.Exp(LogWeight(candidates, column, norm, d) + Term(candidates, column, unit, d) - max);
            double score = max + Math.Log(sum), relative = 0.0;
            for (int d = 0; d < _s; d++)
            {
                double psi = PredictiveChecks.Influence(LogWeight(candidates, column, norm, d), Term(candidates, column, unit, d), score);
                relative += psi * psi;
            }
            return (score, relative);
        }

        private double LogWeight(double[,] candidates, int column, double norm, int d)
            => (_logWeights is null ? 0.0 : _logWeights[d]) + candidates[_rows[d], column] - norm;

        private double Term(double[,] candidates, int column, int unit, int d) => unit < 0 ? candidates[_rows[d], column] : _ll[_rows[d], unit];

        // By fused multiply-adds in draw order, four candidates to a vector: own[b] = Σ_d w̃·p, each lane the product's entry for a
        // unit whose column equals the candidate's, bit for bit; first[b] = Σ_d w̃²·p and second[b] = Σ_d w̃²·p², the error's sums
        // as the other two products take them.
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private static void OwnSums(ReadOnlySpan<double> w, ReadOnlySpan<double> w2, ReadOnlySpan<double> p, int s, int count,
                                    Span<double> own, Span<double> first, Span<double> second)
        {
            ref double wr = ref MemoryMarshal.GetReference(w);
            ref double w2r = ref MemoryMarshal.GetReference(w2);
            ref double pr = ref MemoryMarshal.GetReference(p);
            int b = 0;
            for (; b + 4 <= count; b += 4)
            {
                Vector256<double> acc = Vector256<double>.Zero, a1 = Vector256<double>.Zero, a2 = Vector256<double>.Zero;
                for (int d = 0; d < s; d++)
                {
                    nuint at = (nuint)(d * count + b);
                    Vector256<double> pv = Vector256.LoadUnsafe(ref pr, at), w2v = Vector256.LoadUnsafe(ref w2r, at);
                    acc = Vector256.FusedMultiplyAdd(Vector256.LoadUnsafe(ref wr, at), pv, acc);
                    a1 = Vector256.FusedMultiplyAdd(w2v, pv, a1);
                    a2 = Vector256.FusedMultiplyAdd(w2v, pv * pv, a2);
                }
                acc.CopyTo(own.Slice(b, 4));
                a1.CopyTo(first.Slice(b, 4));
                a2.CopyTo(second.Slice(b, 4));
            }
            for (; b < count; b++)
            {
                double acc = 0.0, a1 = 0.0, a2 = 0.0;
                for (int d = 0; d < s; d++)
                {
                    double pv = p[d * count + b], w2v = w2[d * count + b];
                    acc = Math.FusedMultiplyAdd(w[d * count + b], pv, acc);
                    a1 = Math.FusedMultiplyAdd(w2v, pv, a1);
                    a2 = Math.FusedMultiplyAdd(w2v, pv * pv, a2);
                }
                own[b] = acc;
                first[b] = a1;
                second[b] = a2;
            }
        }
    }
}
