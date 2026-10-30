using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Exact draws from BASS's prior restricted to F (D02, D05): the chains' starts (the stage-4 plan) and the b = 0 member
/// of the likelihood path, which an exact reference rung draws from. A sampler of the posterior's b = 0 end, by
/// rejection from the globally truncated measure, which is exact.
/// <list type="bullet">
/// <item><b>The proposal.</b> M from the count law with its orderings, p(M) ∝ Γ(h₁ + M)/(M!(h₂ + 1)^M) on
/// {0, …, M_max}: the law of λ ~ Gamma(h₁, rate h₂), M | λ ~ Poisson(λ), which D02 integrates out, drawn by inversion
/// and truncated at the cap, which F's cap rule would refuse anyway. Then M terms independently from D01's term measure
/// (<see cref="BassStructurePrior.DrawTerm"/>), and τ from its prior when it is sampled.</item>
/// <item><b>The acceptance.</b> The collection is kept only when it lies in F: no term the intercept (a draw whose
/// parts are all empty), distinct terms, the support and D06's gates as the design admits them, and under the Jeffreys
/// σ² prior no exact fit (<see cref="BassTarget.Refusal"/>). Anything else is drawn again.</item>
/// </list>
/// Its account (<see cref="Attempts"/>, <see cref="Accepted"/>, refusals by reason) is the share of the prior F holds.
/// A draw refuses after <see cref="MaxAttempts"/> attempts. Not thread-safe: it keeps its own design, one per chain.
/// </summary>
public sealed class BassPriorSampler : IReferenceSampler<BassState>
{
    /// <summary>The attempts a draw makes before it refuses, by default.</summary>
    public const long DefaultMaxAttempts = 1_000_000;

    private readonly BassModel _model;
    private readonly BassDesign _design;
    private readonly double[] _countCdf;   // the truncated count law's CDF over M = 0, …, M_max

    /// <summary>Exact draws from <paramref name="model"/>'s prior on F, each refusing after
    /// <paramref name="maxAttempts"/> attempts.</summary>
    public BassPriorSampler(BassModel model, long maxAttempts = DefaultMaxAttempts)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (maxAttempts < 1) throw new ArgumentOutOfRangeException(nameof(maxAttempts), "A draw needs at least one attempt.");
        _model = model;
        _design = model.NewDesign();
        MaxAttempts = maxAttempts;

        int cap = model.MaxBasis;
        var log = new double[cap + 1];
        double max = double.NegativeInfinity;
        for (int m = 0; m <= cap; m++)
        {
            log[m] = model.Prior.LogCollectionWeight(m) - SpecialFunctions.LogGamma(m + 1.0);
            max = Math.Max(max, log[m]);
        }
        _countCdf = new double[cap + 1];
        double total = 0.0;
        for (int m = 0; m <= cap; m++) _countCdf[m] = total += Math.Exp(log[m] - max);
        for (int m = 0; m < cap; m++) _countCdf[m] /= total;
        _countCdf[cap] = 1.0;
    }

    /// <summary>The likelihood path: this samples its b = 0 member, the prior.</summary>
    public TemperingPath Path => TemperingPath.Likelihood;

    /// <summary>The attempts a draw makes before it refuses.</summary>
    public long MaxAttempts { get; }

    /// <summary>Collections drawn so far, accepted ones included.</summary>
    public long Attempts { get; private set; }

    /// <summary>Collections accepted, one per draw.</summary>
    public long Accepted { get; private set; }

    /// <summary>The share of attempts accepted, the prior mass of F against the truncated count law's; NaN before the
    /// first attempt.</summary>
    public double AcceptanceRate => Attempts > 0 ? (double)Accepted / Attempts : double.NaN;

    /// <summary>Refused collections by reason: the design's (repeated terms as numerical, support, numerical, variance
    /// inflation) and the target's exact fits.</summary>
    public Tally<BassRefusal> Refusals { get; } = new();

    /// <summary>Collections refused because a term came out the intercept, every part empty; only fits of several parts
    /// have them.</summary>
    public long InterceptTerms { get; private set; }

    // The account's counters, for a ladder's restart state (BassLadder): the draws depend on the stream alone.
    internal BassPriorSamplerCounts SaveCounts() => new(Attempts, Accepted, Refusals.Copy(), InterceptTerms);

    internal void RestoreCounts(BassPriorSamplerCounts counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        Attempts = counts.Attempts;
        Accepted = counts.Accepted;
        Refusals.Assign(counts.Refusals);
        InterceptTerms = counts.InterceptTerms;
    }

    /// <summary>One exact draw of (structure, τ) from the prior on F, with its evaluation under the model's target.</summary>
    public (BassState State, ITargetEvaluation Evaluation) Draw(Xoshiro256PlusPlus rng)
    {
        (BassState state, BassEvaluation evaluation) = DrawState(rng);
        return (state, evaluation);
    }

    /// <summary><see cref="Draw"/>, its evaluation typed.</summary>
    public (BassState State, BassEvaluation Evaluation) DrawState(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        BassStructurePrior prior = _model.Prior;
        BassTarget target = _model.Target;
        long intercepts = InterceptTerms;
        Tally<BassRefusal> before = Refusals.Copy();
        for (long attempt = 0; attempt < MaxAttempts; attempt++)
        {
            Attempts++;
            var terms = new BassTerm[DrawCount(rng)];
            double logTermPrior = 0.0;
            bool intercept = false;
            for (int j = 0; j < terms.Length && !intercept; j++)
            {
                if (prior.DrawTerm(rng) is not { } term) intercept = true;
                else
                {
                    terms[j] = term;
                    logTermPrior += prior.LogPrior(term);
                }
            }
            if (intercept)
            {
                InterceptTerms++;
                continue;
            }
            if (_design.Load(terms) is { } refusal)
            {
                Refusals.Add(refusal);
                continue;
            }
            var state = new BassState([.. _design.Terms], BassStructureFit.Of(_design, logTermPrior), DrawLogTau(rng));
            BassEvaluation evaluation = target.Evaluate(state);
            if (evaluation.Refusal is { } exact)
            {
                Refusals.Add(exact);
                continue;
            }
            Accepted++;
            return (state, evaluation);
        }
        var refused = Refusals.Entries.Select(e => (e.Category, Count: e.Count - before.Count(e.Category))).Where(e => e.Count > 0);
        throw new InvalidOperationException(
            $"No collection in F after {MaxAttempts} draws from the prior ({InterceptTerms - intercepts} with an intercept term; refusals " +
            $"{string.Join(", ", refused.Select(e => $"{e.Category} {e.Count}"))}): F holds too little of the prior's mass. Loosen the " +
            "support rule or the caps, or move the count prior towards fewer terms.");
    }

    // M from the count law truncated at the cap, by inversion.
    private int DrawCount(Xoshiro256PlusPlus rng)
    {
        double u = rng.NextDouble();
        int lo = 0, hi = _countCdf.Length - 1;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (u < _countCdf[mid]) hi = mid;
            else lo = mid + 1;
        }
        return lo;
    }

    // u = ln τ: an exact draw from τ's prior, redrawn in the far tails where τ leaves the doubles; ln of a held τ. An
    // improper prior (hyper-g at a ≤ 2) has no draw: a chain's start takes ln τ₀, its family's centre, and the reference
    // rung and SMC, which need the draw, refuse it (BassConfig.Validate).
    private double DrawLogTau(Xoshiro256PlusPlus rng)
    {
        ResolvedTauPrior scale = _model.Target.Scale;
        if (scale.HeldTau is double held) return Math.Log(held);
        return scale.IsProper ? scale.SampleLogTau(rng) : Math.Log(scale.UnitInformation);
    }
}

/// <summary>A <see cref="BassPriorSampler"/>'s account: attempts, acceptances, refusals by reason and intercept terms.</summary>
internal sealed record BassPriorSamplerCounts(long Attempts, long Accepted, Tally<BassRefusal> Refusals, long InterceptTerms);
