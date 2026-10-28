using System;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// D01's structure proposals over a <see cref="BassStructurePrior"/>'s term space, given the current structure's
/// <see cref="BassTermCounts"/> (Nott, Kuk and Duc, extended by BASS to several parts):
/// <list type="bullet">
/// <item><b>Birth.</b> For each part, its order from I(K) ∝ w₁ + (terms of that order); for K = 1 its variable
/// uniformly; for K ≥ 2 its variables by weighted sampling without replacement with z_u ∝ w₂ + (terms using u), the
/// set's probability Wallenius'; signs, knots and subsets from the prior. A draw with every part empty is the
/// intercept: a failed birth, as in R.</item>
/// <item><b>Death.</b> A term uniformly (the caller's draw); its reverse is the removed term's birth under the reduced
/// structure's counts, which w₁ > 0 and w₂ > 0 keep positive.</item>
/// <item><b>Change.</b> A term uniformly (the caller's draw), then one of its factors uniformly, whose sign and knot,
/// or subset, is redrawn from the prior.</item>
/// </list>
/// Each move reports its within-move log ratio, ln q(reverse)/q(forward), which a reversible-jump chain adds to its
/// target ratio (the target holding the structure prior) and its move-selection ratio. A birth's and a death's carry the
/// 1/M of the death's pick; with the prior's ratio they make D02's count factor (M + h₁)/((M + 1)(h₂ + 1)). Immutable,
/// with no scratch state, so one serves every chain of a fit.
/// </summary>
public sealed class BassTermProposal
{
    private const int StackLimit = 256;

    /// <summary>Proposals over <paramref name="prior"/>'s term space with the pseudo-counts w₁ and w₂, both positive.</summary>
    public BassTermProposal(BassStructurePrior prior, double orderWeight, double variableWeight)
    {
        ArgumentNullException.ThrowIfNull(prior);
        if (!(orderWeight > 0.0 && double.IsFinite(orderWeight)))
            throw new ArgumentOutOfRangeException(nameof(orderWeight), "w₁ must be positive and finite, or a death could have no reverse birth.");
        if (!(variableWeight > 0.0 && double.IsFinite(variableWeight)))
            throw new ArgumentOutOfRangeException(nameof(variableWeight), "w₂ must be positive and finite, or a death could have no reverse birth.");
        Prior = prior;
        OrderWeight = orderWeight;
        VariableWeight = variableWeight;
    }

    /// <summary>The proposals <paramref name="config"/> sets over <paramref name="prior"/>'s term space.</summary>
    public static BassTermProposal From(BassConfig config, BassStructurePrior prior)
    {
        ArgumentNullException.ThrowIfNull(config);
        return new(prior, config.Structure.OrderWeight, config.Structure.VariableWeight);
    }

    /// <summary>The prior whose term space the proposals range over.</summary>
    public BassStructurePrior Prior { get; }

    /// <summary>w₁, the order pseudo-count.</summary>
    public double OrderWeight { get; }

    /// <summary>w₂, the variable pseudo-count.</summary>
    public double VariableWeight { get; }

    /// <summary>ln q_b(term | counts): the log probability that a birth into the structure <paramref name="counts"/>
    /// describes draws <paramref name="term"/>.</summary>
    public double LogBirth(BassTerm term, BassTermCounts counts) => LogBirth(term, counts, 0);

    /// <summary>The log probability that a birth draws every part empty, the intercept; −∞ for a fit of one part.</summary>
    public double LogEmptyBirth(BassTermCounts counts)
    {
        Check(counts);
        double lq = 0.0;
        for (int index = 0; index < 3; index++)
        {
            var kind = (BassPart)index;
            if (Prior[kind] is not { } part) continue;
            if (part.MinimumOrder > 0) return double.NegativeInfinity;
            lq += LogOrder(part, counts.Orders(kind)[0], counts.Terms);
        }
        return lq;
    }

    /// <summary>
    /// The within-move log ratio of a birth of <paramref name="born"/> into the structure <paramref name="counts"/>
    /// describes, M terms: ln[1/(M + 1)] − ln q_b(born | counts), the reverse death's pick over the birth's draw.
    /// </summary>
    public double LogBirthRatio(BassTerm born, BassTermCounts counts)
        => -Math.Log(counts.Terms + 1.0) - LogBirth(born, counts, 0);

    /// <summary>
    /// The within-move log ratio of a death of <paramref name="dying"/> from the structure <paramref name="counts"/>
    /// describes, M terms with it among them: ln q_b(dying | the counts without it) − ln(1/M).
    /// </summary>
    public double LogDeathRatio(BassTerm dying, BassTermCounts counts)
        => LogBirth(dying, counts, 1) + Math.Log(counts.Terms);

    /// <summary>A birth's term, or null for the empty draw (the intercept, a failed birth).</summary>
    public BassTerm? DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Check(counts);
        Span<int> orders = stackalloc int[3];
        bool empty = true;
        for (int index = 0; index < 3; index++)
        {
            var kind = (BassPart)index;
            if (Prior[kind] is not { } part) continue;
            orders[index] = DrawOrder(rng, part, counts.Orders(kind));
            empty &= orders[index] == 0;
        }
        if (empty) return null;

        HingeFactor[]? continuous = null, functional = null;
        CategoricalFactor[]? categorical = null;
        for (int index = 0; index < 3; index++)
        {
            var kind = (BassPart)index;
            int k = orders[index];
            if (k == 0) continue;
            BassStructurePrior.Part part = Prior[kind]!;
            int[] variables = DrawVariables(rng, part, counts.AppearanceCounts(kind), k);
            if (kind == BassPart.Categorical)
            {
                categorical = new CategoricalFactor[k];
                for (int f = 0; f < k; f++) categorical[f] = DrawSubset(rng, variables[f], part.Sizes[variables[f]]);
                continue;
            }
            var hinges = new HingeFactor[k];
            for (int f = 0; f < k; f++) hinges[f] = DrawHinge(rng, part, variables[f]);
            if (kind == BassPart.Continuous) continuous = hinges;
            else functional = hinges;
        }
        return new BassTerm(continuous, categorical, functional);
    }

    /// <summary>
    /// A change of <paramref name="term"/>: one of its factors uniformly, its sign and knot, or its subset, redrawn
    /// from the prior. <paramref name="logProposalRatio"/> is ln q(back)/q(forward), the redrawn factor's prior
    /// probability before over after: 0 for a hinge, and the prior's ratio cancels it for a subset, so the move's
    /// acceptance is the likelihood ratio alone.
    /// </summary>
    public BassTerm DrawChange(Xoshiro256PlusPlus rng, BassTerm term, out double logProposalRatio)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Prior.Validate(term);
        int continuous = term.Continuous.Count, categorical = term.Categorical.Count;
        int pick = BassDraws.Below(rng, continuous + categorical + term.Functional.Count);
        logProposalRatio = 0.0;
        if (pick < continuous)
        {
            HingeFactor[] hinges = term.Continuous.ToArray();
            hinges[pick] = DrawHinge(rng, Prior[BassPart.Continuous]!, hinges[pick].Variable);
            return new BassTerm(hinges, term.Categorical, term.Functional);
        }
        pick -= continuous;
        if (pick < categorical)
        {
            BassStructurePrior.Part part = Prior[BassPart.Categorical]!;
            CategoricalFactor[] subsets = term.Categorical.ToArray();
            CategoricalFactor old = subsets[pick];
            subsets[pick] = DrawSubset(rng, old.Variable, part.Sizes[old.Variable]);
            double[] log = part.LogSubset![old.Variable];
            logProposalRatio = log[old.Levels.Count] - log[subsets[pick].Levels.Count];
            return new BassTerm(term.Continuous, subsets, term.Functional);
        }
        pick -= categorical;
        HingeFactor[] outputs = term.Functional.ToArray();
        outputs[pick] = DrawHinge(rng, Prior[BassPart.Functional]!, outputs[pick].Variable);
        return new BassTerm(term.Continuous, term.Categorical, outputs);
    }

    // ln q_b(term) under the counts less `shift` copies of the term itself: 0 for a birth into them, 1 for the reverse
    // of a death from them.
    private double LogBirth(BassTerm term, BassTermCounts counts, int shift)
    {
        Check(counts);
        Prior.Validate(term);
        int terms = counts.Terms - shift;
        if (terms < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
        double lq = 0.0;
        for (int index = 0; index < 3; index++)
        {
            var kind = (BassPart)index;
            if (Prior[kind] is not { } part) continue;
            int k = BassStructurePrior.Order(term, kind);
            int count = counts.Orders(kind)[k] - shift;
            if (count < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
            lq += LogOrder(part, count, terms);
            if (k == 0) continue;
            lq += LogVariables(term, kind, part, counts.AppearanceCounts(kind), shift);
            lq += part.LogFactors(term);
        }
        return lq;
    }

    // ln I(K) = ln(w₁ + c_K) − ln Σⱼ(w₁ + cⱼ); every term takes exactly one order in the part, so Σⱼ cⱼ = M.
    private double LogOrder(BassStructurePrior.Part part, int count, int terms)
        => Math.Log(OrderWeight + count) - Math.Log(OrderWeight * part.OrderCount + terms);

    // The variable set's log probability: uniform over the p variables for one, Wallenius' for several.
    private double LogVariables(BassTerm term, BassPart kind, BassStructurePrior.Part part, int[] appearances, int shift)
    {
        int p = part.Sizes.Length, k = BassStructurePrior.Order(term, kind);
        Span<int> set = k <= StackLimit ? stackalloc int[k] : new int[k];
        int f = 0;
        foreach (int v in BassTermCounts.Variables(term, kind))
        {
            if (appearances[v] - shift < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
            set[f++] = v;
        }
        if (k == 1) return -Math.Log(p);

        Span<double> z = p <= StackLimit ? stackalloc double[p] : new double[p];
        for (int u = 0; u < p; u++) z[u] = VariableWeight + appearances[u];
        foreach (int v in set) z[v] = VariableWeight + (appearances[v] - shift);
        return Wallenius.LogProbability(z, set);
    }

    private int DrawOrder(Xoshiro256PlusPlus rng, BassStructurePrior.Part part, int[] counts)
    {
        double total = 0.0;
        for (int k = part.MinimumOrder; k <= part.Cap; k++) total += OrderWeight + counts[k];
        double target = rng.NextDouble() * total, cumulative = 0.0;
        for (int k = part.MinimumOrder; k < part.Cap; k++)
        {
            cumulative += OrderWeight + counts[k];
            if (target < cumulative) return k;
        }
        return part.Cap;
    }

    private int[] DrawVariables(Xoshiro256PlusPlus rng, BassStructurePrior.Part part, int[] appearances, int k)
    {
        int p = part.Sizes.Length;
        if (k == 1) return [BassDraws.Below(rng, p)];
        Span<double> z = p <= StackLimit ? stackalloc double[p] : new double[p];
        for (int u = 0; u < p; u++) z[u] = VariableWeight + appearances[u];
        var set = new int[k];
        Wallenius.Sample(rng, z, set);
        return set;
    }

    private static HingeFactor DrawHinge(Xoshiro256PlusPlus rng, BassStructurePrior.Part part, int variable)
        => new(variable, BassDraws.Sign(rng), BassDraws.Below(rng, part.Sizes[variable]));

    // A proper, nonempty subset of n levels from the prior's measure: R's, its size uniform on {1, …, n − 1} and then
    // a uniform subset of that size (a partial shuffle); or the paper's, each level in with probability ½, the empty
    // and full sets refused.
    private CategoricalFactor DrawSubset(Xoshiro256PlusPlus rng, int variable, int n)
    {
        Span<int> pool = n <= StackLimit ? stackalloc int[n] : new int[n];
        int size;
        if (Prior.Subsets == CategoricalSubsetMeasure.SizeUniform)
        {
            size = 1 + BassDraws.Below(rng, n - 1);
            for (int i = 0; i < n; i++) pool[i] = i;
            for (int i = 0; i < size; i++)
            {
                int j = i + BassDraws.Below(rng, n - i);
                (pool[i], pool[j]) = (pool[j], pool[i]);
            }
        }
        else
        {
            do
            {
                size = 0;
                ulong word = 0;
                for (int level = 0; level < n; level++)
                {
                    if ((level & 63) == 0) word = rng.NextUInt64();
                    if ((word & 1) != 0) pool[size++] = level;
                    word >>= 1;
                }
            } while (size == 0 || size == n);
        }
        return new CategoricalFactor(variable, pool[..size].ToArray());
    }

    private void Check(BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (!ReferenceEquals(counts.Prior, Prior)) throw new ArgumentException("The counts are over another prior's term space.", nameof(counts));
    }
}

/// <summary>BASS's exact discrete draws.</summary>
internal static class BassDraws
{
    /// <summary>
    /// A uniform index in [0, count), exactly: Lemire's multiply-shift with rejection. The generator's own
    /// <c>NextInt</c> carries a bias of count/2³², which knot draws over thousands of distinct values would feel in
    /// the proposal law the acceptance ratio assumes.
    /// </summary>
    public static int Below(Xoshiro256PlusPlus rng, int count)
    {
        if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
        ulong bound = (ulong)count;
        ulong high = Math.BigMul(rng.NextUInt64(), bound, out ulong low);
        if (low < bound)
        {
            ulong threshold = unchecked(0UL - bound) % bound;
            while (low < threshold) high = Math.BigMul(rng.NextUInt64(), bound, out low);
        }
        return (int)high;
    }

    /// <summary>A sign, +1 or −1 with probability ½ each, from the word's top bit.</summary>
    public static int Sign(Xoshiro256PlusPlus rng) => (rng.NextUInt64() >> 63) == 0 ? -1 : 1;
}
