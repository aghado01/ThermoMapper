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
/// or subset, is redrawn from the prior; with the local steps on (F-M2's lever 2), with probability ρ stepped locally
/// instead (<see cref="ProposeChange"/>).</item>
/// </list>
/// Each move reports its within-move log ratio, ln q(reverse)/q(forward), which a reversible-jump chain adds to its
/// target ratio (the target holding the structure prior) and its move-selection ratio. A birth's and a death's carry the
/// 1/M of the death's pick; with the prior's ratio they make D02's count factor (M + h₁)/((M + 1)(h₂ + 1)). Immutable,
/// with no scratch state, so one serves every chain of a fit.
/// </summary>
public sealed class BassTermProposal
{
    private const int StackLimit = 256;

    private readonly int[]?[] _windows = new int[]?[3];   // the local steps' w by hinge variable, per hinge part

    /// <summary>Proposals over <paramref name="prior"/>'s term space with the pseudo-counts w₁ and w₂, both positive, and
    /// the local steps <paramref name="steps"/> inside a change (null for D01's redraw alone).</summary>
    public BassTermProposal(BassStructurePrior prior, double orderWeight, double variableWeight, LocalSteps? steps = null)
    {
        ArgumentNullException.ThrowIfNull(prior);
        if (!(orderWeight > 0.0 && double.IsFinite(orderWeight)))
            throw new ArgumentOutOfRangeException(nameof(orderWeight), "w₁ must be positive and finite, or a death could have no reverse birth.");
        if (!(variableWeight > 0.0 && double.IsFinite(variableWeight)))
            throw new ArgumentOutOfRangeException(nameof(variableWeight), "w₂ must be positive and finite, or a death could have no reverse birth.");
        if (steps is not null && !(steps.Probability > 0.0 && steps.Probability <= 1.0 && steps.WindowFraction > 0.0 && steps.WindowFraction <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(steps), "The local steps' ρ and ω lie in (0, 1].");
        Prior = prior;
        OrderWeight = orderWeight;
        VariableWeight = variableWeight;
        Steps = steps;
        if (steps is null) return;
        foreach (BassPart kind in (ReadOnlySpan<BassPart>)[BassPart.Continuous, BassPart.Functional])
        {
            if (prior[kind] is not { } part) continue;
            var w = new int[part.Sizes.Length];
            for (int v = 0; v < w.Length; v++) w[v] = steps.Window(part.Sizes[v]);
            _windows[(int)kind] = w;
        }
    }

    /// <summary>The proposals <paramref name="config"/> sets over <paramref name="prior"/>'s term space.</summary>
    public static BassTermProposal From(BassConfig config, BassStructurePrior prior)
    {
        ArgumentNullException.ThrowIfNull(config);
        return new(prior, config.Structure.OrderWeight, config.Structure.VariableWeight, config.Sampling.LocalSteps);
    }

    /// <summary>The local steps inside a change; null for D01's redraw alone.</summary>
    public LocalSteps? Steps { get; }

    /// <summary>The local steps' window w over hinge variable <paramref name="variable"/> of <paramref name="part"/>
    /// (<see cref="LocalSteps.Window"/>); refused with the steps off.</summary>
    public int Window(BassPart part, int variable)
        => (_windows[(int)part] ?? throw new InvalidOperationException($"No local steps over a {part} part."))[variable];

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
                for (int f = 0; f < k; f++) categorical[f] = part.DrawSubset(rng, variables[f], Prior.Subsets);
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
        return Redraw(rng, term, BassDraws.Below(rng, FactorCount(term)), out logProposalRatio);
    }

    /// <summary>
    /// A change of <paramref name="term"/>, with the local steps when they are on: one of its factors uniformly, then with
    /// probability ρ its local step (<see cref="Stepped"/> for a hinge, <see cref="Toggled"/> for a subset), else its
    /// redraw from the prior as <see cref="DrawChange"/> makes it. With the steps off it draws exactly what
    /// <see cref="DrawChange"/> draws.
    /// </summary>
    public BassTermMove ProposeChange(Xoshiro256PlusPlus rng, BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Prior.Validate(term);
        int pick = BassDraws.Below(rng, FactorCount(term));
        if (Steps is { } steps && rng.NextDouble() < steps.Probability)
        {
            (BassPart kind, int index) = FactorAt(term, pick);
            if (kind == BassPart.Categorical)
                return Toggled(term, pick, BassDraws.Below(rng, Prior[BassPart.Categorical]!.Sizes[term.Categorical[index].Variable]));
            int w = Window(kind, BassStructurePrior.Hinges(term, kind)[index].Variable);
            int size = 1 + BassDraws.Below(rng, w);
            return Stepped(term, pick, BassDraws.Sign(rng) * size);
        }
        BassTerm changed = Redraw(rng, term, pick, out double logProposalRatio);
        return new BassTermMove(BassMoveKind.Redraw, changed, logProposalRatio);
    }

    /// <summary>
    /// Factor <paramref name="pick"/> of <paramref name="term"/> (its factors counted continuous, categorical, then
    /// functional), a hinge, with its knot stepped by <paramref name="delta"/> ≠ 0 and its sign kept: lever 2's hinge
    /// step. A knot off the variable's distinct values is a null move. The step is symmetric, so its log ratio is 0 and
    /// the move's acceptance is Δℓ alone, the knot prior being uniform.
    /// </summary>
    internal BassTermMove Stepped(BassTerm term, int pick, int delta)
    {
        (BassPart kind, int index) = FactorAt(term, pick);
        if (kind == BassPart.Categorical) throw new ArgumentException("A subset toggles; it does not step.", nameof(pick));
        if (delta == 0) throw new ArgumentOutOfRangeException(nameof(delta), "A step moves the knot.");
        HingeFactor[] hinges = [.. BassStructurePrior.Hinges(term, kind)];
        HingeFactor h = hinges[index];
        int knot = h.Knot + delta;
        if (knot < 0 || knot >= Prior[kind]!.Sizes[h.Variable]) return new BassTermMove(BassMoveKind.LocalStep, null, double.NaN);
        hinges[index] = h with { Knot = knot };
        BassTerm stepped = kind == BassPart.Continuous ? new BassTerm(hinges, term.Categorical, term.Functional) : new BassTerm(term.Continuous, term.Categorical, hinges);
        return new BassTermMove(BassMoveKind.LocalStep, stepped, 0.0);
    }

    /// <summary>
    /// Factor <paramref name="pick"/> of <paramref name="term"/>, a subset D, with level <paramref name="level"/>
    /// toggled, D Δ {l}: lever 2's categorical step. An empty or full subset is a null move. The toggle is symmetric, so
    /// its log ratio is 0; the target's prior ratio carries the subset measure's ln π_v(D′) − ln π_v(D).
    /// </summary>
    internal BassTermMove Toggled(BassTerm term, int pick, int level)
    {
        (BassPart kind, int index) = FactorAt(term, pick);
        if (kind != BassPart.Categorical) throw new ArgumentException("A hinge steps; it does not toggle.", nameof(pick));
        CategoricalFactor[] subsets = [.. term.Categorical];
        CategoricalFactor old = subsets[index];
        int levels = Prior[BassPart.Categorical]!.Sizes[old.Variable];
        if ((uint)level >= (uint)levels) throw new ArgumentOutOfRangeException(nameof(level));
        int[] toggled = old.Contains(level) ? [.. old.Levels.Where(l => l != level)] : [.. old.Levels, level];
        if (toggled.Length == 0 || toggled.Length == levels) return new BassTermMove(BassMoveKind.Toggle, null, double.NaN);
        subsets[index] = new CategoricalFactor(old.Variable, toggled);
        return new BassTermMove(BassMoveKind.Toggle, new BassTerm(term.Continuous, subsets, term.Functional), 0.0);
    }

    /// <summary>
    /// A redraw's within-move log ratio, the redrawn factor's prior probability before over after: 0 for a hinge, whose
    /// sign and knot are uniform, and the subset measure's ln π_v(D) − ln π_v(D′) for a subset, which the prior's ratio
    /// cancels. <paramref name="from"/> and <paramref name="to"/> differ in factor <paramref name="pick"/> alone.
    /// </summary>
    internal double LogRedrawRatio(BassTerm from, int pick, BassTerm to)
    {
        (BassPart kind, int index) = FactorAt(from, pick);
        if (kind != BassPart.Categorical) return 0.0;
        double[] log = Prior[BassPart.Categorical]!.LogSubset![from.Categorical[index].Variable];
        return log[from.Categorical[index].Levels.Count] - log[to.Categorical[index].Levels.Count];
    }

    // Factor `pick` redrawn from the prior: DrawChange's draws, in its order.
    private BassTerm Redraw(Xoshiro256PlusPlus rng, BassTerm term, int pick, out double logProposalRatio)
    {
        int continuous = term.Continuous.Count, categorical = term.Categorical.Count;
        logProposalRatio = 0.0;
        if (pick < continuous)
        {
            HingeFactor[] hinges = term.Continuous.ToArray();
            hinges[pick] = DrawHinge(rng, Prior[BassPart.Continuous]!, hinges[pick].Variable);
            return new BassTerm(hinges, term.Categorical, term.Functional);
        }
        if (pick < continuous + categorical)
        {
            BassStructurePrior.Part part = Prior[BassPart.Categorical]!;
            CategoricalFactor[] subsets = term.Categorical.ToArray();
            subsets[pick - continuous] = part.DrawSubset(rng, subsets[pick - continuous].Variable, Prior.Subsets);
            var changed = new BassTerm(term.Continuous, subsets, term.Functional);
            logProposalRatio = LogRedrawRatio(term, pick, changed);
            return changed;
        }
        HingeFactor[] outputs = term.Functional.ToArray();
        outputs[pick - continuous - categorical] = DrawHinge(rng, Prior[BassPart.Functional]!, outputs[pick - continuous - categorical].Variable);
        return new BassTerm(term.Continuous, term.Categorical, outputs);
    }

    private static int FactorCount(BassTerm term) => term.Continuous.Count + term.Categorical.Count + term.Functional.Count;

    // Factor `pick`'s part and its index there, the factors counted continuous, categorical, then functional.
    private static (BassPart Part, int Index) FactorAt(BassTerm term, int pick)
    {
        int continuous = term.Continuous.Count, categorical = term.Categorical.Count;
        if ((uint)pick >= (uint)FactorCount(term)) throw new ArgumentOutOfRangeException(nameof(pick));
        if (pick < continuous) return (BassPart.Continuous, pick);
        if (pick < continuous + categorical) return (BassPart.Categorical, pick - continuous);
        return (BassPart.Functional, pick - continuous - categorical);
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

    private void Check(BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (!ReferenceEquals(counts.Prior, Prior)) throw new ArgumentException("The counts are over another prior's term space.", nameof(counts));
    }
}

/// <summary>A structure move's kind, as the moves' tallies count it: D01's birth, death and redraw, lever 2's local step
/// and toggle, and the fixed-M family's prune, grow and variable change (F-M2's Phase 1).</summary>
public enum BassMoveKind
{
    /// <summary>A term born (D01).</summary>
    Birth,

    /// <summary>A term removed (D01).</summary>
    Death,

    /// <summary>A factor's sign and knot, or subset, redrawn from the prior: D01's change.</summary>
    Redraw,

    /// <summary>A hinge's knot stepped locally, its sign kept (lever 2).</summary>
    LocalStep,

    /// <summary>One level of a subset toggled (lever 2).</summary>
    Toggle,

    /// <summary>A factor pruned from a term's part (the fixed-M family).</summary>
    Prune,

    /// <summary>A factor grown into a term's part (the fixed-M family).</summary>
    Grow,

    /// <summary>A factor's variable exchanged for one its part does not use (the fixed-M family).</summary>
    VariableChange,
}

/// <summary>
/// One term's proposal at fixed M: its kind, the term proposed in its place, and the within-move log ratio
/// ln q(reverse)/q(forward). <see cref="Term"/> is null when nothing is proposed: a null move (a step off the knots, a
/// toggle to an empty or full subset, a prune from an empty part, a grow into a full one), or, with
/// <see cref="LeavesIntercept"/>, a prune of a term's only factor, whose result is the intercept, outside F.
/// </summary>
public readonly record struct BassTermMove(BassMoveKind Kind, BassTerm? Term, double LogProposalRatio, bool LeavesIntercept = false)
{
    /// <summary>Whether the move is a null move: nothing proposed, and nothing refused.</summary>
    public bool IsNull => Term is null && !LeavesIntercept;
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
