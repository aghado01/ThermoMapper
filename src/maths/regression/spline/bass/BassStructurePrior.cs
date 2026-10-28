using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>A term's parts (D01): its continuous hinges, its categorical indicators and its functional hinges.</summary>
public enum BassPart
{
    /// <summary>Hinges on the continuous predictors.</summary>
    Continuous,

    /// <summary>Indicators of subsets of the categorical predictors' levels.</summary>
    Categorical,

    /// <summary>Hinges on the output variables of a direct functional response.</summary>
    Functional,
}

/// <summary>
/// BASS's prior on structures (D01, D02). A collection of M distinct terms has weight
/// Γ(h₁ + M)/(Γ(h₁)(h₂ + 1)^M) · Πₘ π(termₘ), restricted globally to F: the Poisson count prior with its M! orderings,
/// λ ~ Gamma(h₁, rate h₂) integrated out, times the term measure, with no per-order or per-λ normalizer.
/// <para>The term measure, for each part a fit has: the order uniform (on {1, …, cap} when the fit has one part; on
/// {0, …, cap} for each part when it has several, the term whose parts are all empty being the intercept, outside F);
/// given the order K ≥ 1, the variable set uniform over the C(p, K) subsets; each hinge's sign uniform and its knot
/// uniform over its variable's distinct values; each categorical subset by the declared measure, R's (its size uniform
/// on {1, …, n − 1}, then uniform among the subsets of that size) or the paper's (uniform over the 2ⁿ − 2 proper
/// subsets). One part and K ≤ K_max make this Francom et al.'s c_K/K_max.</para>
/// <para><see cref="LogPrior"/> scores the term measure and <see cref="DrawTerm"/> draws from it. Immutable, so one
/// serves every chain of a fit.</para>
/// </summary>
public sealed class BassStructurePrior
{
    private readonly Part?[] _parts = new Part?[3];

    /// <summary>
    /// The prior on <paramref name="inputs"/>' predictors and output variables, with each part's interaction cap (0
    /// where the fit has no variables of that kind), the categorical subset measure and the count prior.
    /// </summary>
    public BassStructurePrior(BassInputs inputs, int maxInteraction, int maxCategoricalInteraction, int maxFunctionalInteraction,
                              CategoricalSubsetMeasure subsets, BassCountPrior count)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        ArgumentNullException.ThrowIfNull(count);
        if (!(count.Shape > 0.0 && double.IsFinite(count.Shape))) throw new ArgumentOutOfRangeException(nameof(count), "h₁ must be positive and finite.");
        if (!(count.Rate > 0.0 && double.IsFinite(count.Rate))) throw new ArgumentOutOfRangeException(nameof(count), "h₂ must be positive and finite.");
        if (subsets is not (CategoricalSubsetMeasure.SizeUniform or CategoricalSubsetMeasure.Uniform))
            throw new ArgumentOutOfRangeException(nameof(subsets));
        Count = count;
        Subsets = subsets;

        CaseInputs cases = inputs.Cases;
        int present = (maxInteraction > 0 ? 1 : 0) + (maxCategoricalInteraction > 0 ? 1 : 0) + (maxFunctionalInteraction > 0 ? 1 : 0);
        if (present == 0) throw new ArgumentException("A structure prior needs a part with a positive cap.");
        int minimum = present > 1 ? 0 : 1;

        var knots = new int[cases.ContinuousCount];
        for (int v = 0; v < knots.Length; v++) knots[v] = inputs.ContinuousKnots(v).Length;
        _parts[(int)BassPart.Continuous] = Part.Hinges(BassPart.Continuous, knots, maxInteraction, minimum, nameof(maxInteraction));

        var levels = new int[cases.CategoricalCount];
        for (int v = 0; v < levels.Length; v++) levels[v] = inputs.LevelCount(v);
        _parts[(int)BassPart.Categorical] = Part.Indicators(levels, maxCategoricalInteraction, minimum, subsets);

        var outputKnots = new int[inputs.Outputs?.Count ?? 0];
        for (int v = 0; v < outputKnots.Length; v++) outputKnots[v] = inputs.FunctionalKnots(v).Length;
        _parts[(int)BassPart.Functional] = Part.Hinges(BassPart.Functional, outputKnots, maxFunctionalInteraction, minimum, nameof(maxFunctionalInteraction));
    }

    /// <summary>The prior <paramref name="config"/> sets, with the caps <paramref name="resolved"/> holds.</summary>
    public static BassStructurePrior From(BassConfig config, BassResolved resolved, BassInputs inputs)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(resolved);
        return new(inputs, resolved.MaxInteraction, resolved.MaxCategoricalInteraction, resolved.MaxFunctionalInteraction,
                   config.Structure.CategoricalSubsets, config.Count);
    }

    /// <summary>The count prior's λ hyperparameters, h₁ and h₂ (D02).</summary>
    public BassCountPrior Count { get; }

    /// <summary>The categorical subset measure.</summary>
    public CategoricalSubsetMeasure Subsets { get; }

    /// <summary>Whether the fit has <paramref name="part"/>: a positive cap.</summary>
    public bool Has(BassPart part) => _parts[(int)part] is not null;

    /// <summary><paramref name="part"/>'s interaction cap, 0 where the fit lacks it.</summary>
    public int Cap(BassPart part) => _parts[(int)part]?.Cap ?? 0;

    /// <summary><paramref name="part"/>'s smallest order: 0 when the fit has several parts, else 1.</summary>
    public int MinimumOrder(BassPart part) => _parts[(int)part]?.MinimumOrder ?? 0;

    /// <summary><paramref name="part"/>'s variables, 0 where the fit lacks it.</summary>
    public int Variables(BassPart part) => _parts[(int)part]?.Sizes.Length ?? 0;

    internal Part? this[BassPart part] => _parts[(int)part];

    /// <summary>ln π(term), the term measure's log probability.</summary>
    public double LogPrior(BassTerm term)
    {
        Validate(term);
        double lp = 0.0;
        foreach (Part? p in _parts)
        {
            if (p is null) continue;
            int k = Order(term, p.Kind);
            lp += p.LogOrder;
            if (k == 0) continue;
            lp += p.LogSets[k];
            lp += p.LogFactors(term);
        }
        return lp;
    }

    /// <summary>
    /// A term drawn from the term measure, the law <see cref="LogPrior"/> scores: each part's order uniform, its
    /// variables a uniform subset of that size, each hinge's sign and knot uniform, each subset from the declared
    /// measure. Null is the intercept, the draw whose parts all came out empty; it has mass only when the fit has several
    /// parts, and it lies outside F.
    /// </summary>
    public BassTerm? DrawTerm(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Span<int> orders = stackalloc int[3];
        bool empty = true;
        for (int index = 0; index < 3; index++)
        {
            if (_parts[index] is not { } p) continue;
            orders[index] = p.MinimumOrder + BassDraws.Below(rng, p.OrderCount);
            empty &= orders[index] == 0;
        }
        if (empty) return null;

        HingeFactor[]? continuous = null, functional = null;
        CategoricalFactor[]? categorical = null;
        for (int index = 0; index < 3; index++)
        {
            int k = orders[index];
            if (k == 0) continue;
            Part p = _parts[index]!;
            int[] variables = UniformSubset(rng, p.Sizes.Length, k);
            if (p.Kind == BassPart.Categorical)
            {
                categorical = new CategoricalFactor[k];
                for (int f = 0; f < k; f++) categorical[f] = p.DrawSubset(rng, variables[f], Subsets);
                continue;
            }
            var hinges = new HingeFactor[k];
            for (int f = 0; f < k; f++)
                hinges[f] = new HingeFactor(variables[f], BassDraws.Sign(rng), BassDraws.Below(rng, p.Sizes[variables[f]]));
            if (p.Kind == BassPart.Continuous) continuous = hinges;
            else functional = hinges;
        }
        return new BassTerm(continuous, categorical, functional);
    }

    // k of the p variables, uniformly among the C(p, k) subsets: a partial Fisher–Yates shuffle.
    private static int[] UniformSubset(Xoshiro256PlusPlus rng, int p, int k)
    {
        var pool = new int[p];
        for (int v = 0; v < p; v++) pool[v] = v;
        for (int i = 0; i < k; i++)
        {
            int j = i + BassDraws.Below(rng, p - i);
            (pool[i], pool[j]) = (pool[j], pool[i]);
        }
        return pool[..k];
    }

    /// <summary>
    /// ln of the count prior's weight on a collection of <paramref name="terms"/> distinct terms with its orderings,
    /// ln Γ(h₁ + M) − ln Γ(h₁) − M ln(h₂ + 1). A structure's log prior is this plus its terms' <see cref="LogPrior"/>,
    /// on F.
    /// </summary>
    public double LogCollectionWeight(int terms)
    {
        if (terms < 0) throw new ArgumentOutOfRangeException(nameof(terms));
        return SpecialFunctions.LogGamma(Count.Shape + terms) - SpecialFunctions.LogGamma(Count.Shape)
               - terms * Math.Log(Count.Rate + 1.0);
    }

    /// <summary>
    /// ln p(M + 1)/p(M) = ln[(M + h₁)/((M + 1)(h₂ + 1))], the count law's ratio with λ integrated out: a birth's count
    /// factor, and the prior-cancelling schedule's input (D01, D02).
    /// </summary>
    public double LogCountRatio(int terms)
    {
        if (terms < 0) throw new ArgumentOutOfRangeException(nameof(terms));
        return Math.Log(terms + Count.Shape) - Math.Log(terms + 1.0) - Math.Log(Count.Rate + 1.0);
    }

    /// <summary>
    /// Refuses a term outside the term space: a factor in a part the fit lacks, an order above its part's cap, a
    /// variable, knot or level beyond the inputs', or a subset that is not proper.
    /// </summary>
    public void Validate(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        for (int index = 0; index < 3; index++)
        {
            var kind = (BassPart)index;
            int k = Order(term, kind);
            Part? p = _parts[(int)kind];
            if (p is null)
            {
                if (k > 0) throw new ArgumentException($"The fit has no {Name(kind)} part.", nameof(term));
                continue;
            }
            if (k > p.Cap) throw new ArgumentException($"The term's {Name(kind)} order {k} exceeds the cap {p.Cap}.", nameof(term));
            if (kind == BassPart.Categorical)
            {
                foreach (CategoricalFactor c in term.Categorical)
                {
                    if (c.Variable >= p.Sizes.Length) throw new ArgumentException($"Categorical variable {c.Variable} does not exist.", nameof(term));
                    int levels = p.Sizes[c.Variable];
                    if (c.Levels[^1] >= levels) throw new ArgumentException($"Categorical variable {c.Variable} has {levels} levels.", nameof(term));
                    if (c.Levels.Count >= levels) throw new ArgumentException($"A subset of categorical variable {c.Variable} must leave a level out.", nameof(term));
                }
                continue;
            }
            foreach (HingeFactor h in Hinges(term, kind))
            {
                if (h.Variable >= p.Sizes.Length) throw new ArgumentException($"The {Name(kind)} variable {h.Variable} does not exist.", nameof(term));
                if (h.Knot >= p.Sizes[h.Variable]) throw new ArgumentException($"Knot {h.Knot} of {Name(kind)} variable {h.Variable} does not exist.", nameof(term));
            }
        }
    }

    internal static int Order(BassTerm term, BassPart part) => part switch
    {
        BassPart.Continuous => term.Continuous.Count,
        BassPart.Categorical => term.Categorical.Count,
        _ => term.Functional.Count,
    };

    internal static IReadOnlyList<HingeFactor> Hinges(BassTerm term, BassPart part)
        => part == BassPart.Continuous ? term.Continuous : term.Functional;

    private static string Name(BassPart part) => part switch
    {
        BassPart.Continuous => "continuous",
        BassPart.Categorical => "categorical",
        _ => "functional",
    };

    /// <summary>One part's measure: its cap and smallest order, and each variable's size (knots or levels) and
    /// factor log probabilities.</summary>
    internal sealed class Part
    {
        private Part(BassPart kind, int[] sizes, int cap, int minimum)
        {
            Kind = kind;
            Sizes = sizes;
            Cap = cap;
            MinimumOrder = minimum;
            LogOrder = -Math.Log(cap - minimum + 1);
            LogSets = new double[cap + 1];
            for (int k = 1; k <= cap; k++) LogSets[k] = -LogChoose(sizes.Length, k);
        }

        public BassPart Kind { get; }

        /// <summary>Each variable's knots (hinge parts) or training levels (the categorical part).</summary>
        public int[] Sizes { get; }

        public int Cap { get; }

        public int MinimumOrder { get; }

        /// <summary>The orders this part's I(K) ranges over, Cap − MinimumOrder + 1.</summary>
        public int OrderCount => Cap - MinimumOrder + 1;

        /// <summary>ln of the order's prior probability, the same for every order.</summary>
        public double LogOrder { get; }

        /// <summary>−ln C(p, K) by K, the variable set's log probability.</summary>
        public double[] LogSets { get; }

        /// <summary>A hinge's log probability by variable, −ln 2 − ln n_v.</summary>
        public double[]? LogHinge { get; private init; }

        /// <summary>A subset's log probability by variable and size.</summary>
        public double[][]? LogSubset { get; private init; }

        public static Part? Hinges(BassPart kind, int[] knots, int cap, int minimum, string capName)
        {
            if (cap < 0 || cap > knots.Length)
                throw new ArgumentOutOfRangeException(capName, $"The cap is at least 0 and at most the {knots.Length} variables.");
            if (cap == 0) return null;
            var log = new double[knots.Length];
            for (int v = 0; v < log.Length; v++) log[v] = -Math.Log(2.0) - Math.Log(knots[v]);
            return new Part(kind, knots, cap, minimum) { LogHinge = log };
        }

        public static Part? Indicators(int[] levels, int cap, int minimum, CategoricalSubsetMeasure measure)
        {
            if (cap < 0 || cap > levels.Length)
                throw new ArgumentOutOfRangeException(nameof(cap), $"The cap is at least 0 and at most the {levels.Length} variables.");
            if (cap == 0) return null;
            var log = new double[levels.Length][];
            for (int v = 0; v < log.Length; v++)
            {
                int n = levels[v];
                if (n < 2) throw new ArgumentException($"Categorical variable {v} has fewer than two levels.", nameof(levels));
                log[v] = new double[n];
                // The paper's uniform over the 2ⁿ − 2 proper subsets: ln(2ⁿ − 2) = n ln 2 + ln(1 − 2^(1−n)).
                double uniform = -(n * Math.Log(2.0) + Math.Log(1.0 - Math.Pow(2.0, 1 - n)));
                for (int s = 1; s < n; s++)
                    log[v][s] = measure == CategoricalSubsetMeasure.Uniform ? uniform : -Math.Log(n - 1) - LogChoose(n, s);
            }
            return new Part(BassPart.Categorical, levels, cap, minimum) { LogSubset = log };
        }

        /// <summary>The term's factors in this part: signs and knots, or subsets.</summary>
        public double LogFactors(BassTerm term)
        {
            double lp = 0.0;
            if (LogSubset is { } subset)
                foreach (CategoricalFactor c in term.Categorical) lp += subset[c.Variable][c.Levels.Count];
            else
                foreach (HingeFactor h in BassStructurePrior.Hinges(term, Kind)) lp += LogHinge![h.Variable];
            return lp;
        }

        /// <summary>
        /// A proper, nonempty subset of <paramref name="variable"/>'s levels from <paramref name="measure"/>: R's, its
        /// size uniform on {1, …, n − 1} and then a uniform subset of that size (a partial shuffle); or the paper's,
        /// each level in with probability ½, the empty and full sets refused. The categorical part only.
        /// </summary>
        public CategoricalFactor DrawSubset(Xoshiro256PlusPlus rng, int variable, CategoricalSubsetMeasure measure)
        {
            int n = Sizes[variable];
            Span<int> pool = n <= 256 ? stackalloc int[n] : new int[n];
            int size;
            if (measure == CategoricalSubsetMeasure.SizeUniform)
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

        private static double LogChoose(int n, int k)
            => SpecialFunctions.LogGamma(n + 1.0) - SpecialFunctions.LogGamma(k + 1.0) - SpecialFunctions.LogGamma(n - k + 1.0);
    }
}

/// <summary>
/// A structure's counts, which D01's birth proposal reads: for each part, how many terms take each order, and how
/// many terms use each variable. Kept beside the structure, a term added or removed with it.
/// </summary>
public sealed class BassTermCounts
{
    private readonly int[][] _orders = new int[3][];
    private readonly int[][] _appearances = new int[3][];

    /// <summary>The empty structure's counts under <paramref name="prior"/>.</summary>
    public BassTermCounts(BassStructurePrior prior)
    {
        ArgumentNullException.ThrowIfNull(prior);
        Prior = prior;
        for (int p = 0; p < 3; p++)
        {
            _orders[p] = new int[prior.Cap((BassPart)p) + 1];
            _appearances[p] = new int[prior.Variables((BassPart)p)];
        }
    }

    private BassTermCounts(BassTermCounts source)
    {
        Prior = source.Prior;
        Terms = source.Terms;
        for (int p = 0; p < 3; p++)
        {
            _orders[p] = (int[])source._orders[p].Clone();
            _appearances[p] = (int[])source._appearances[p].Clone();
        }
    }

    /// <summary>The prior whose term space these count.</summary>
    public BassStructurePrior Prior { get; }

    /// <summary>The structure's terms, M.</summary>
    public int Terms { get; private set; }

    /// <summary>The terms whose <paramref name="part"/> has order <paramref name="order"/>.</summary>
    public int Order(BassPart part, int order) => _orders[(int)part][order];

    /// <summary>The terms whose <paramref name="part"/> uses <paramref name="variable"/>.</summary>
    public int Appearances(BassPart part, int variable) => _appearances[(int)part][variable];

    internal int[] Orders(BassPart part) => _orders[(int)part];

    internal int[] AppearanceCounts(BassPart part) => _appearances[(int)part];

    /// <summary>Counts <paramref name="term"/> in.</summary>
    public void Add(BassTerm term)
    {
        Prior.Validate(term);
        Shift(term, +1);
        Terms++;
    }

    /// <summary>Counts <paramref name="term"/> out; refused when the counts hold no such term.</summary>
    public void Remove(BassTerm term)
    {
        Prior.Validate(term);
        if (Terms == 0 || !Holds(term)) throw new ArgumentException("The counts hold no such term.", nameof(term));
        Shift(term, -1);
        Terms--;
    }

    /// <summary>An independent copy.</summary>
    public BassTermCounts Clone() => new(this);

    private bool Holds(BassTerm term)
    {
        for (int p = 0; p < 3; p++)
        {
            var part = (BassPart)p;
            if (!Prior.Has(part)) continue;
            if (_orders[p][BassStructurePrior.Order(term, part)] == 0) return false;
            foreach (int v in Variables(term, part))
                if (_appearances[p][v] == 0) return false;
        }
        return true;
    }

    private void Shift(BassTerm term, int by)
    {
        for (int p = 0; p < 3; p++)
        {
            var part = (BassPart)p;
            if (!Prior.Has(part)) continue;
            _orders[p][BassStructurePrior.Order(term, part)] += by;
            foreach (int v in Variables(term, part)) _appearances[p][v] += by;
        }
    }

    internal static IEnumerable<int> Variables(BassTerm term, BassPart part)
    {
        if (part == BassPart.Categorical)
            foreach (CategoricalFactor c in term.Categorical) yield return c.Variable;
        else
            foreach (HingeFactor h in BassStructurePrior.Hinges(term, part)) yield return h.Variable;
    }
}
