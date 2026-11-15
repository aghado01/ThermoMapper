using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 2 gates on directions in the structure chain (the BPPR brief; decisions 4 and 5):
/// <list type="bullet">
/// <item>the knot's range against R's type-7 quantile by hand, on the selection buffer, the sort and an index with no fraction, and a range
/// without width, whose terms lie outside the support;</item>
/// <item>the term measure's log density against its closed form, its refusals, and its draws against their laws;</item>
/// <item>the birth's law: its (a, 𝒥) mass summing to one, its draws against it, Wallenius' term by quadrature apart from the engine's;</item>
/// <item>each move's ratio, birth, death and change, against a direct evaluation: the dense marginal and the hand prior and proposal;</item>
/// <item>the collapsed law given frozen directions against the dense marginal;</item>
/// <item>Geweke's test with the prior alone, since directions are continuous: the chain at b = 0 against the prior sampler's draws;</item>
/// <item>the config's and the resolution's refusals, and the readouts that refuse a directional term with their reasons.</item>
/// </list>
/// The hand columns are <see cref="DirectionalHingeTests"/>'s, from the raw inputs and a compensated standardization.
/// </summary>
public sealed class DirectionalChainTests(ITestOutputHelper output)
{
    private static readonly ProjectionDirection[] Directions =
    [
        new([1], [1.0]),
        new([2], [-1.0]),
        new([0, 1], [0.6, 0.8]),
        new([1, 3], [-0.8, 0.6]),
        new([0, 2, 3], [0.48, -0.6, 0.64]),
    ];

    // Four inputs on [0, 1], each declaring that domain, on golden-ratio lattices; y the single index exp(0.6x₀ + 0.8x₁) and noise. With
    // `tied`, x₃ is 0 on all but five cases, so +e₃'s projections tie at their least value up to the knot's upper end.
    internal static BassData Data(int n = 120, bool tied = false, int seed = 5)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] alpha = [0.6180339887498949, 0.4142135623730950, 0.7320508075688772, 0.2360679774997897];
        var x = new double[4][];
        for (int j = 0; j < 4; j++)
        {
            x[j] = new double[n];
            for (int i = 0; i < n; i++) x[j][i] = (0.5 + i * alpha[j]) % 1.0;
        }
        if (tied)
            for (int i = 0; i < n; i++) x[3][i] = i % 24 == 0 ? 0.2 * (i / 24 + 1) : 0.0;
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = Math.Exp(0.6 * x[0][i] + 0.8 * x[1][i]) + 0.05 * Normal.SampleStandard(rng);
        return new BassData([.. Enumerable.Range(0, 4).Select(j => new ContinuousPredictor($"x{j}", x[j], new InputDomain(0.0, 1.0)))], new ScalarResponse(y));
    }

    // A directional basis under proper priors that read no response, so b = 0 exists: the zero centre, σ² ~ IG(3, 2), τ held.
    internal static BassConfig Config(int maxBasis = 8, int degree = 1, double? activeCases = null) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 3, MaxBasis = maxBasis, Degree = degree, Directional = new DirectionalBasis { ActiveCases = activeCases } },
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.InverseGamma(3, 2),
        Coefficients = BassCoefficientPrior.G(TauPrior.Held()),
        Sampling = new BassSampling { Chains = 1, MasterSeed = 11, BurnIn = 100, Samples = 100 },
    };

    private static BassTerm Term(ProjectionDirection direction, double knot) => new(directional: new DirectionalHinge(direction, knot));

    // A knot at fraction `at` of its range from L.
    private static BassTerm Within(DirectionalTermSpace space, ProjectionDirection direction, double at)
    {
        KnotRange r = space.KnotRangeOf(direction)!.Value;
        return Term(direction, r.Lower + at * r.Width);
    }

    // θᵀz at the training cases, z from the raw values by the compensated mean and sample deviation.
    private static double[] Projections(BassData data, ProjectionDirection d)
    {
        double[][] raw = DirectionalHingeTests.Raw(data);
        (double Mean, double Deviation)[] m = DirectionalHingeTests.Moments(data);
        return [.. Enumerable.Range(0, raw[0].Length).Select(i =>
            Enumerable.Range(0, d.Count).Sum(k => d.Weights[k] * (raw[d.Variables[k]][i] - m[d.Variables[k]].Mean) / m[d.Variables[k]].Deviation))];
    }

    // ---- the knot's range ----------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData(120, null)]    // q = 0.9: 13 order statistics, the insertion buffer
    [InlineData(120, 80.0)]    // 81 order statistics, the sort
    [InlineData(11, 1.1)]      // 1 + (n − 1)q = 10 exactly: no fraction
    public void KnotRange_IsTheTypeSevenQuantileRange_ByHand(int n, double? active)
    {
        BassData data = Data(n);
        DirectionalTermSpace space = BassModel.Create(Config(activeCases: active), data).Prior.Directional!;
        double q = 1.0 - (active ?? Math.Min(20.0, 0.1 * n)) / n;
        Assert.Equal(q, space.KnotLevel);
        Assert.Equal(2.0 / 3.0, space.LocalProbability);
        double worst = 0.0;
        foreach (ProjectionDirection d in Directions)
        {
            double[] p = [.. Projections(data, d).Order()];
            double index = 1.0 + (n - 1) * q;
            int lo = (int)Math.Floor(index);
            double h = index - lo;
            if (active == 1.1) Assert.Equal(0.0, h);
            double upper = h > 0.0 && p[lo] != p[lo - 1] ? (1.0 - h) * p[lo - 1] + h * p[lo] : p[lo - 1];
            double lower = upper - (upper - p[0]) / (2.0 / 3.0);
            KnotRange range = space.KnotRangeOf(d)!.Value;
            double scale = Math.Max(1.0, Math.Max(Math.Abs(upper), Math.Abs(lower)));
            worst = Math.Max(worst, Math.Max(Math.Abs(range.Upper - upper), Math.Abs(range.Lower - lower)) / scale);
            Assert.Equal(range, space.KnotRangeOf(new ProjectionDirection(d.Variables, d.Weights))!.Value);   // an equal direction, computed afresh
        }
        output.WriteLine($"n {n}, q {q:R}: worst |Δ| over the ranges {worst:E2}");
        Assert.True(worst <= 1e-13, $"{worst:E2}");
    }

    [Fact]
    public void ARangeWithoutWidth_HasNoKnot_AndItsTermsLieOutsideTheSupport()
    {
        BassData data = Data(tied: true);
        BassModel model = BassModel.Create(Config(), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        var up = new ProjectionDirection([3], [1.0]);
        var down = new ProjectionDirection([3], [-1.0]);
        Assert.Null(space.KnotRangeOf(up));
        Assert.NotNull(space.KnotRangeOf(down));
        Assert.Equal(double.NegativeInfinity, model.Prior.LogPrior(Term(up, 0.0)));

        // From −e₃ the fair coin proposes +e₃ about half the time: a null move, nothing proposed; −e₃ again with a fresh knot otherwise.
        var counts = new BassTermCounts(model.Prior);
        BassTerm from = Within(space, down, 0.5);
        counts.Add(from);
        var rng = new Xoshiro256PlusPlus(3611);
        int nulls = 0;
        for (int k = 0; k < 2000; k++)
        {
            BassTermMove move = model.Kernel.ProposeChange(rng, from, counts);
            Assert.Equal(BassMoveKind.DirectionChange, move.Kind);
            if (move.Term is null) nulls++;
            else Assert.Equal(down, move.Term.Directional!.Direction);
        }
        output.WriteLine($"null moves {nulls} of 2000");
        Assert.InRange(nulls, 900, 1100);
    }

    // ---- the term measure ----------------------------------------------------------------------------------------------------------

    [Fact]
    public void LogPrior_IsTheTermMeasure_AndTheSpaceRefusesWhatIsNotInIt()
    {
        BassModel model = BassModel.Create(Config(), Data());
        BassStructurePrior prior = model.Prior;
        DirectionalTermSpace space = prior.Directional!;
        Assert.True(prior.Has(BassPart.Directional));
        Assert.False(prior.Has(BassPart.Continuous));
        Assert.Equal(3, prior.Cap(BassPart.Directional));
        Assert.Equal(4, prior.Variables(BassPart.Directional));
        double[] logSphere = [0.0, -Math.Log(2.0), -Math.Log(2.0 * Math.PI), -Math.Log(4.0 * Math.PI)];   // −ln|S^(a−1)|
        int[] choose = [1, 4, 6, 4];
        foreach (ProjectionDirection d in Directions)
        {
            KnotRange r = space.KnotRangeOf(d)!.Value;
            int a = d.Count;
            double expected = -Math.Log(3.0) - Math.Log(choose[a]) + logSphere[a] - Math.Log(r.Upper - r.Lower);
            Assert.Equal(expected, prior.LogPrior(Within(space, d, 0.3)), 12);
            Assert.Equal(logSphere[a], space.LogDirection(a), 14);
            Assert.Equal(double.NegativeInfinity, prior.LogPrior(Term(d, r.Upper + 1e-6 * r.Width)));
            Assert.Equal(double.NegativeInfinity, prior.LogPrior(Term(d, r.Lower - 1e-6 * r.Width)));
        }
        Assert.Throws<ArgumentException>(() => prior.Validate(new BassTerm([new HingeFactor(0, 1, 3)])));
        Assert.Throws<ArgumentException>(() => prior.Validate(Term(new ProjectionDirection([0, 1, 2, 3], [0.5, 0.5, 0.5, 0.5]), 0.0)));
        Assert.Throws<ArgumentException>(() => prior.Validate(Term(new ProjectionDirection([4], [1.0]), 0.0)));

        // A hinge fit's prior refuses a directional term, rather than scoring it as an empty one.
        BassConfig hinges = Config() with { Structure = Config().Structure with { Directional = null } };
        ArgumentException refusal = Assert.Throws<ArgumentException>(() => BassModel.Create(hinges, Data()).Prior.Validate(Term(Directions[0], 0.0)));
        Assert.Contains("no directional part", refusal.Message);
    }

    [Fact]
    public void DrawTerm_DrawsTheTermMeasure()
    {
        BassModel model = BassModel.Create(Config(), Data());
        BassStructurePrior prior = model.Prior;
        DirectionalTermSpace space = prior.Directional!;
        var rng = new Xoshiro256PlusPlus(3602);
        const int draws = 60_000;
        var byActive = new int[4];
        var included = new int[4, 4];
        var squares = new List<double>[4];
        var positive = new int[4];
        var positions = new double[draws];
        for (int a = 1; a <= 3; a++) squares[a] = [];
        for (int k = 0; k < draws; k++)
        {
            DirectionalHinge h = prior.DrawTerm(rng)!.Directional!;
            int a = h.Direction.Count;
            byActive[a]++;
            foreach (int v in h.Direction.Variables) included[a, v]++;
            squares[a].Add(h.Direction.Weights[0] * h.Direction.Weights[0]);
            if (h.Direction.Weights[0] > 0.0) positive[a]++;
            KnotRange r = space.KnotRangeOf(h.Direction)!.Value;
            positions[k] = (r.Upper - h.Knot) / r.Width;
        }
        var z = new List<(string, double)>();
        for (int a = 1; a <= 3; a++)
        {
            z.Add(($"P(a = {a})", (byActive[a] - draws / 3.0) / Math.Sqrt(draws / 3.0 * (2.0 / 3.0))));
            for (int v = 0; v < 4; v++)
            {
                double p = a / 4.0;
                z.Add(($"P(x{v} | a = {a})", (included[a, v] - byActive[a] * p) / Math.Sqrt(byActive[a] * p * (1.0 - p))));
            }
            if (a == 1) Assert.All(squares[a], s => Assert.Equal(1.0, s));   // S⁰: θ₁ = ±1
            else
            {
                double mean = squares[a].Average(), sd = Math.Sqrt(squares[a].Sum(s => (s - mean) * (s - mean)) / (squares[a].Count - 1));
                z.Add(($"E[θ₁² | a = {a}]", (mean - 1.0 / a) / (sd / Math.Sqrt(squares[a].Count))));
            }
            z.Add(($"P(θ₁ > 0 | a = {a})", (positive[a] - 0.5 * byActive[a]) / Math.Sqrt(0.25 * byActive[a])));
        }
        Array.Sort(positions);
        double d = positions.Select((r, i) => Math.Max(Math.Abs(r - (double)i / draws), Math.Abs(r - (i + 1.0) / draws))).Max();
        foreach ((string name, double value) in z) output.WriteLine($"{name,-22} z {value,6:F2}");
        output.WriteLine($"the knot's position on its range: √n·D {Math.Sqrt(draws) * d:F2}");
        Assert.True(z.Max(e => Math.Abs(e.Item2)) < 4.0, "a draw's law departs from the measure's");
        Assert.True(Math.Sqrt(draws) * d < 1.8, "the knot is not uniform on its range");
    }

    // ---- the birth's law -----------------------------------------------------------------------------------------------------------

    // Wallenius' probability that weighted sampling without replacement draws the set `set`, by quadrature of its integral
    // ∫₀¹ Π_{j∈set}(1 − t^{ω_j/D}) dt, D the weight left out, under t = s^m so that the integrand is smooth: apart from the engine's.
    private static double WalleniusByQuadrature(double[] weights, int[] set)
    {
        double left = weights.Where((_, u) => !set.Contains(u)).Sum();
        double[] c = [.. set.Select(j => weights[j] / left)];
        const int m = 16, intervals = 20_000;
        double F(double s)
        {
            double value = m * Math.Pow(s, m - 1);
            foreach (double cj in c) value *= 1.0 - Math.Pow(s, m * cj);
            return value;
        }
        double h = 1.0 / intervals, sum = F(0.0) + F(1.0);
        for (int i = 1; i < intervals; i++) sum += (i % 2 == 1 ? 4.0 : 2.0) * F(i * h);
        return sum * h / 3.0;
    }

    private static IEnumerable<int[]> Subsets(int p, int a)
    {
        if (a == 0) { yield return []; yield break; }
        for (int first = 0; first < p; first++)
            foreach (int[] rest in Subsets(p, a - 1))
                if (rest.Length == 0 || rest[0] > first) yield return [first, .. rest];
    }

    // The order and variable laws' log probability by hand: q(a) = (w₁ + c_a)/(A·w₁ + M), q(𝒥 | a) = 1/p at a = 1, Wallenius' by quadrature
    // otherwise, with z_u = w₂ + (terms using u) less `shift` copies of the term's own.
    private static double HandLogOrderAndVariables(BassTermCounts counts, int[] set, int shift, double w1 = 5.0, double w2 = 5.0)
    {
        int a = set.Length, terms = counts.Terms - shift;
        double lq = Math.Log(w1 + counts.Order(BassPart.Directional, a) - shift) - Math.Log(3.0 * w1 + terms);
        if (a == 1) return lq - Math.Log(4.0);
        double[] z = [.. Enumerable.Range(0, 4).Select(u => w2 + counts.Appearances(BassPart.Directional, u) - (shift > 0 && set.Contains(u) ? shift : 0))];
        return lq + Math.Log(WalleniusByQuadrature(z, set));
    }

    [Fact]
    public void TheBirth_DrawsItsLaw_AndLogBirthIsItsExactProbability()
    {
        BassModel model = BassModel.Create(Config(), Data());
        var kernel = Assert.IsType<DirectionalTermProposal>(model.Kernel);
        DirectionalTermSpace space = model.Prior.Directional!;
        var counts = new BassTermCounts(model.Prior);
        foreach (BassTerm t in new[] { Within(space, Directions[2], 0.4), Within(space, Directions[0], 0.6), Within(space, Directions[4], 0.2) }) counts.Add(t);

        // Every (a, 𝒥): the engine's order and variable mass against the hand law, and the masses summing to one.
        var cells = new List<(int[] Set, double P)>();
        double worst = 0.0, total = 0.0;
        for (int a = 1; a <= 3; a++)
            foreach (int[] set in Subsets(4, a))
            {
                var d = new ProjectionDirection(set, Enumerable.Repeat(1.0 / Math.Sqrt(a), a));
                BassTerm term = Within(space, d, 0.5);
                double engine = kernel.LogBirth(term, counts) - space.LogDirection(a) - space.LogKnot(term.Directional!);
                double hand = HandLogOrderAndVariables(counts, set, 0);
                worst = Math.Max(worst, Math.Abs(engine - hand));
                total += Math.Exp(engine);
                cells.Add((set, Math.Exp(engine)));
            }
        output.WriteLine($"the order and variable law: worst |Δ ln q| against the hand law {worst:E2}; total mass {total:R}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.Equal(1.0, total, 12);

        // The draws' (a, 𝒥) against that mass, and their θ and knot from the priors.
        var rng = new Xoshiro256PlusPlus(3603);
        const int draws = 50_000;
        var seen = new int[cells.Count];
        for (int k = 0; k < draws; k++)
        {
            ProjectionDirection d = kernel.DrawBirth(rng, counts)!.Directional!.Direction;
            seen[cells.FindIndex(c => c.Set.SequenceEqual(d.Variables))]++;
        }
        double chi2 = cells.Select((c, i) => Math.Pow(seen[i] - draws * c.P, 2) / (draws * c.P)).Sum();
        int df = cells.Count - 1;
        output.WriteLine($"births' (a, 𝒥) against the law over {cells.Count} cells: χ² {chi2:F1} on {df}");
        Assert.True(chi2 <= df + 6 * Math.Sqrt(2.0 * df), $"χ² = {chi2:F1} on {df}");
    }

    // ---- the dense marginal and each move's ratio ----------------------------------------------------------------------------------

    // ln p(y | B) under the zero centre, held τ and σ² ~ IG(a, b), by dense algebra: y | σ² ~ N(0, σ²(I + B(τBᵀB)⁻¹Bᵀ)) with B = [1, columns],
    // σ² integrated: ln Γ(a + n/2) − ln Γ(a) + a ln b − (n/2) ln 2π − ½ ln|Σ| − (a + n/2) ln(b + yᵀΣ⁻¹y/2).
    internal static double DenseLogMarginal(double[] y, IReadOnlyList<double[]> columns, double tau, double shape, double scale)
    {
        int n = y.Length, k = columns.Count + 1;
        double[] Column(int j) => j == 0 ? [.. Enumerable.Repeat(1.0, n)] : columns[j - 1];
        var g = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++)
            {
                double[] a = Column(i), b = Column(j);
                double s = 0.0;
                for (int r = 0; r < n; r++) s += a[r] * b[r];
                g[i, j] = s;
            }
        double[,] gInverse = Inverse(g);
        var sigma = new double[n, n];
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
            {
                double s = 0.0;
                for (int i = 0; i < k; i++)
                    for (int j = 0; j < k; j++) s += Column(i)[r] * gInverse[i, j] * Column(j)[c];
                sigma[r, c] = (r == c ? 1.0 : 0.0) + s / tau;
            }
        double[,] l = Cholesky(sigma);
        double logDet = 0.0;
        for (int i = 0; i < n; i++) logDet += 2.0 * Math.Log(l[i, i]);
        var w = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = y[i];
            for (int j = 0; j < i; j++) s -= l[i, j] * w[j];
            w[i] = s / l[i, i];
        }
        double quad = w.Sum(v => v * v);
        return SpecialFunctions.LogGamma(shape + 0.5 * n) - SpecialFunctions.LogGamma(shape) + shape * Math.Log(scale) - 0.5 * n * Math.Log(2.0 * Math.PI)
             - 0.5 * logDet - (shape + 0.5 * n) * Math.Log(scale + 0.5 * quad);
    }

    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = a[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        return l;
    }

    private static double[,] Inverse(double[,] a)
    {
        int n = a.GetLength(0);
        double[,] l = Cholesky(a);
        var inverse = new double[n, n];
        for (int c = 0; c < n; c++)
        {
            var x = new double[n];
            for (int i = 0; i < n; i++)
            {
                double s = i == c ? 1.0 : 0.0;
                for (int j = 0; j < i; j++) s -= l[i, j] * x[j];
                x[i] = s / l[i, i];
            }
            for (int i = n - 1; i >= 0; i--)
            {
                double s = x[i];
                for (int j = i + 1; j < n; j++) s -= l[j, i] * x[j];
                x[i] = s / l[i, i];
            }
            for (int i = 0; i < n; i++) inverse[i, c] = x[i];
        }
        return inverse;
    }

    private static double Dense(BassModel model, BassData data, IEnumerable<BassTerm> terms)
    {
        double[][] raw = DirectionalHingeTests.Raw(data);
        int degree = model.Config.Structure.Degree;
        double[][] columns = [.. terms.Select(t => DirectionalHingeTests.Hand(data, raw, t.Directional!.Direction, t.Directional.Knot, degree, HingeNormalization.Box).Values)];
        double[] y = ((ScalarResponse)data.Response).Values.ToArray();
        return DenseLogMarginal(y, columns, model.Target.HeldTau!.Value, 3.0, 2.0);
    }

    private static BassEvaluation Evaluate(BassModel model, BassTerm[] terms)
    {
        var chain = new BassChain(model, terms, Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        return model.Target.Evaluate(chain.Current);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void TheCollapsedLaw_GivenFrozenDirections_IsTheDenseMarginal(int degree)
    {
        BassData data = Data();
        BassModel model = BassModel.Create(Config(degree: degree), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        BassTerm[] terms = [Within(space, Directions[2], 0.7), Within(space, Directions[1], 0.4), Within(space, Directions[4], 0.9)];
        double worst = 0.0;
        for (int m = 0; m <= terms.Length; m++)
        {
            BassEvaluation e = Evaluate(model, terms[..m]);
            double dense = Dense(model, data, terms[..m]);
            worst = Math.Max(worst, Math.Abs(e.LogMarginal - dense) / Math.Max(1.0, Math.Abs(dense)));
            output.WriteLine($"M = {m}: engine {e.LogMarginal:R}, dense {dense:R}");
        }
        output.WriteLine($"degree {degree}: worst relative |Δ ln p(y | B)| {worst:E2}");
        Assert.True(worst <= 1e-11, $"{worst:E2}");
    }

    [Fact]
    public void EachMovesRatio_MatchesADirectEvaluation()
    {
        BassData data = Data();
        BassModel model = BassModel.Create(Config(), data);
        BassStructurePrior prior = model.Prior;
        DirectionalTermSpace space = prior.Directional!;
        IBassProposalKernel kernel = model.Kernel;
        BassCountPrior count = model.Config.Count;
        double Target(BassTerm[] s)
        {
            BassEvaluation e = Evaluate(model, s);
            return e.LogPrior + e.LogMarginal;
        }

        BassTerm[] s = [Within(space, Directions[2], 0.7), Within(space, Directions[1], 0.4)];
        var counts = new BassTermCounts(prior);
        foreach (BassTerm t in s) counts.Add(t);
        double worst = 0.0;
        foreach (ProjectionDirection d in Directions)
        {
            // The birth of `born` into S, M = 2: the engine's target ratio and within-move ratio against the dense marginal's ratio, the count
            // law's (M + h₁)/((M + 1)(h₂ + 1)), and π(a)π(𝒥|a)/(q(a)q(𝒥|a)) by hand; θ's and the knot's densities cancel.
            BassTerm born = Within(space, d, 0.55);
            BassTerm[] grown = [.. s, born];
            int a = d.Count, m = s.Length;
            double engine = Target(grown) - Target(s) + kernel.LogBirthRatio(born, counts);
            double direct = Dense(model, data, grown) - Dense(model, data, s)
                          + Math.Log((m + count.Shape) / ((m + 1.0) * (count.Rate + 1.0)))
                          - Math.Log(3.0) - Math.Log(new[] { 1, 4, 6, 4 }[a]) - HandLogOrderAndVariables(counts, [.. d.Variables], 0);
            // The reverse death from S + born: its ratio is the birth's negated.
            BassTermCounts after = counts.Clone();
            after.Add(born);
            double death = Target(s) - Target(grown) + kernel.LogDeathRatio(born, after);
            worst = Math.Max(worst, Math.Max(Math.Abs(engine - direct), Math.Abs(death + direct)));
            output.WriteLine($"birth of {born}: engine {engine:F12}, direct {direct:F12}; death {death:F12}");
        }

        // The change: a and 𝒥 kept, θ by the power spherical law (a fair sign at a = 1), the knot from its prior given θ*; the engine's
        // ratio with its within-move ratio against the dense marginal's ratio alone.
        var rng = new Xoshiro256PlusPlus(3604);
        foreach (int index in new[] { 0, 1 })
            for (int k = 0; k < 4; k++)
            {
                BassTermMove move = kernel.ProposeChange(rng, s[index], counts);
                BassTerm moved = move.Term!;
                Assert.Equal(s[index].Directional!.Direction.Variables, moved.Directional!.Direction.Variables);
                BassTerm[] changed = [.. s];
                changed[index] = moved;
                double engine = Target(changed) - Target(s) + move.LogProposalRatio;
                double direct = Dense(model, data, changed) - Dense(model, data, s);
                worst = Math.Max(worst, Math.Abs(engine - direct));
                output.WriteLine($"change of {s[index]} to {moved}: engine {engine:F12}, direct {direct:F12}");
            }
        output.WriteLine($"worst |Δ| over the moves {worst:E2}");
        Assert.True(worst <= 1e-9, $"{worst:E2}");
    }

    // ---- Geweke's test with the prior alone ----------------------------------------------------------------------------------------

    private static readonly string[] Statistics =
        ["M", "a = 1", "a = 2", "a = 3", "uses x0", "uses x1", "uses x2", "uses x3", "θ mass x0", "θ mass x1", "θ mass x2", "θ mass x3",
         "θ₁ > 0", "knot position", "knot position²"];

    private static double[] Functions(BassState state, DirectionalTermSpace space)
    {
        var f = new double[Statistics.Length];
        f[0] = state.Count;
        foreach (BassTerm t in state.Terms)
        {
            DirectionalHinge h = t.Directional!;
            f[h.Direction.Count]++;
            for (int k = 0; k < h.Direction.Count; k++)
            {
                f[4 + h.Direction.Variables[k]]++;
                f[8 + h.Direction.Variables[k]] += h.Direction.Weights[k] * h.Direction.Weights[k];
            }
            if (h.Direction.Weights[0] > 0.0) f[12]++;
            KnotRange r = space.KnotRangeOf(h.Direction)!.Value;
            double position = (r.Upper - h.Knot) / r.Width;
            f[13] += position;
            f[14] += position * position;
        }
        return f;
    }

    private static double BatchError(double[] x, int batches)
    {
        int size = x.Length / batches;
        double[] means = [.. Enumerable.Range(0, batches).Select(b => x.Skip(b * size).Take(size).Average())];
        double grand = means.Average();
        return Math.Sqrt(means.Sum(m => (m - grand) * (m - grand)) / (batches - 1) / batches);
    }

    [Theory]
    [InlineData(3605)]
    [InlineData(3615)]
    public void Geweke_ThePriorChain_LeavesThePriorInvariant(int seed)
    {
        BassData data = Data(80);
        BassModel model = BassModel.Create(Config(maxBasis: 4), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        var rng = new Xoshiro256PlusPlus(seed);
        var sampler = new BassPriorSampler(model);
        const int draws = 60_000, burn = 5_000, steps = 300_000, batches = 50;
        var direct = new double[draws][];
        for (int k = 0; k < draws; k++) direct[k] = Functions(sampler.DrawState(rng).State, space);

        var chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(seed + 1), new Xoshiro256PlusPlus(seed + 2));
        chain.Moves.SetTemperature(new Temperature(0.0, TemperingPath.Likelihood));
        for (int t = 0; t < burn; t++) chain.Step();
        var successive = new double[steps][];
        for (int t = 0; t < steps; t++) successive[t] = Functions(chain.Step(), space);

        double worst = 0.0;
        for (int g = 0; g < Statistics.Length; g++)
        {
            double[] a = [.. direct.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(BatchError(b, batches), 2));
            double z = (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{Statistics[g],-16} prior {meanA,9:F5}  chain {meanB,9:F5}  z {z,6:F2}");
        }
        output.WriteLine($"worst |z| {worst:F2}; moves {string.Join(", ", chain.Structure.Kinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}; "
                         + $"prior sampler acceptance {sampler.AcceptanceRate:F3}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- what a directional basis refuses ------------------------------------------------------------------------------------------

    [Fact]
    public void TheConfigAndTheResolution_RefuseWhatADirectionalBasisDoesNotRead()
    {
        BassConfig good = Config();
        good.Validate();
        BassStructureSpec s = good.Structure;
        BassConfig With(BassSampling sampling) => good with { Sampling = sampling };
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Legendre = new LegendreBasis() } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Directional = new DirectionalBasis { LocalProbability = 0.0 } } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Directional = new DirectionalBasis { LocalProbability = 1.5 } } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Directional = new DirectionalBasis { ActiveCases = 0.0 } } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Directional = new DirectionalBasis { Concentration = -1.0 } } }).Validate());
        Assert.Throws<ArgumentException>(() => With(good.Sampling with { LocalSteps = new LocalSteps() }).Validate());
        Assert.Throws<ArgumentException>(() => With(good.Sampling with { FactorMoves = new FactorMoves() }).Validate());
        Assert.Throws<ArgumentException>(() => With(good.Sampling with { Births = new BirthMixture { History = new HistoryBirth() } }).Validate());
        Assert.Throws<ArgumentException>(() => With(good.Sampling with { Births = new BirthMixture { Residual = new ResidualBirth() } }).Validate());

        // Data a directional basis does not serve, ActiveCases that leave no case below U, KHAOS's modified g-prior.
        BassData categorical = new([new ContinuousPredictor("x", [0, 0.5, 1, 0.2]), new CategoricalPredictor("z", [0, 1, 0, 1], ["a", "b"])], new ScalarResponse([1, 2, 3, 2.5]));
        Assert.Contains("decision 7", Assert.Throws<ArgumentException>(() => BassModel.Create(good, categorical)).Message);
        Assert.Throws<ArgumentException>(() => BassModel.Create(good with { Structure = s with { Directional = new DirectionalBasis { ActiveCases = 120.0 } } }, Data()));
        Assert.Throws<ArgumentException>(() => BassModel.Create(good with { Coefficients = BassCoefficientPrior.ModifiedG(TauPrior.Held()) }, Data()));

        BassModel model = BassModel.Create(good, Data());
        Assert.IsType<DirectionalTermProposal>(model.Kernel);
        Assert.Contains("E36", Assert.Throws<InvalidOperationException>(() => model.Proposal).Message);
        Assert.Equal(0, model.Resolved.MaxCategoricalInteraction);
        Assert.Contains(model.Resolved.AutoPicks, p => p.Setting == "Directional");

        // The config round-trips through JSON with its basis, and a hinge config's form leaves it out.
        Assert.Equal(good.Structure.Directional, BassConfig.FromJson(good.ToJson()).Structure.Directional);
        Assert.DoesNotContain("Directional", (good with { Structure = s with { Directional = null } }).ToJson());
    }
}
