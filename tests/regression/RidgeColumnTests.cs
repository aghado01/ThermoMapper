using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 3 gates on a ridge block's columns (the BPPR brief; <see cref="RidgeColumn"/>), against the paper's eq. (spline_expansion)
/// and eq. (eta_-0) by hand:
/// <list type="bullet">
/// <item>the knots, the type-7 quantiles of the projected training cases above t₀, against R's arrangement by hand, and the ties that leave a
/// column undefined;</item>
/// <item>the columns at training cases and past the range on every side, against the paper's formulas in decimal, so the cubes' difference
/// keeps its digits where the engine reads the linear form;</item>
/// <item>a seeded sweep of knots, columns and far points, values and derivatives against the decimal hand formulas;</item>
/// <item>derivatives against central differences through the basis;</item>
/// <item>through a fit: prediction and the gradient readout against the hand computation, the moment service's refusal, and E8's sampled C
/// rank one along θ's raw-unit image;</item>
/// <item>the column's and the term's refusals, and value equality.</item>
/// </list>
/// </summary>
public sealed class RidgeColumnTests(ITestOutputHelper output)
{
    private static readonly ProjectionDirection Dense = new([0, 1, 2], [0.6, -0.48, 0.64]);
    private static readonly ProjectionDirection Pair = new([2, 0], [-0.6, 0.8]);
    private static readonly ProjectionDirection Single = new([1], [1.0]);

    private static BassConfig Config() => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 3, MaxBasis = 6, Directional = new DirectionalBasis() },
        Sampling = new BassSampling { Chains = 1, MasterSeed = 36, BurnIn = 20, Samples = 2 },
    };

    private static DirectionalTermSpace Space(BassData data) => BassModel.Create(Config(), data).Prior.Directional!;

    // θᵀz at raw points, z by the compensated training mean and sample deviation.
    private static double[] Projections(BassData data, double[][] raw, ProjectionDirection d)
    {
        (double Mean, double Deviation)[] m = DirectionalHingeTests.Moments(data);
        return [.. Enumerable.Range(0, raw[0].Length).Select(i =>
            Enumerable.Range(0, d.Count).Sum(k => d.Weights[k] * (raw[d.Variables[k]][i] - m[d.Variables[k]].Mean) / m[d.Variables[k]].Deviation))];
    }

    // R's quantile(x, probs, type = 7) on sorted x, probs from seq(0, 1, length.out = K + 1).
    private static double[] HandKnots(double[] projections, double t0, int columns)
    {
        double[] above = [.. projections.Where(p => p > t0).Order()];
        var knots = new double[columns + 2];
        knots[0] = t0;
        for (int l = 0; l <= columns; l++)
        {
            double probability = l == 0 ? 0.0 : l == columns ? 1.0 : l * (1.0 / columns);
            double index = 1.0 + (above.Length - 1) * probability;
            int lo = (int)Math.Floor(index);
            double h = index - lo;
            knots[l + 1] = h > 0.0 && above[lo] != above[lo - 1] ? (1.0 - h) * above[lo - 1] + h * above[lo] : above[lo - 1];
        }
        return knots;
    }

    // An exact rational, for the hand formulas: every double is a binary fraction, so differences, cubes and quotients of the knots and u
    // carry no rounding until the result is read back as a double.
    private readonly record struct Exact(BigInteger Num, BigInteger Den)
    {
        public static Exact Of(double d)
        {
            long bits = BitConverter.DoubleToInt64Bits(d);
            bool negative = bits < 0;
            int exponent = (int)((bits >> 52) & 0x7FF);
            long mantissa = bits & 0xFFFFFFFFFFFFFL;
            if (exponent == 0) exponent++; else mantissa |= 1L << 52;
            exponent -= 1075;
            var m = new BigInteger(negative ? -mantissa : mantissa);
            return exponent >= 0 ? new(m << exponent, BigInteger.One) : new(m, BigInteger.One << -exponent);
        }
        public static Exact operator +(Exact a, Exact b) => new(a.Num * b.Den + b.Num * a.Den, a.Den * b.Den);
        public static Exact operator -(Exact a, Exact b) => new(a.Num * b.Den - b.Num * a.Den, a.Den * b.Den);
        public static Exact operator *(Exact a, Exact b) => new(a.Num * b.Num, a.Den * b.Den);
        public static Exact operator /(Exact a, Exact b) => b.Num.Sign < 0 ? new(-a.Num * b.Den, a.Den * -b.Num) : new(a.Num * b.Den, a.Den * b.Num);
        public Exact Positive => Num.Sign > 0 ? this : new(BigInteger.Zero, BigInteger.One);
        public Exact Abs => Num.Sign < 0 ? new(-Num, Den) : this;
        public double ToDouble()
        {
            if (Num.IsZero) return 0.0;
            BigInteger n = BigInteger.Abs(Num);
            int shift = Math.Max(0, (int)(Den.GetBitLength() - n.GetBitLength()) + 64);
            double q = (double)((n << shift) / Den);
            return (Num.Sign < 0 ? -1.0 : 1.0) * Math.ScaleB(q, -shift);
        }
    }

    // The paper's columns and their derivatives by exact arithmetic: b₁ = (u − t₀)₊, b_ℓ = d_{ℓ−1} − d_K, d_j = ((u − t_j)₊³ − (u − t_{K+1})₊³)/(t_{K+1} − t_j);
    // `Scale` is |d_{ℓ−1}| + |d_K| and `SlopeScale` |d′_{ℓ−1}| + |d′_K|, the sizes of what each difference cancels.
    internal static (double Value, double Derivative, double Scale, double SlopeScale) Paper(double[] knots, int column, double u)
    {
        Exact x = Exact.Of(u), t0 = Exact.Of(knots[0]);
        if (column == 1) return ((x - t0).Positive.ToDouble(), u > knots[0] ? 1.0 : u == knots[0] ? 0.5 : 0.0, Math.Abs(u), 1.0);
        int k = knots.Length - 2;
        Exact top = Exact.Of(knots[k + 1]), three = Exact.Of(3.0);
        (Exact D, Exact Slope) Divided(int j)
        {
            Exact t = Exact.Of(knots[j]), a = (x - t).Positive, b = (x - top).Positive;
            return ((a * a * a - b * b * b) / (top - t), three * (a * a - b * b) / (top - t));
        }
        (Exact first, Exact firstSlope) = Divided(column - 1);
        (Exact last, Exact lastSlope) = Divided(k);
        return ((first - last).ToDouble(), (firstSlope - lastSlope).ToDouble(), (first.Abs + last.Abs).ToDouble(), (firstSlope.Abs + lastSlope.Abs).ToDouble());
    }

    private static BassTerm Column(ProjectionDirection d, double[] knots, int column) => new(ridge: new RidgeColumn(d, knots, column));

    // ---- the knots -----------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Knots_AreTheTypeSevenQuantilesAboveTheFirstKnot_ByHand()
    {
        BassData data = DirectionalHingeTests.Data();
        DirectionalTermSpace space = Space(data);
        double[][] raw = DirectionalHingeTests.Raw(data);
        double worst = 0.0;
        foreach (ProjectionDirection d in new[] { Dense, Pair, Single })
        {
            double[] p = Projections(data, raw, d);
            double[] sorted = [.. p.Order()];
            // t₀ below every case, then between two cases near the bottom and halfway up: off every case, so the strict cut reads the same
            // cases whether the projections are the engine's or the hand's.
            foreach (double at in new[] { -1.0, 0.0, 0.5 })
            foreach (int columns in new[] { 1, 2, 4, 6 })
            {
                int lower = (int)(at * (sorted.Length - 2));
                double t0 = at < 0.0 ? sorted[0] - 1.0 : 0.5 * (sorted[lower] + sorted[lower + 1]);
                double[] hand = HandKnots(p, t0, columns);
                double[] engine = space.RidgeKnots(new DirectionalHinge(d, t0), columns)!;
                Assert.Equal(columns + 2, engine.Length);
                Assert.Equal(t0, engine[0]);
                for (int j = 1; j < hand.Length; j++) worst = Math.Max(worst, Math.Abs(engine[j] - hand[j]) / Math.Max(1.0, Math.Abs(hand[j])));
            }
        }
        output.WriteLine($"worst |knot − hand| / max(1, |hand|) = {worst:E2}");
        Assert.True(worst <= 1e-14, $"{worst:E2}");

        // Ties leave a column undefined: on x₃ tied at 0 on 115 of 120 cases, −e₃'s projections tie at their top, so above any t₀ below the
        // tie the upper quantiles sit in it and t_K = t_{K+1} for K ≥ 2; the hinge alone (K = 1) is defined; above the top there are no cases.
        BassData tied = DirectionalChainTests.Data(tied: true);
        DirectionalTermSpace tiedSpace = Space(tied);
        var down = new ProjectionDirection([3], [-1.0]);
        double[] q = [.. Projections(tied, DirectionalHingeTests.Raw(tied), down).Order()];
        double below = 0.5 * (q[0] + q[^1]);
        Assert.Null(tiedSpace.RidgeKnots(new DirectionalHinge(down, below), 4));
        Assert.Null(tiedSpace.RidgeKnots(new DirectionalHinge(down, q[0] - 1.0), 4));   // the five low cases reach only the low quantiles
        Assert.NotNull(tiedSpace.RidgeKnots(new DirectionalHinge(down, below), 1));
        Assert.Null(tiedSpace.RidgeKnots(new DirectionalHinge(down, q[^1]), 4));
    }

    // ---- the columns ---------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Columns_MatchThePapersExpansion_ByHand()
    {
        BassData data = DirectionalHingeTests.Data();
        (BassInputTransform transform, _, BassBasis basis) = DirectionalHingeTests.Setup(data, 1);
        DirectionalTermSpace space = Space(data);
        double[][] raw = DirectionalHingeTests.Raw(data);
        double[][] points = DirectionalHingeTests.Points(data);
        CaseInputs at = DirectionalHingeTests.Scaled(transform, points);
        double worst = 0.0;
        foreach ((ProjectionDirection d, double fraction) in new[] { (Dense, 0.2), (Pair, 0.5), (Single, 0.0) })
        {
            double[] sorted = [.. Projections(data, raw, d).Order()];
            double t0 = fraction == 0.0 ? sorted[0] - 0.5 : sorted[(int)(fraction * (sorted.Length - 1))];
            foreach (int columns in new[] { 1, 4, 6 })
            {
                double[] knots = space.RidgeKnots(new DirectionalHinge(d, t0), columns)!;
                double[] u = Projections(data, points, d);
                for (int l = 1; l <= columns; l++)
                {
                    var column = new double[at.Cases];
                    basis.CaseColumn(Column(d, knots, l), at, column);
                    for (int i = 0; i < u.Length; i++)
                    {
                        (double value, _, double scale, _) = Paper(knots, l, u[i]);
                        worst = Math.Max(worst, Math.Abs(column[i] - value) / Math.Max(1.0, scale));
                    }
                }
            }
        }
        output.WriteLine($"worst |column − paper| / max(1, |d_ℓ₋₁| + |d_K|) = {worst:E2} (training cases and points past the range)");
        Assert.True(worst <= 1e-13, $"{worst:E2}");
    }

    [Fact]
    public void SeededSweep_ValuesAndDerivatives_MatchThePapersFormulas()
    {
        var rng = new Xoshiro256PlusPlus(3641);
        double worstValue = 0.0, worstSlope = 0.0;
        int cases = 0;
        for (int c = 0; c < 400; c++)
        {
            int columns = 1 + (int)(rng.NextDouble() * 6);
            double scale = Math.Pow(10.0, -1.0 + 3.0 * rng.NextDouble());   // knot spreads 0.1 to 100
            double t0 = scale * (2.0 * rng.NextDouble() - 1.0);
            var rest = new double[columns + 1];
            for (int j = 0; j < rest.Length; j++) rest[j] = t0 + scale * (1e-3 + rng.NextDouble());
            Array.Sort(rest);
            double[] knots = [t0, .. rest];
            var direction = new ProjectionDirection([0], [1.0]);
            for (int l = 1; l <= columns; l++)
            {
                var r = new RidgeColumn(direction, knots, l);
                for (int k = 0; k < 12; k++)
                {
                    // Inside the knots, at them, and far past them on both sides.
                    double u = k switch
                    {
                        < 6 => knots[0] + (knots[^1] - knots[0]) * (1.2 * rng.NextDouble() - 0.1),
                        6 => knots[1 + (int)(rng.NextDouble() * (columns + 1))],
                        7 => knots[0],
                        _ => (k % 2 == 0 ? 1.0 : -1.0) * scale * Math.Pow(10.0, 3.0 * rng.NextDouble()) + knots[^1],
                    };
                    (double value, double slope, double size, double slopeSize) = Paper(knots, l, u);
                    worstValue = Math.Max(worstValue, Math.Abs(r.Value(u) - value) / Math.Max(1.0, size));
                    worstSlope = Math.Max(worstSlope, Math.Abs(r.Derivative(u) - slope) / Math.Max(1.0, slopeSize));
                    cases++;
                }
            }
        }
        output.WriteLine($"{cases} points: worst value {worstValue:E2}, worst derivative {worstSlope:E2}");
        Assert.True(worstValue <= 1e-14, $"value {worstValue:E2}");
        Assert.True(worstSlope <= 1e-14, $"derivative {worstSlope:E2}");
    }

    [Fact]
    public void Derivatives_MatchCentralDifferences_ThroughTheBasis()
    {
        BassData data = DirectionalHingeTests.Data();
        (BassInputTransform transform, BassInputs inputs, BassBasis basis) = DirectionalHingeTests.Setup(data, 1);
        DirectionalTermSpace space = Space(data);
        double[] sorted = [.. Projections(data, DirectionalHingeTests.Raw(data), Dense).Order()];
        double[] knots = space.RidgeKnots(new DirectionalHinge(Dense, 0.5 * (sorted[20] + sorted[21])), 4)!;   // b₁'s kink off every case
        CaseInputs at = inputs.Cases;
        const double h = 1e-6;
        double worst = 0.0, largest = 0.0;
        for (int l = 1; l <= 4; l++)
        {
            BassTerm term = Column(Dense, knots, l);
            for (int v = 0; v < 3; v++)
            {
                var derivative = new double[at.Cases];
                basis.CaseDerivative(term, v, at, derivative);
                double[][] up = [.. Enumerable.Range(0, 3).Select(j => at.Continuous(j).ToArray())], down = [.. up.Select(c => (double[])c.Clone())];
                for (int i = 0; i < at.Cases; i++) { up[v][i] += h; down[v][i] -= h; }
                var above = new double[at.Cases];
                var belowColumn = new double[at.Cases];
                basis.CaseColumn(term, new CaseInputs(up, []), above);
                basis.CaseColumn(term, new CaseInputs(down, []), belowColumn);
                for (int i = 0; i < at.Cases; i++)
                {
                    largest = Math.Max(largest, Math.Abs(derivative[i]));
                    worst = Math.Max(worst, Math.Abs((above[i] - belowColumn[i]) / (2.0 * h) - derivative[i]));
                }
            }
        }
        output.WriteLine($"worst |difference − derivative| / (max |derivative| + 1) = {worst / (largest + 1.0):E2}");
        Assert.True(worst / (largest + 1.0) <= 1e-7, $"{worst / (largest + 1.0):E2}");
    }

    // ---- through a fit -------------------------------------------------------------------------------------------------------------

    [Fact]
    public void ThroughAFit_PredictionAndTheGradientReadTheBlock_TheMomentServiceRefuses_AndE8Samples()
    {
        BassData data = DirectionalHingeTests.Data();
        BassFit run = Bass.Fit(new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 6 },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 36, BurnIn = 20, Samples = 2 },
        }, data);
        DirectionalTermSpace space = Space(data);
        double[] sorted = [.. Projections(data, DirectionalHingeTests.Raw(data), Dense).Order()];
        double[] knots = space.RidgeKnots(new DirectionalHinge(Dense, sorted[15]), 4)!;
        var table = new BassTermTable();
        int[] ids = [.. Enumerable.Range(1, 4).Select(l => table.Intern(Column(Dense, knots, l)))];
        double[] beta = [0.4, 1.5, -0.7, 2.2, 0.9];
        var store = new DrawStore<BassDraw>();
        store.Add(new BassDraw(new BassStructure(ids), 1.0, 1.0, beta, 1.0));
        BassFit fit = run with { Terms = table, Draws = [store], DrawWeights = null };

        double[][] points = DirectionalHingeTests.Points(data);
        PredictionCases cases = PredictionCases.Of([.. points.Select((column, j) => (BassPredictor)new ContinuousPredictor(data.Predictors[j].Name, column))]);
        double[,,] latent = fit.Predict(data, cases).Latent;
        BassGradient gradient = fit.Gradient(data, cases);
        (double Mean, double Deviation)[] moments = DirectionalHingeTests.Moments(data);
        double[] u = Projections(data, points, Dense);
        double worstValue = 0.0, worstSlope = 0.0;
        for (int i = 0; i < u.Length; i++)
        {
            double f = beta[0], outer = 0.0;
            for (int l = 1; l <= 4; l++)
            {
                (double value, double slope, _, _) = Paper(knots, l, u[i]);
                f += beta[l] * value;
                outer += beta[l] * slope;
            }
            worstValue = Math.Max(worstValue, Math.Abs(latent[0, i, 0] - f) / Math.Max(1.0, Math.Abs(f)));
            for (int j = 0; j < 3; j++)
            {
                double slope = outer * Dense.WeightOn(j) / moments[j].Deviation;
                worstSlope = Math.Max(worstSlope, Math.Abs(gradient.Of(data.Predictors[j].Name)[0, i, 0] - slope) / Math.Max(1.0, Math.Abs(slope)));
            }
        }
        output.WriteLine($"prediction {worstValue:E2}, gradient {worstSlope:E2} (relative to max(1, |·|))");
        Assert.True(worstValue <= 1e-12, $"prediction {worstValue:E2}");
        Assert.True(worstSlope <= 1e-12, $"gradient {worstSlope:E2}");

        // Sobol's exact route refuses up front, naming the term, and its automatic route takes the sampled route (the sampled-sensitivity
        // brief, decision 1).
        ArgumentException refusal = Assert.Throws<ArgumentException>(() => fit.Sensitivity(data, new SensitivitySpec { Route = SensitivityRoute.Exact }));
        Assert.Contains("not a product of one-variable factors", refusal.Message);
        Assert.Contains(table[ids[0]].ToString(), refusal.Message);
        Assert.Equal(SensitivityRoute.Sampled, fit.Sensitivity(data, new SensitivitySpec { Points = 64 }).Route);
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { Route = ActiveSubspaceRoute.Exact }));
        Assert.Contains("ridge block's column", Assert.ThrowsAny<Exception>(() => BassSensitivity.FromCoefficients(fit, data)).Message);

        // One block: every gradient is a multiple of θ's raw-unit image v_j = θ_j/sd_j, so the sampled C is rank one along v.
        ActiveSubspaceAnalysis analysis = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Metric = InputMetric.Raw, Points = 256 });
        Assert.Equal(ActiveSubspaceRoute.Sampled, analysis.Route);
        Assert.Null(analysis.Poincare);
        double[,] meanC = analysis.Posterior.MeanC;
        double[] v = [.. Enumerable.Range(0, 3).Select(j => Dense.WeightOn(j) / moments[j].Deviation)];
        double norm = Math.Sqrt(v.Sum(x => x * x));
        double largest = 0.0, residual = 0.0, along = 0.0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) along += v[a] * meanC[a, b] * v[b] / (norm * norm);
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
            {
                largest = Math.Max(largest, Math.Abs(meanC[a, b]));
                residual = Math.Max(residual, Math.Abs(meanC[a, b] - along * v[a] * v[b] / (norm * norm)));
            }
        output.WriteLine($"C off its rank-one image along θ's raw-unit image: {residual / largest:E2} of its largest entry");
        Assert.True(residual / largest <= 1e-13, $"{residual / largest:E2}");
    }

    // ---- the refusals --------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TheColumnAndTheTerm_RefuseWhatIsNotOne_AndCompareByValue()
    {
        double[] knots = [-1.0, -0.5, 0.0, 0.5, 1.0, 2.0];   // K = 4
        Assert.Throws<ArgumentException>(() => new RidgeColumn(Dense, [0.0, 1.0], 1));                    // K = 0
        Assert.Throws<ArgumentOutOfRangeException>(() => new RidgeColumn(Dense, knots, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new RidgeColumn(Dense, knots, 5));
        Assert.Throws<ArgumentException>(() => new RidgeColumn(Dense, [-0.5, -0.5, 0.0, 1.0], 1));       // t₀ not below t₁
        Assert.Throws<ArgumentException>(() => new RidgeColumn(Dense, [-1.0, 0.5, 0.0, 1.0], 1));        // decreasing
        Assert.Throws<ArgumentException>(() => new RidgeColumn(Dense, [-1.0, 0.0, 1.0, 1.0], 2));        // t_K = t_{K+1}, K = 2
        Assert.Throws<ArgumentException>(() => new RidgeColumn(Dense, [-1.0, 0.0, double.NaN, 1.0], 1));
        Assert.NotNull(new RidgeColumn(Dense, [-1.0, 0.0, 0.0], 1));                                      // K = 1 reads t₀ alone

        var ridge = new RidgeColumn(Dense, knots, 2);
        Assert.Throws<ArgumentException>(() => new BassTerm([new HingeFactor(0, 1, 3)], ridge: ridge));
        Assert.Throws<ArgumentException>(() => new BassTerm(directional: new DirectionalHinge(Dense, 0.0), ridge: ridge));
        Assert.Throws<ArgumentException>(() => new BassTerm(legendre: [new LegendreFactor(0, 2)], ridge: ridge));
        Assert.Equal(Column(Dense, knots, 2), Column(new ProjectionDirection([2, 1, 0], [0.64, -0.48, 0.6]), [.. knots], 2));
        Assert.Equal(Column(Dense, knots, 2).GetHashCode(), Column(Dense, [.. knots], 2).GetHashCode());
        Assert.NotEqual(Column(Dense, knots, 2), Column(Dense, knots, 3));
        Assert.NotEqual(Column(Dense, knots, 2), Column(Dense, [-1.0, -0.5, 0.0, 0.5, 1.0, 2.0000000000000004], 2));
        Assert.NotEqual(Column(Dense, knots, 1), new BassTerm(directional: new DirectionalHinge(Dense, -1.0)));

        // A structure prior refuses a column term: a hinge fit's, and a directional fit's, whose structure terms are (θ, t₀).
        BassData data = DirectionalHingeTests.Data();
        BassTerm column = Column(Dense, knots, 1);
        Assert.Throws<ArgumentException>(() => BassModel.Create(Config() with { Structure = Config().Structure with { Directional = null } }, data).Prior.Validate(column));
        Assert.Throws<ArgumentException>(() => BassModel.Create(Config(), data).Prior.Validate(column));
    }
}
