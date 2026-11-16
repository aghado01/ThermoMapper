using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 1 gates on the directional hinge (θᵀz − t)₊^q, BASS's fourth case-side term type (the BPPR brief):
/// <list type="bullet">
/// <item>the direction and the term refuse what is not one, and a directional hinge takes categorical factors beside it alone (decision 6;
/// slice 5's gate, decision 7);</item>
/// <item>columns against a hand computation from the raw inputs, its own mean and deviation (decision 3's metric), and c at the box's
/// maximizing corner; and a seeded sweep of directions, knots, degrees and normalizations, far points included;</item>
/// <item>derivatives against central differences, and at the knot at degree one the mean of the one-sided limits;</item>
/// <item>D06 admits θ and −θ and several knots on one θ, and refuses a near-duplicate;</item>
/// <item>through a fit: prediction and the gradient readout against the hand computation, the moment service's refusal, and E8's
/// automatic route taking the sampled route, where C is rank one along θ's raw-unit image.</item>
/// </list>
/// Existing fits are bit for bit with the type off: the pinned tier and the whole assembly carry that gate.
/// </summary>
public sealed class DirectionalHingeTests(ITestOutputHelper output)
{
    private static readonly ProjectionDirection Dense = new([0, 1, 2], [0.6, -0.48, 0.64]);
    private static readonly ProjectionDirection Single = new([1], [1.0]);
    private static readonly ProjectionDirection Pair = new([2, 0], [-0.6, 0.8]);

    // Three continuous inputs on unequal ranges: x0 on [0, 1], x1 on [−2, 5] skewed toward −2, x2 on [10, 13] inside a declared
    // domain [9, 14], so the box, the training range and the standardization all differ.
    internal static BassData Data()
    {
        const int n = 80;
        var x0 = new double[n];
        var x1 = new double[n];
        var x2 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = i * 0.6180339887498949 % 1.0, v = i * 0.4142135623730950 % 1.0, w = i * 0.7320508075688772 % 1.0;
            x0[i] = u;
            x1[i] = -2.0 + 7.0 * v * v;
            x2[i] = 10.0 + 3.0 * w;
            y[i] = Math.Sin(2.0 * (0.6 * x0[i] + 0.1 * x1[i])) + 0.3 * x2[i] + 0.01 * Math.Cos(17.0 * i);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2, new InputDomain(9.0, 14.0))],
                            new ScalarResponse(y));
    }

    internal static double[][] Raw(BassData data) => [.. data.Predictors.Cast<ContinuousPredictor>().Select(p => p.Values.ToArray())];

    // Raw points: the training cases' first ten, and points past the training range on every side, x2's inside its declared domain.
    internal static double[][] Points(BassData data)
    {
        double[][] raw = Raw(data);
        double[][] extra = [[-0.3, 1.4, 0.5, 0.0], [-3.0, 6.5, 1.0, 5.0], [9.2, 13.8, 11.0, 12.9]];
        return [.. raw.Select((column, j) => column.Take(10).Concat(extra[j]).ToArray())];
    }

    internal static CaseInputs Scaled(BassInputTransform transform, double[][] raw)
        => new([.. raw.Select((column, j) => column.Select(x => ((RangeScale)transform.Predictors[j]).Scale(x)).ToArray())], []);

    internal static (BassInputTransform Transform, BassInputs Inputs, BassBasis Basis) Setup(BassData data, int degree, HingeNormalization normalization = HingeNormalization.Box)
    {
        BassInputTransform transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        return (transform, inputs, new BassBasis(inputs, degree, normalization));
    }

    private static BassTerm Term(ProjectionDirection direction, double knot) => new(directional: new DirectionalHinge(direction, knot));

    // The training mean and sample deviation of each column, compensated, apart from the engine's pairwise sums.
    internal static (double Mean, double Deviation)[] Moments(BassData data)
        => [.. Raw(data).Select(column =>
        {
            double mean = Compensated(column) / column.Length;
            return (mean, Math.Sqrt(Compensated([.. column.Select(x => (x - mean) * (x - mean))]) / (column.Length - 1)));
        })];

    internal static double Compensated(IEnumerable<double> values)
    {
        double sum = 0.0, carry = 0.0;
        foreach (double v in values)
        {
            double next = sum + v;
            carry += Math.Abs(sum) >= Math.Abs(v) ? sum - next + v : v - next + sum;
            sum = next;
        }
        return sum + carry;
    }

    // The hand computation: z from the raw values by the column's own mean and deviation, c at the corner of each variable's
    // normalization interval in raw units that maximizes θᵀz, and the guard c = 1 where the hinge is zero over the box.
    internal static (double[] Values, double C) Hand(BassData data, double[][] raw, ProjectionDirection d, double knot, int degree, HingeNormalization normalization)
    {
        (double Mean, double Deviation)[] moments = Moments(data);
        double[][] training = Raw(data);
        double top = 0.0;
        for (int k = 0; k < d.Count; k++)
        {
            int j = d.Variables[k];
            var predictor = (ContinuousPredictor)data.Predictors[j];
            double lower = training[j].Min(), upper = training[j].Max();
            if (normalization == HingeNormalization.Box && predictor.Domain is { } domain) (lower, upper) = (domain.Lower, domain.Upper);
            top += d.Weights[k] * ((d.Weights[k] > 0.0 ? upper : lower) - moments[j].Mean) / moments[j].Deviation;
        }
        double c = Math.Pow(Math.Max(top - knot, 0.0), degree);
        if (c == 0.0) c = 1.0;
        var values = new double[raw[0].Length];
        for (int i = 0; i < values.Length; i++)
        {
            double p = 0.0;
            for (int k = 0; k < d.Count; k++)
            {
                int j = d.Variables[k];
                p += d.Weights[k] * (raw[j][i] - moments[j].Mean) / moments[j].Deviation;
            }
            values[i] = Math.Pow(Math.Max(p - knot, 0.0), degree) / c;
        }
        return (values, c);
    }

    [Fact]
    public void DirectionAndTerm_RefuseWhatIsNotOne()
    {
        Assert.Throws<ArgumentException>(() => new ProjectionDirection([], []));
        Assert.Throws<ArgumentException>(() => new ProjectionDirection([0, 1], [1.0]));
        Assert.Throws<ArgumentException>(() => new ProjectionDirection([1, 1], [0.6, 0.8]));
        Assert.Throws<ArgumentException>(() => new ProjectionDirection([0, 1], [1.0, 0.0]));
        Assert.Throws<ArgumentException>(() => new ProjectionDirection([0, 1], [0.6, 0.81]));
        Assert.Throws<ArgumentException>(() => new ProjectionDirection([0], [double.NaN]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ProjectionDirection([-1], [1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DirectionalHinge(Single, double.PositiveInfinity));

        // Canonical order, value equality bit for bit, and a sign that is the direction's own.
        Assert.Equal(new[] { 0, 2 }, Pair.Variables);
        Assert.Equal(new[] { 0.8, -0.6 }, Pair.Weights);
        Assert.Equal(-0.6, Pair.WeightOn(2));
        Assert.Equal(0.0, Pair.WeightOn(1));
        Assert.Equal(Term(new ProjectionDirection([0, 2], [0.8, -0.6]), 0.25), Term(Pair, 0.25));
        Assert.Equal(Term(new ProjectionDirection([0, 2], [0.8, -0.6]), 0.25).GetHashCode(), Term(Pair, 0.25).GetHashCode());
        Assert.NotEqual(Term(Pair, 0.25), Term(new ProjectionDirection([0, 2], [-0.8, 0.6]), 0.25));
        Assert.NotEqual(Term(Pair, 0.25), Term(Pair, 0.2500000000000001));
        Assert.NotEqual(Term(Pair, 0.25), new BassTerm([new HingeFactor(0, 1, 3)]));
        Assert.Equal("(0.8·z0 − 0.6·z2 − 0.25)₊", Term(Pair, 0.25).ToString());
        Assert.Equal("(−1·z1 + 0.5)₊", Term(new ProjectionDirection([1], [-1.0]), -0.5).ToString());

        // Beside categorical factors alone on the case side, which gate it (slice 5); an output side is its own.
        var hinge = new DirectionalHinge(Single, 0.0);
        Assert.Throws<ArgumentException>(() => new BassTerm([new HingeFactor(0, 1, 3)], directional: hinge));
        BassTerm gated = new(categorical: [new CategoricalFactor(0, [1])], directional: hinge);
        Assert.Equal(hinge, gated.Directional);
        Assert.Equal([new CategoricalFactor(0, [1])], gated.Categorical);
        Assert.Throws<ArgumentException>(() => new BassTerm(legendre: [new LegendreFactor(0, 2)], directional: hinge));
        Assert.Equal(hinge, new BassTerm(functional: [new HingeFactor(0, 1, 2)], directional: hinge).Directional);
    }

    [Theory]
    [InlineData(1, HingeNormalization.Box)]
    [InlineData(2, HingeNormalization.Box)]
    [InlineData(3, HingeNormalization.Box)]
    [InlineData(1, HingeNormalization.TrainingRange)]
    [InlineData(2, HingeNormalization.TrainingRange)]
    public void Columns_MatchAHandComputation(int degree, HingeNormalization normalization)
    {
        BassData data = Data();
        (BassInputTransform transform, _, BassBasis basis) = Setup(data, degree, normalization);
        IDesignKernel<CaseInputs> kernel = basis.Kernels.Case[3];
        double[][] points = Points(data);
        CaseInputs at = Scaled(transform, points);
        double worst = 0.0;
        foreach ((ProjectionDirection d, double knot) in new[] { (Dense, -0.5), (Dense, 0.7), (Single, 0.0), (Pair, -1.2), (Pair, 0.3) })
        {
            BassTerm term = Term(d, knot);
            var column = new double[at.Cases];
            basis.CaseColumn(term, at, column);
            (double[] hand, double c) = Hand(data, points, d, knot, degree, normalization);
            for (int i = 0; i < column.Length; i++) worst = Math.Max(worst, Math.Abs(column[i] - hand[i]) / Math.Max(1.0, Math.Abs(hand[i])));
            Assert.True(Math.Abs(kernel.Scale(term) * c - 1.0) <= 1e-13, $"{term}: 1/c {kernel.Scale(term):R} against c {c:R}");
            Assert.True(hand.Max() > 0.0, $"{term} should be nonzero somewhere");
        }
        output.WriteLine($"degree {degree}, {normalization}: worst |column − hand| / max(1, |hand|) = {worst:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");

        // A hinge zero over the whole box: zero everywhere it is read, and c left at 1.
        BassTerm dead = Term(Single, 10.0);
        var zero = new double[at.Cases];
        basis.CaseColumn(dead, at, zero);
        Assert.All(zero.Take(10), v => Assert.Equal(0.0, v));
        Assert.Equal(1.0, kernel.Scale(dead));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void Derivatives_MatchCentralDifferences_AndTheKinkIsTheMeanOfItsLimits(int degree)
    {
        BassData data = Data();
        (BassInputTransform transform, BassInputs inputs, BassBasis basis) = Setup(data, degree);
        CaseInputs at = Scaled(transform, Points(data));
        double worst = 0.0;
        foreach ((ProjectionDirection d, double knot) in new[] { (Dense, -0.5), (Single, 0.0), (Pair, 0.3) })
        {
            BassTerm term = Term(d, knot);
            for (int v = 0; v < 3; v++)
            {
                var derivative = new double[at.Cases];
                basis.CaseDerivative(term, v, at, derivative);
                if (d.WeightOn(v) == 0.0)
                {
                    Assert.All(derivative, x => Assert.Equal(0.0, x));
                    continue;
                }
                const double h = 1e-6;
                var up = new double[at.Cases];
                var down = new double[at.Cases];
                basis.CaseColumn(term, Shifted(at, v, h), up);
                basis.CaseColumn(term, Shifted(at, v, -h), down);
                double scale = derivative.Max(Math.Abs) + 1.0;
                for (int i = 0; i < at.Cases; i++) worst = Math.Max(worst, Math.Abs((up[i] - down[i]) / (2 * h) - derivative[i]) / scale);
            }
        }
        output.WriteLine($"degree {degree}: worst |difference − derivative| / (max |derivative| + 1) = {worst:E2}");
        Assert.True(worst <= 1e-6, $"{worst:E2}");

        // A knot on case 0's own projection, accumulated as the kernel does: degree one's derivative there is half its slope, θ_v·b_v/(2c);
        // degree two's is 0.
        double p0 = 0.0;
        for (int k = 0; k < Dense.Count; k++)
        {
            (double a, double b) = inputs.Standardizing(Dense.Variables[k])!.Value;
            p0 += Dense.Weights[k] * (a + b * at.Continuous(Dense.Variables[k])[0]);
        }
        BassTerm kinked = Term(Dense, p0);
        double inverseC = basis.Kernels.Case[3].Scale(kinked);
        for (int v = 0; v < 3; v++)
        {
            var derivative = new double[at.Cases];
            basis.CaseDerivative(kinked, v, at, derivative);
            double expected = degree == 1 ? 0.5 * Dense.WeightOn(v) * inputs.Standardizing(v)!.Value.Slope * inverseC : 0.0;
            Assert.Equal(expected, derivative[0], 14);
        }
    }

    // The primitives' seeded sweep: random sparse directions (a from 1 to 3, random variables, weights and signs), knots inside and far
    // outside the projections' range, degrees one to three under both normalizations, at points that include values a thousand times
    // past the box. Columns against the hand computation, and derivatives against the hand derivative q(p − t)₊^(q−1)·θ_v·b_v/c.
    [Fact]
    public void SeededSweep_ColumnsAndDerivatives_MatchTheHandComputation()
    {
        BassData data = Data();
        var rng = new Maths.Rng.Xoshiro256PlusPlus(3601);
        (double Mean, double Deviation)[] moments = Moments(data);
        double[][] points = [.. Points(data).Select((column, j) => column.Concat([1e3 * (j + 1), -2e3 * (j + 1), 0.5 * (j + 1)]).ToArray())];
        double worstValue = 0.0, worstSlope = 0.0;
        for (int c = 0; c < 300; c++)
        {
            int degree = 1 + rng.NextInt(3);
            HingeNormalization normalization = rng.NextDouble() < 0.5 ? HingeNormalization.Box : HingeNormalization.TrainingRange;
            (BassInputTransform transform, BassInputs inputs, BassBasis basis) = Setup(data, degree, normalization);
            int a = 1 + rng.NextInt(3);
            int[] variables = [.. Enumerable.Range(0, 3).OrderBy(_ => rng.NextDouble()).Take(a)];
            double[] weights = [.. variables.Select(_ => (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * (0.05 + rng.NextDouble()))];
            double norm = Math.Sqrt(weights.Sum(w => w * w));
            var d = new ProjectionDirection(variables, weights.Select(w => w / norm));
            double knot = rng.NextDouble() < 0.8 ? -2.5 + 5.0 * rng.NextDouble() : (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * (5.0 + 50.0 * rng.NextDouble());
            BassTerm term = Term(d, knot);
            CaseInputs at = Scaled(transform, points);
            (double[] hand, double normalizer) = Hand(data, points, d, knot, degree, normalization);
            var column = new double[at.Cases];
            basis.CaseColumn(term, at, column);
            for (int i = 0; i < column.Length; i++) worstValue = Math.Max(worstValue, Math.Abs(column[i] - hand[i]) / Math.Max(1.0, Math.Abs(hand[i])));
            for (int v = 0; v < 3; v++)
            {
                var derivative = new double[at.Cases];
                basis.CaseDerivative(term, v, at, derivative);
                var scale = (RangeScale)transform.Predictors[v];
                for (int i = 0; i < at.Cases; i++)
                {
                    double p = 0.0;
                    for (int k = 0; k < d.Count; k++) p += d.Weights[k] * (points[d.Variables[k]][i] - moments[d.Variables[k]].Mean) / moments[d.Variables[k]].Deviation;
                    double outer = p > knot ? degree * Math.Pow(p - knot, degree - 1) : 0.0;
                    double expected = outer * d.WeightOn(v) * (scale.Max - scale.Min) / moments[v].Deviation / normalizer;
                    worstSlope = Math.Max(worstSlope, Math.Abs(derivative[i] - expected) / Math.Max(1.0, Math.Abs(expected)));
                }
            }
        }
        output.WriteLine($"300 cases: worst column {worstValue:E2}, worst derivative {worstSlope:E2} (relative to max(1, |·|))");
        Assert.True(worstValue <= 1e-12, $"columns {worstValue:E2}");
        Assert.True(worstSlope <= 1e-12, $"derivatives {worstSlope:E2}");
    }

    private static CaseInputs Shifted(CaseInputs at, int variable, double by)
    {
        double[][] continuous = [.. Enumerable.Range(0, at.ContinuousCount).Select(v => at.Continuous(v).ToArray())];
        for (int i = 0; i < continuous[variable].Length; i++) continuous[variable][i] += by;
        return new CaseInputs(continuous, []);
    }

    [Fact]
    public void Admission_TakesBothSignsAndSeveralKnots_AndRefusesANearDuplicate()
    {
        BassData data = Data();
        (_, _, BassBasis basis) = Setup(data, 1);
        var design = new BassDesign(basis, data, new BassAdmission(10, SupportKind.Count, 3, null, new BassStructureSpec().AngleFloor, null));
        var negated = new ProjectionDirection([0, 1, 2], [-0.6, 0.48, -0.64]);
        foreach (BassTerm term in new[] { Term(Dense, 0.0), Term(negated, 0.0), Term(Dense, -0.8), Term(Dense, 0.9), Term(negated, 0.4) })
        {
            DesignProposal p = design.ProposeBirth(term);
            output.WriteLine($"{term}: {(p.Admitted ? "admitted" : p.Refusal.ToString())}, newest angle {p.NewestAngle:E3}");
            Assert.True(p.Admitted, $"{term}: {p.Refusal}");
            design.Accept();
        }

        // The same term again, and one whose direction moved by 10⁻⁹: both fall under the angle floor.
        Assert.Equal(BassRefusal.Numerical, design.ProposeBirth(Term(Dense, 0.9)).Refusal);
        double[] moved = [0.6 + 1e-9, -0.48, 0.64];
        double norm = Math.Sqrt(moved.Sum(w => w * w));
        DesignProposal near = design.ProposeBirth(Term(new ProjectionDirection([0, 1, 2], moved.Select(w => w / norm)), 0.9));
        output.WriteLine($"near-duplicate: {near.Refusal}, newest angle {near.NewestAngle:E3}");
        Assert.Equal(BassRefusal.Numerical, near.Refusal);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void ThroughAFit_PredictionAndTheGradientReadTheHinge_TheMomentServiceRefuses_AndE8Samples(int degree)
    {
        BassData data = Data();
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 6, Degree = degree },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 36, BurnIn = 20, Samples = 2 },
        };
        BassFit run = Bass.Fit(config, data);
        var table = new BassTermTable();
        (ProjectionDirection Direction, double Knot)[] hinges = [(Dense, -0.5), (Pair, 0.3)];
        int[] ids = [.. hinges.Select(h => table.Intern(Term(h.Direction, h.Knot)))];
        double[] beta = [0.4, 1.5, -0.7];
        var store = new DrawStore<BassDraw>();
        store.Add(new BassDraw(new BassStructure(ids), 1.0, 1.0, beta, 1.0));
        BassFit fit = run with { Terms = table, Draws = [store], DrawWeights = null };

        double[][] points = Points(data);
        PredictionCases cases = PredictionCases.Of([.. points.Select((column, j) => (BassPredictor)new ContinuousPredictor(data.Predictors[j].Name, column))]);
        double[,,] latent = fit.Predict(data, cases).Latent;
        BassGradient gradient = fit.Gradient(data, cases);
        (double Mean, double Deviation)[] moments = Moments(data);
        double worstValue = 0.0, worstSlope = 0.0;
        for (int i = 0; i < points[0].Length; i++)
        {
            double f = beta[0];
            var slope = new double[3];
            for (int m = 0; m < hinges.Length; m++)
            {
                (ProjectionDirection d, double knot) = hinges[m];
                (double[] hand, double c) = Hand(data, points, d, knot, degree, HingeNormalization.Box);
                f += beta[m + 1] * hand[i];
                double p = 0.0;
                for (int k = 0; k < d.Count; k++) p += d.Weights[k] * (points[d.Variables[k]][i] - moments[d.Variables[k]].Mean) / moments[d.Variables[k]].Deviation;
                double outer = p > knot ? degree * Math.Pow(p - knot, degree - 1) : 0.0;
                for (int j = 0; j < 3; j++) slope[j] += beta[m + 1] * outer * d.WeightOn(j) / moments[j].Deviation / c;
            }
            worstValue = Math.Max(worstValue, Math.Abs(latent[0, i, 0] - f) / Math.Max(1.0, Math.Abs(f)));
            for (int j = 0; j < 3; j++)
                worstSlope = Math.Max(worstSlope, Math.Abs(gradient.Of(data.Predictors[j].Name)[0, i, 0] - slope[j]) / Math.Max(1.0, Math.Abs(slope[j])));
        }
        output.WriteLine($"degree {degree}: prediction {worstValue:E2}, gradient {worstSlope:E2} (relative to max(1, |·|))");
        Assert.True(worstValue <= 1e-12, $"prediction {worstValue:E2}");
        Assert.True(worstSlope <= 1e-12, $"gradient {worstSlope:E2}");

        // Sobol's exact route refuses up front, with its reason, and its automatic route takes the sampled route (the sampled-sensitivity
        // brief, decision 1); E8's automatic route takes the sampled route, the exact route refuses up front.
        ArgumentException refusal = Assert.Throws<ArgumentException>(() => fit.Sensitivity(data, new SensitivitySpec { Route = SensitivityRoute.Exact }));
        Assert.Contains("not a product of one-variable factors", refusal.Message);
        Assert.Contains(table[ids[0]].ToString(), refusal.Message);
        Assert.Equal(SensitivityRoute.Sampled, fit.Sensitivity(data, new SensitivitySpec { Points = 64 }).Route);
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { Route = ActiveSubspaceRoute.Exact }));

        // One hinge alone: every gradient is a multiple of θ's raw-unit image v_j = θ_j/sd_j, so the sampled C is rank one along v.
        var lone = new DrawStore<BassDraw>();
        lone.Add(new BassDraw(new BassStructure([ids[0]]), 1.0, 1.0, [0.4, 1.5], 1.0));
        ActiveSubspaceAnalysis analysis = (fit with { Draws = [lone] }).ActiveSubspace(data, new ActiveSubspaceSpec { Metric = InputMetric.Raw, Points = 256 });
        Assert.Equal(ActiveSubspaceRoute.Sampled, analysis.Route);
        Assert.Equal(0, analysis.Poincare!.Violations);   // E8 slice 8a′: the totals are Sobol's sampled route's, with their errors
        Assert.NotNull(analysis.Poincare.TotalError);
        double[,] meanC = analysis.Posterior.MeanC;
        double[] v = [.. Enumerable.Range(0, 3).Select(j => Dense.WeightOn(j) / moments[j].Deviation)];
        double vNorm = Math.Sqrt(v.Sum(x => x * x));
        for (int j = 0; j < 3; j++) v[j] /= vNorm;
        double along = 0.0, largest = 0.0, residual = 0.0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) along += v[a] * meanC[a, b] * v[b];
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
            {
                largest = Math.Max(largest, Math.Abs(meanC[a, b]));
                residual = Math.Max(residual, Math.Abs(meanC[a, b] - along * v[a] * v[b]));
            }
        output.WriteLine($"degree {degree}: C off θ's image {residual / largest:E2} of its largest entry");
        Assert.True(largest > 0.0 && residual <= 1e-12 * largest, $"{residual / largest:E2}");
    }
}
