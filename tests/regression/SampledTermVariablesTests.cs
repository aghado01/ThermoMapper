using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampled route's rule for a term's variables (<c>BassSensitivity.TermVariablesOf</c>, which the basis fits' route reads too): every
/// variable the term's value depends on, since an input no held term reads has its block skipped and reads zero. For each kind of term and
/// each input the rule leaves out, the term's column on A's rows equals its column on the route's hybrid moving that input alone (an output
/// variable moved alone at its own locations), bit for bit. A term type whose reads are understated fails here, not as a silent zero index.
/// </summary>
public sealed class SampledTermVariablesTests(ITestOutputHelper output)
{
    private const int N = 64;

    private static readonly double[] Knots = [-1.0, -0.4, 0.3, 1.2];

    // One term of each kind on the scalar fixture's x0, x1 and z, each leaving an input out.
    private static (string Kind, BassTerm Term)[] HandTerms() =>
    [
        ("hinge", new([new HingeFactor(0, 1, 3)])),
        ("hinge pair", new([new HingeFactor(0, -1, 6), new HingeFactor(1, 1, 2)])),
        ("hinge with a categorical factor", new([new HingeFactor(1, -1, 8)], [new CategoricalFactor(0, [0, 1])])),
        ("categorical", new(categorical: [new CategoricalFactor(0, [0])])),
        ("Legendre", new(legendre: [new LegendreFactor(0, 2)])),
        ("Legendre pair", new(legendre: [new LegendreFactor(0, 1), new LegendreFactor(1, 3)])),
        ("directional hinge", new(directional: new DirectionalHinge(new ProjectionDirection([0, 1], [0.6, 0.8]), 0.1))),
        ("directional hinge on one input", new(directional: new DirectionalHinge(new ProjectionDirection([1], [1.0]), -0.2))),
        ("gated directional hinge", new(categorical: [new CategoricalFactor(0, [0, 2])], directional: new DirectionalHinge(new ProjectionDirection([0], [1.0]), 0.0))),
        ("ridge column", new(ridge: new RidgeColumn(new ProjectionDirection([1], [1.0]), Knots, 2))),
        ("gated ridge column", new(categorical: [new CategoricalFactor(0, [1])], ridge: new RidgeColumn(new ProjectionDirection([0], [1.0]), Knots, 1))),
    ];

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void EachKindOfTerm_IsUnmovedByEveryInputItsVariablesLeaveOut(int degree)
    {
        BassData data = SensitivityFixture.Scalar();
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        var measure = ProductMeasure.Resolve(transform, inputs);
        (string Kind, BassTerm Term)[] terms = HandTerms();
        var pairs = new Dictionary<string, int>();
        List<string> failures = Check($"degree {degree}", transform, inputs, measure, degree, HingeNormalization.Box, terms, pairs, 4101, varies: true);
        foreach ((string kind, _) in terms) output.WriteLine($"degree {degree}, {kind}: {pairs.GetValueOrDefault(kind)} inputs left out");
        Assert.True(failures.Count == 0, string.Join(Environment.NewLine, failures));
        Assert.All(terms, t => Assert.True(pairs.GetValueOrDefault(t.Kind) > 0, $"{t.Kind} leaves no input out, so nothing is tested"));
    }

    [Fact]
    public void TheFixturesFittedTerms_AreUnmovedByEveryInputTheirVariablesLeaveOut()
    {
        var pairs = new Dictionary<string, int>();
        var failures = new List<string>();
        foreach ((string label, BassFit fit, BassData data) in new[]
                 {
                     ("scalar", SensitivityFixture.FitScalar(), SensitivityFixture.Scalar()),
                     ("functional, integrated", SensitivityFixture.FitFunctional(), SensitivityFixture.Functional()),
                     ("two outputs, integrated", SensitivityFixture.FitTwoOutputs(), SensitivityFixture.TwoOutputs()),
                 })
        {
            BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
            string[] outputs = [.. fit.Transform.OutputVariables.Select(o => o.Name)];
            var measure = ProductMeasure.Resolve(fit.Transform, inputs, integratedOutputs: outputs.Length == 0 ? null : outputs);
            var held = new SortedSet<int>();
            for (int c = 0; c < fit.Draws.Count; c++)
                for (int i = 0; i < fit.Draws[c].Count; i++) held.UnionWith(fit.Draws[c][i].Structure.TermIds);
            failures.AddRange(Check(label, fit.Transform, inputs, measure, fit.Config.Structure.Degree, fit.Config.Structure.ResolvedHingeNormalization,
                                    held.Select(id => (label, fit.Terms[id])), pairs, 4102, varies: false));
            output.WriteLine($"{label}: {held.Count} held terms, {pairs.GetValueOrDefault(label)} (term, input left out) pairs");
            Assert.True(pairs.GetValueOrDefault(label) > 0, $"{label}: no held term leaves an input out, so nothing is tested");
        }
        Assert.True(failures.Count == 0, string.Join(Environment.NewLine, failures));
    }

    // Each term's column on A's rows against the route's hybrid moving one left-out input from B; an output variable is moved alone at
    // locations drawn uniformly over its scaled range. With `varies`, a term's column must vary over A's rows, so that it reads something.
    private static List<string> Check(string label, BassInputTransform transform, BassInputs inputs, ProductMeasure measure, int degree,
                                      HingeNormalization hinges, IEnumerable<(string Kind, BassTerm Term)> terms, Dictionary<string, int> pairs,
                                      int seed, bool varies)
    {
        var moments = BassMoments.On(inputs, degree, measure, hinges);
        var kernels = new DesignKernels(inputs, degree, hinges);
        var basis = new BassBasis(inputs, degree, hinges);
        Func<BassTerm, int[]> of = BassSensitivity.TermVariablesOf(moments, kernels);
        PredictionCases drawn = BassActiveSubspaces.Points(transform, measure, 2 * N, seed);
        var rng = new Xoshiro256PlusPlus(seed);
        double[] Uniforms() => [.. Enumerable.Range(0, N).Select(_ => rng.NextDouble())];
        double[][] locations = [.. transform.OutputVariables.Select(_ => Uniforms())];
        CaseInputs CasesMoving(params int[] moved) =>
            BassSensitivity.Stacked(drawn, measure, N, [(false, moved)]).Prepare(transform, inputs, UnseenLevel.Refuse).Inputs;
        CaseInputs atA = CasesMoving();

        var failures = new List<string>();
        foreach ((string kind, BassTerm term) in terms)
        {
            basis.Validate(term);
            int[] listed = of(term);
            var a = new double[N];
            basis.CaseColumn(term, atA, a);
            if (varies && a.Distinct().Count() < 2) failures.Add($"{label}, {kind} ({term}): its column is constant over A's rows, so it reads nothing");
            for (int v = 0; v < measure.Count; v++)
            {
                if (Array.IndexOf(listed, v) >= 0) continue;
                pairs[kind] = pairs.GetValueOrDefault(kind) + 1;
                double[] before = a, after = new double[N];
                if (measure.OutputOf(v) is int o and >= 0)
                {
                    before = new double[N];
                    basis.OutputColumn(term, new OutputInputs(locations), before);
                    basis.OutputColumn(term, new OutputInputs([.. locations.Select((l, k) => k == o ? Uniforms() : l)]), after);
                }
                else basis.CaseColumn(term, CasesMoving(v), after);
                int row = Enumerable.Range(0, N).FirstOrDefault(r => !(before[r] == after[r]), -1);
                if (row >= 0)
                    failures.Add($"{label}, {kind} ({term}): reads [{string.Join(", ", listed.Select(x => measure.Names[x]))}], yet moving " +
                                 $"{measure.Names[v]} alone moves row {row} from {before[row]:R} to {after[row]:R}");
            }
        }
        return failures;
    }
}
