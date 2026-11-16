using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Sensitivity;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampled-sensitivity brief's slice 3 gates (decision 8): the sampled route carried through a fit of several outputs (E34), a basis
/// fit (D07) and a direct functional fit's forms, each against its exact route on separable fits within the stated errors, and a basis fit
/// without factors against quadrature of its predicted output, a method apart from the estimator.
/// </summary>
public sealed class SampledOutputsTests(ITestOutputHelper output)
{
    private void Report(string name, SampledSensitivityTests.Tally t) =>
        output.WriteLine($"{name}: {t.Cells} cells, {t.Exact} exact zeros; worst √R·|mean z| {t.WorstBias:F2} ({t.WorstBiasCell}); pooled sd of z "
                         + $"{t.PooledSd:F3}, per cell {t.LowestSd:F2} to {t.HighestSd:F2}; beyond 3: {t.Beyond3:P2}; worst z {t.WorstZ:F2} ({t.WorstZCell})");

    private static void Gate(SampledSensitivityTests.Tally t)
    {
        Assert.True(t.WorstBias <= SampledSensitivityTests.BiasCut(t.Cells), $"{t.WorstBiasCell}: √R·|mean z| {t.WorstBias:F2} against {SampledSensitivityTests.BiasCut(t.Cells):F2}");
        Assert.True(t.PooledSd >= SampledSensitivityTests.PooledLow && t.PooledSd <= SampledSensitivityTests.PooledHigh, $"pooled sd of z {t.PooledSd:F3}");
        Assert.True(t.Beyond3 <= 0.01, $"{t.Beyond3:P2} beyond 3");
    }

    // Every draw's pieces of a sampled decomposition against an exact one's, keyed by a label.
    private static void Add(Dictionary<string, SampledSensitivityTests.Piece> cells, string label, SobolDecomposition exact, SobolDecomposition sampled)
    {
        Assert.Equal(exact.Subsets.Select(u => string.Join(",", u)), sampled.Subsets.Select(u => string.Join(",", u)));
        for (int d = 0; d < exact.Count; d++)
        {
            Assert.Equal(exact.IsDefined(d), sampled.IsDefined(d));
            int draw = d;
            double Truth(string piece, int[]? u) => piece switch
            {
                "V" => exact.Variance[draw],
                "V_u" => exact.MainVariance[draw, Index(exact, u!)],
                "S_u" => exact.MainIndex(draw, Index(exact, u!)),
                "T" => exact.TotalVariance[draw, u![0]],
                _ => exact.TotalIndex(draw, u![0]),
            };
            foreach (SampledSensitivityTests.Piece p in SampledSensitivityTests.Pieces(sampled, d, Truth)) cells[$"{label} draw {d} {p.Name}"] = p;
        }
    }

    private static int Index(SobolDecomposition s, int[] u)
    {
        for (int k = 0; k < s.Subsets.Count; k++)
            if (s.Subsets[k].SequenceEqual(u)) return k;
        throw new ArgumentException($"No subset {string.Join(",", u)}.");
    }

    // A cell's scale for the exact-zero test: 1 for an index, else its decomposition's variance at that draw.
    private static Func<string, double> Scales(Func<string, int, double> variance) => cell =>
    {
        string[] parts = cell.Split(' ');
        int at = Array.IndexOf(parts, "draw");
        return parts[at + 2].StartsWith('S') ? 1.0 : variance(string.Join(" ", parts.Take(at)), int.Parse(parts[at + 1]));
    };

    public static TheoryData<MainVarianceForm> Forms() => [MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference];

    // ── Several outputs (E34) ───────────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [MemberData(nameof(Forms))]
    public void SeveralOutputs_SampledMeetsExact_WithinTheStatedError(MainVarianceForm form)
    {
        const int Seeds = 24, N = 2048;
        (BassFit fit, BassData data) = MultiOutputReadoutTests.Fixed();
        double[][] functionals = [[1.0, -0.4], [0.0, 1.0]];
        SobolDecomposition[] exact = [.. functionals.Select(a => BassMultiOutput.Sensitivity(fit, data, a, new SensitivitySpec { Route = SensitivityRoute.Exact }).Decompositions[0])];
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(341, Seeds))
        {
            var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
            for (int f = 0; f < functionals.Length; f++)
            {
                SobolAnalysis sampled = BassMultiOutput.Sensitivity(fit, data, functionals[f], new SensitivitySpec { Route = SensitivityRoute.Sampled, Points = N, Seed = seed, MainVarianceForm = form });
                Assert.Equal(SensitivityRoute.Sampled, sampled.Route);
                Add(cells, $"output {f}", exact[f], sampled.Decompositions[0]);
            }
            runs.Add(cells);
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, Scales((label, d) => exact[label.EndsWith('0') ? 0 : 1].Variance[d]));
        Report($"several outputs, {form}", t);
        Gate(t);
    }

    [Fact]
    public void SeveralOutputs_TheProjectedRoute_IsTheScalarReaders_BitForBit()
    {
        // A directional fit of two outputs from the chain, and E34's separable fixture: the sampled route on a declared output functional
        // reads the scalar fit whose draws carry the projected coefficients, at the same points, bit for bit.
        BassData data = MultiOutputCoefficientTests.Data(2);
        BassConfig directional = MultiOutputCoefficientTests.Config("raw") with
        {
            Structure = new BassStructureSpec { MaxBasis = 4, MaxInteraction = 2, Directional = new DirectionalBasis() },
            Sampling = new BassSampling { Chains = 1, BurnIn = 200, Samples = 10, Thin = 5, MasterSeed = 7 },
        };
        BassFit chain = Bass.Fit(directional, data);
        Assert.Contains(chain.Draws.SelectMany(c => c).SelectMany(d => d.Structure.TermIds), id => chain.Terms[id].Directional is not null);
        foreach ((BassFit fit, BassData d) in new[] { (chain, data), MultiOutputReadoutTests.Fixed() })
        {
            double[] functional = [1.0, -0.4];
            var spec = new SensitivitySpec { Route = SensitivityRoute.Sampled, Points = 512, Seed = 3 };
            SobolDecomposition projected = BassMultiOutput.Sensitivity(fit, d, functional, spec).Decompositions[0];
            if (ReferenceEquals(fit, chain))
            {
                foreach (MainVarianceForm form in new[] { MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference })
                {
                    SobolAnalysis auto = BassMultiOutput.Sensitivity(fit, d, functional, spec with { Route = SensitivityRoute.Auto, MainVarianceForm = form });
                    SobolAnalysis sampled = BassMultiOutput.Sensitivity(fit, d, functional, spec with { MainVarianceForm = form });
                    Assert.Equal(SensitivityRoute.Sampled, auto.Route);
                    Assert.Equal(sampled.Decompositions[0].MainVariance, auto.Decompositions[0].MainVariance);
                    Assert.Equal(sampled.Decompositions[0].MainIndexError, auto.Decompositions[0].MainIndexError);
                    if (form == MainVarianceForm.TwoDifference) Assert.Equal(projected.MainVariance, sampled.Decompositions[0].MainVariance);
                }
            }
            var response = (MultiOutputResponse)d.Response;
            var scalarData = new BassData(d.Predictors, new ScalarResponse([.. Enumerable.Range(0, response.Cases).Select(i => functional[0] * response[i, 0] + functional[1] * response[i, 1])]));
            BassFit scalar = Bass.Fit(fit.Config with { OutputCovariance = null, Sigma2 = Sigma2Prior.InverseGamma(3, 1), Sampling = new BassSampling { Chains = 1, BurnIn = 2, Samples = 1 } }, scalarData);
            var store = new DrawStore<BassDraw>();
            foreach (BassDraw draw in fit.Draws.SelectMany(c => c))
                store.Add(draw with { Coefficients = BassMultiOutput.Project(draw, functional), OutputCoefficients = null, Sigma = null, Sigma2 = 1, OutputRates = null });
            scalar = scalar with { Terms = fit.Terms, Draws = [store], DrawWeights = fit.DrawWeights };
            SobolDecomposition expected = scalar.Sensitivity(scalarData, spec).Decompositions[0];
            Assert.Equal(expected.MainVariance, projected.MainVariance);
            Assert.Equal(expected.TotalVariance, projected.TotalVariance);
            Assert.Equal(expected.MainIndexError, projected.MainIndexError);
            Assert.Equal(expected.TotalVarianceError, projected.TotalVarianceError);
        }
    }

    // ── Basis fits (D07) ──────────────────────────────────────────────────────────────────────────────────────────────────────

    // The basis views are bilinear in two components, so their cells move together and the pooled statistics have about the seeds'
    // degrees of freedom: 24 seeds read the two-difference form's pooled sd at 1.134, and 120 fresh seeds read the same cells at 1.026
    // (probe artifacts/chip-sensitivity/probes/BasisCalibrationProbe.cs); the gate reads 80.
    [Theory]
    [MemberData(nameof(Forms))]
    public void BasisFit_SampledMeetsExact_AtEachLocationAndIntegrated(MainVarianceForm form)
    {
        const int Seeds = 80, N = 1024;
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        var spec = new SensitivitySpec { Draws = [.. Enumerable.Range(0, 20)] };
        BasisSobolAnalysis exact = fit.Sensitivity(spec with { Route = SensitivityRoute.Exact });
        int[] locations = [0, 4, 7, 11, 14];
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(342, Seeds))
        {
            BasisSobolAnalysis sampled = fit.Sensitivity(spec with { Route = SensitivityRoute.Sampled, Points = N, Seed = seed, MainVarianceForm = form });
            Assert.Equal(SensitivityRoute.Sampled, sampled.Route);
            Assert.Equal(N, sampled.Points);
            var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
            foreach (int l in locations) Add(cells, $"location {l}", exact.AtLocation(l), sampled.AtLocation(l));
            Assert.Equal(exact.Integrated.Inputs, sampled.Integrated.Inputs);
            Add(cells, "integrated", exact.Integrated, sampled.Integrated);
            runs.Add(cells);
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, Scales((label, d) =>
            label == "integrated" ? exact.Integrated.Variance[d] : exact.AtLocation(int.Parse(label.Split(' ')[1])).Variance[d]));
        Report($"basis fit, {fit.Components.Count} components, {form}", t);
        Gate(t);
    }

    // A basis fit whose components are directional fits from the chain, on the basis tests' problem.
    private static readonly Lazy<BassBasisFit> DirectionalComponents = new(() =>
    {
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.95, ScaleByColumnSd = true });
        BassConfig config = BassOutputComponentsTests.Config(chains: 1, samples: 8) with
        {
            Structure = new BassStructureSpec { MaxBasis = 5, MaxInteraction = 2, Directional = new DirectionalBasis() },
        };
        return BassComponents.Fit(d, predictors, config);
    });

    // Quadrature of the predicted output over the inputs' boxes, 3-node Gauss–Legendre on `panels` panels per input: per draw and
    // location the closed variance of every subset of the inputs and, integrated, of every subset of the inputs and the location.
    internal sealed class OutputQuadrature
    {
        private readonly double[,,] _latent;   // draw × case × location, the cases the tensor nodes, input 0 slowest
        private readonly double[] _w;
        private readonly int _q;

        public OutputQuadrature(BassBasisFit fit, IReadOnlyList<int> draws, int panels)
        {
            var boxes = fit.Components[0].Transform.Predictors.Select(p => (RangeScale)p).ToArray();
            QuadratureRule[] rules = [.. boxes.Select(b => QuadratureRule.Uniform(b.Min, b.Max, 3, panels))];
            _q = rules[0].Nodes.Length;
            _w = rules[0].Weights;
            double[] a = [.. Enumerable.Range(0, _q * _q).Select(k => rules[0].Nodes[k / _q])], b = [.. Enumerable.Range(0, _q * _q).Select(k => rules[1].Nodes[k % _q])];
            _latent = fit.Predict(PredictionCases.Of([new ContinuousPredictor(boxes[0].Name, a), new ContinuousPredictor(boxes[1].Name, b)]),
                                  new BassPredictOptions { Draws = draws, Predictive = false }).Latent;
        }

        public int Locations => _latent.GetLength(2);

        // τ_u at draw d, u ⊆ {0, 1} at location `at`, or u ⊆ {0, 1, 2} integrated (2 the location, equally weighted) when `at` < 0:
        // one pass sums the dropped coordinates into the kept ones' cells, then the cells' weighted squared deviations.
        public double Closed(int d, int[] u, int at)
        {
            int[] ls = at >= 0 ? [at] : [.. Enumerable.Range(0, Locations)];
            double wl = 1.0 / ls.Length;
            bool k0 = u.Contains(0), k1 = u.Contains(1), kl = at < 0 && u.Contains(2);
            int n0 = k0 ? _q : 1, n1 = k1 ? _q : 1, nl = kl ? ls.Length : 1;
            var cells = new double[n0 * n1 * nl];
            double mean = 0.0;
            for (int i = 0; i < _q; i++)
                for (int j = 0; j < _q; j++)
                    for (int li = 0; li < ls.Length; li++)
                    {
                        double y = _latent[d, i * _q + j, ls[li]];
                        mean += _w[i] * _w[j] * wl * y;
                        double dropped = (k0 ? 1.0 : _w[i]) * (k1 ? 1.0 : _w[j]) * (kl ? 1.0 : wl);
                        cells[((k0 ? i : 0) * n1 + (k1 ? j : 0)) * nl + (kl ? li : 0)] += dropped * y;
                    }
            double tau = 0.0;
            for (int i = 0; i < n0; i++)
                for (int j = 0; j < n1; j++)
                    for (int li = 0; li < nl; li++)
                    {
                        double c = cells[(i * n1 + j) * nl + li] - mean;
                        tau += (k0 ? _w[i] : 1.0) * (k1 ? _w[j] : 1.0) * (kl ? wl : 1.0) * c * c;
                    }
            return tau;
        }

        public double Main(int d, int[] u, int at)
        {
            double s = 0.0;
            for (int mask = 0; mask < 1 << u.Length; mask++)
            {
                int[] v = [.. Enumerable.Range(0, u.Length).Where(k => (mask >> k & 1) != 0).Select(k => u[k])];
                if (v.Length > 0) s += ((u.Length - v.Length) % 2 == 0 ? 1.0 : -1.0) * Closed(d, v, at);
            }
            return s;
        }
    }

    [Fact]
    public void DirectionalBasisFit_SampledMeetsQuadratureOfItsOutput()
    {
        const int Seeds = 30, N = 2048;
        BassBasisFit fit = DirectionalComponents.Value;
        int[] draws = [0, 3, 7];
        var quadrature = new OutputQuadrature(fit, draws, 48);
        var coarse = new OutputQuadrature(fit, draws, 24);
        int[] locations = [0, 5, 10, 14];
        int[] all3 = [0, 1, 2];
        double worstQuadrature = 0.0;
        for (int d = 0; d < draws.Length; d++)
            worstQuadrature = Math.Max(worstQuadrature, Math.Abs(quadrature.Closed(d, all3, -1) - coarse.Closed(d, all3, -1)) / quadrature.Closed(d, all3, -1));
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(343, Seeds))
        {
            BasisSobolAnalysis sampled = fit.Sensitivity(new SensitivitySpec { Draws = draws, Points = N, Seed = seed });
            Assert.Equal(SensitivityRoute.Sampled, sampled.Route);
            var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
            void AddAgainst(string label, SobolDecomposition s, int at)
            {
                int p = at >= 0 ? 2 : 3;
                for (int d = 0; d < s.Count; d++)
                {
                    int draw = d;
                    int[] everything = [.. Enumerable.Range(0, p)];
                    double v = quadrature.Closed(draw, everything, at);
                    double Truth(string piece, int[]? u) => piece switch
                    {
                        "V" => v,
                        "V_u" => quadrature.Main(draw, u!, at),
                        "S_u" => quadrature.Main(draw, u!, at) / v,
                        "T" => v - quadrature.Closed(draw, [.. everything.Where(i => i != u![0])], at),
                        _ => (v - quadrature.Closed(draw, [.. everything.Where(i => i != u![0])], at)) / v,
                    };
                    foreach (SampledSensitivityTests.Piece piece in SampledSensitivityTests.Pieces(s, d, Truth)) cells[$"{label} draw {d} {piece.Name}"] = piece;
                }
            }
            foreach (int l in locations) AddAgainst($"location {l}", sampled.AtLocation(l), l);
            AddAgainst("integrated", sampled.Integrated, -1);
            runs.Add(cells);
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, Scales((label, d) =>
            label == "integrated" ? quadrature.Closed(d, all3, -1) : quadrature.Closed(d, [0, 1], int.Parse(label.Split(' ')[1]))));
        Report($"directional basis fit, {fit.Components.Count} components", t);
        output.WriteLine($"the quadrature's integrated variance moves by at most {worstQuadrature:E2} of itself on halving its panels");
        Assert.True(worstQuadrature <= 1e-3, $"quadrature at {worstQuadrature:E2}");
        Gate(t);
    }

    [Fact]
    public void BasisFit_Route_AndItsRefusals()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        var spec = new SensitivitySpec { Draws = [0, 1, 2] };
        BasisSobolAnalysis auto = fit.Sensitivity(spec), exact = fit.Sensitivity(spec with { Route = SensitivityRoute.Exact });
        Assert.Equal(SensitivityRoute.Exact, auto.Route);
        Assert.Equal(0, auto.Points);
        Assert.Equal(exact.Integrated.MainVariance, auto.Integrated.MainVariance);
        Assert.Equal(exact.AtLocation(3).TotalVariance, auto.AtLocation(3).TotalVariance);

        BasisSobolAnalysis sampled = fit.Sensitivity(spec with { Route = SensitivityRoute.Sampled, Points = 256 });
        Assert.Contains("exact route", Assert.Throws<InvalidOperationException>(() => sampled.ComponentCovariance(0)).Message);
        Assert.NotNull(sampled.AtLocation(0).MainIndexError);
        Assert.Equal(sampled.Integrated.Count, sampled.Count);

        BassBasisFit directional = DirectionalComponents.Value;
        var directionalSpec = new SensitivitySpec { Draws = [0], Points = 64 };
        BasisSobolAnalysis automatic = directional.Sensitivity(directionalSpec);
        Assert.Equal(SensitivityRoute.Sampled, automatic.Route);
        foreach (MainVarianceForm form in new[] { MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference })
        {
            BasisSobolAnalysis autoForm = directional.Sensitivity(directionalSpec with { MainVarianceForm = form });
            BasisSobolAnalysis explicitForm = directional.Sensitivity(directionalSpec with { Route = SensitivityRoute.Sampled, MainVarianceForm = form });
            Assert.Equal(autoForm.Integrated.MainVariance, explicitForm.Integrated.MainVariance);
            Assert.Equal(autoForm.AtLocation(3).MainIndexError, explicitForm.AtLocation(3).MainIndexError);
            if (form == MainVarianceForm.TwoDifference) Assert.Equal(automatic.Integrated.MainVariance, explicitForm.Integrated.MainVariance);
        }
        Assert.Contains("not a product of one-variable factors",
            Assert.Throws<ArgumentException>(() => directional.Sensitivity(new SensitivitySpec { Route = SensitivityRoute.Exact })).Message);
    }

    // ── Basis effects ─────────────────────────────────────────────────────────────────────────────────────────────────────────

    // Every draw's mean and every cell of every effect of a basis view, against the truth's view.
    private static void AddEffects(Dictionary<string, SampledSensitivityTests.Piece> cells, string label, EffectFunctions truth, EffectFunctions sampled)
    {
        Assert.Equal(truth.Effects.Select(e => string.Join(",", e.Inputs)), sampled.Effects.Select(e => string.Join(",", e.Inputs)));
        for (int d = 0; d < truth.Mean.Length; d++)
        {
            cells[$"{label} draw {d} mean"] = new("mean", truth.Mean[d], sampled.Mean[d], sampled.MeanError![d]);
            for (int k = 0; k < truth.Effects.Count; k++)
                for (int c = 0; c < truth.Effects[k].Values.GetLength(1); c++)
                    cells[$"{label} draw {d} f_{string.Join(",", truth.Effects[k].Inputs)}[{c}]"] =
                        new("f", truth.Effects[k].Values[d, c], sampled.Effects[k].Values[d, c], sampled.Effects[k].Errors![d, c]);
        }
    }

    [Fact]
    public void BasisEffects_SampledMeetsExact_PerComponentAtEachLocationAndIntegrated()
    {
        const int Seeds = 24, N = 1024;
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        var spec = new SensitivitySpec { Draws = [.. Enumerable.Range(0, 20)] };
        var grid = new EffectSpec { GridPoints = 6 };
        BasisEffectAnalysis exact = fit.Effects(spec with { Route = SensitivityRoute.Exact }, grid);
        int[] locations = [0, 7, 14];
        double scale = Enumerable.Range(0, fit.Components.Count).Max(k => exact.ComponentEffects(k).Effects.Max(e => e.Values.Cast<double>().Max(Math.Abs)));
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(344, Seeds))
        {
            BasisEffectAnalysis sampled = fit.Effects(spec with { Route = SensitivityRoute.Sampled, Points = N, Seed = seed }, grid);
            Assert.Equal(SensitivityRoute.Sampled, sampled.Route);
            var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
            for (int k = 0; k < fit.Components.Count; k++) AddEffects(cells, $"component {k}", exact.ComponentEffects(k), sampled.ComponentEffects(k));
            foreach (int l in locations) AddEffects(cells, $"location {l}", exact.AtLocation(l), sampled.AtLocation(l));
            AddEffects(cells, "integrated", exact.Integrated(), sampled.Integrated());
            runs.Add(cells);
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, _ => scale);
        Report($"basis effects, {fit.Components.Count} components", t);
        Gate(t);
    }

    [Fact]
    public void DirectionalBasisEffects_SampledMeetsQuadratureOfItsOutput()
    {
        const int Seeds = 30, N = 2048;
        BassBasisFit fit = DirectionalComponents.Value;
        int[] draws = [0, 3, 7];
        var grid = new EffectSpec { GridPoints = 5 };
        var boxes = fit.Components[0].Transform.Predictors.Select(p => (RangeScale)p).ToArray();
        QuadratureRule[] rules = [.. boxes.Select(b => QuadratureRule.Uniform(b.Min, b.Max, 3, 48))];
        int q = rules[0].Nodes.Length;
        double[,,] Predict(double[] a, double[] b) => fit.Predict(PredictionCases.Of([new ContinuousPredictor(boxes[0].Name, a), new ContinuousPredictor(boxes[1].Name, b)]),
                                                                  new BassPredictOptions { Draws = draws, Predictive = false }).Latent;
        // The tensor's mean per draw and location; each input's conditional mean at its grid values; the pair at its cells.
        EffectFunctions shape = fit.Effects(new SensitivitySpec { Draws = draws, Points = 64 }, grid).AtLocation(0);
        double[] g0 = shape.Effects[0].Grids[0].Values, g1 = shape.Effects[1].Grids[0].Values;
        int gm = g0.Length;
        double[,,] full = Predict([.. Enumerable.Range(0, q * q).Select(k => rules[0].Nodes[k / q])], [.. Enumerable.Range(0, q * q).Select(k => rules[1].Nodes[k % q])]);
        double[,,] at0 = Predict([.. Enumerable.Range(0, gm * q).Select(k => g0[k / q])], [.. Enumerable.Range(0, gm * q).Select(k => rules[1].Nodes[k % q])]);
        double[,,] at1 = Predict([.. Enumerable.Range(0, gm * q).Select(k => rules[0].Nodes[k % q])], [.. Enumerable.Range(0, gm * q).Select(k => g1[k / q])]);
        double[,,] pair = Predict([.. Enumerable.Range(0, gm * gm).Select(k => g0[k / gm])], [.. Enumerable.Range(0, gm * gm).Select(k => g1[k % gm])]);
        double Mean(int d, int l) { double s = 0.0; for (int k = 0; k < q * q; k++) s += rules[0].Weights[k / q] * rules[1].Weights[k % q] * full[d, k, l]; return s; }
        double Given(double[,,] at, double[] w, int d, int g, int l) { double s = 0.0; for (int k = 0; k < q; k++) s += w[k] * at[d, g * q + k, l]; return s; }
        EffectFunctions Truth(int l)
        {
            int n = draws.Length;
            var mean = new double[n];
            var e0 = new double[n, gm];
            var e1 = new double[n, gm];
            var e01 = new double[n, gm * gm];
            for (int d = 0; d < n; d++)
            {
                mean[d] = Mean(d, l);
                for (int g = 0; g < gm; g++)
                {
                    e0[d, g] = Given(at0, rules[1].Weights, d, g, l) - mean[d];
                    e1[d, g] = Given(at1, rules[0].Weights, d, g, l) - mean[d];
                }
                for (int c = 0; c < gm * gm; c++) e01[d, c] = pair[d, c, l] - e0[d, c / gm] - e1[d, c % gm] - mean[d];
            }
            return new EffectFunctions(shape.Inputs, [new Effect([0], shape.Effects[0].Grids, e0), new Effect([1], shape.Effects[1].Grids, e1),
                                                       new Effect([0, 1], shape.Effects[2].Grids, e01)], mean);
        }
        int[] locations = [0, 5, 10, 14];
        EffectFunctions[] truths = [.. locations.Select(Truth)];
        double scale = truths.Max(t => t.Effects.Max(e => e.Values.Cast<double>().Max(Math.Abs)));
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(345, Seeds))
        {
            BasisEffectAnalysis sampled = fit.Effects(new SensitivitySpec { Draws = draws, Points = N, Seed = seed }, grid);
            Assert.Equal(SensitivityRoute.Sampled, sampled.Route);
            var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
            for (int i = 0; i < locations.Length; i++) AddEffects(cells, $"location {locations[i]}", truths[i], sampled.AtLocation(locations[i]));
            runs.Add(cells);
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, _ => scale);
        Report($"directional basis effects, {fit.Components.Count} components", t);
        Gate(t);
    }

    [Fact]
    public void BasisEffects_Route()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        var spec = new SensitivitySpec { Draws = [0, 1, 2] };
        var grid = new EffectSpec { GridPoints = 5 };
        BasisEffectAnalysis auto = fit.Effects(spec, grid), exact = fit.Effects(spec with { Route = SensitivityRoute.Exact }, grid);
        Assert.Equal(SensitivityRoute.Exact, auto.Route);
        Assert.Equal(exact.AtLocation(4).Effects[0].Values, auto.AtLocation(4).Effects[0].Values);
        Assert.Null(auto.AtLocation(4).Effects[0].Errors);
        Assert.Null(auto.Integrated().MeanError);
        BassBasisFit directional = DirectionalComponents.Value;
        Assert.Equal(SensitivityRoute.Sampled, directional.Effects(new SensitivitySpec { Draws = [0], Points = 64 }, grid).Route);
        Assert.Contains("not a product of one-variable factors",
            Assert.Throws<ArgumentException>(() => directional.Effects(new SensitivitySpec { Route = SensitivityRoute.Exact }, grid)).Message);
    }

    // ── A direct functional fit's forms ───────────────────────────────────────────────────────────────────────────────────────

    // The functional fixtures' forms: integrated (one decomposition), conditioned on every output variable (one per training coordinate),
    // and the two-output fixture conditioned on t0 alone with t1 integrated.
    private static (BassFit Fit, BassData Data, SensitivitySpec Spec)[] FunctionalForms() =>
    [
        (SensitivityFixture.FitFunctional(), SensitivityFixture.Functional(), new SensitivitySpec { Functional = FunctionalForm.Integrate }),
        (SensitivityFixture.FitFunctional(), SensitivityFixture.Functional(), new SensitivitySpec { Functional = FunctionalForm.Condition }),
        (SensitivityFixture.FitTwoOutputs(), SensitivityFixture.TwoOutputs(), new SensitivitySpec { Functional = FunctionalForm.Condition, ConditionOn = ["t0"] }),
    ];

    [Theory]
    [MemberData(nameof(Forms))]
    public void FunctionalFit_SampledMeetsExact_InEachForm(MainVarianceForm form)
    {
        const int Seeds = 24, N = 1024;
        var cases = FunctionalForms();
        SobolAnalysis[] exact = [.. cases.Select(c => BassSensitivity.Of(c.Fit, c.Data, c.Spec with { Route = SensitivityRoute.Exact, Draws = [.. Enumerable.Range(0, 20)] }))];
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(346, Seeds))
        {
            var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
            for (int i = 0; i < cases.Length; i++)
            {
                SobolAnalysis sampled = cases[i].Fit.Sensitivity(cases[i].Data,
                    cases[i].Spec with { Route = SensitivityRoute.Sampled, Draws = [.. Enumerable.Range(0, 20)], Points = N, Seed = seed, MainVarianceForm = form });
                Assert.Equal(SensitivityRoute.Sampled, sampled.Route);
                Assert.Equal(exact[i].Decompositions.Count, sampled.Decompositions.Count);
                Assert.Equal(exact[i].Form, sampled.Form);
                for (int c = 0; c < exact[i].Decompositions.Count; c++) Add(cells, $"form {i} surface {c}", exact[i].Decompositions[c], sampled.Decompositions[c]);
            }
            runs.Add(cells);
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, Scales((label, d) =>
        {
            string[] parts = label.Split(' ');
            return exact[int.Parse(parts[1])].Decompositions[int.Parse(parts[3])].Variance[d];
        }));
        Report($"functional forms, {form}", t);
        Gate(t);
    }

    [Fact]
    public void FunctionalFit_Effects_SampledMeetExact_InEachForm()
    {
        const int Seeds = 24, N = 1024;
        var cases = FunctionalForms();
        var grid = new EffectSpec { GridPoints = 5 };
        EffectAnalysis[] exact = [.. cases.Select(c => BassSensitivity.Effects(c.Fit, c.Data, c.Spec with { Route = SensitivityRoute.Exact, Draws = [.. Enumerable.Range(0, 20)] }, grid))];
        double scale = exact.Max(a => a.At.Max(at => at.Effects.Max(e => e.Values.Cast<double>().Max(Math.Abs))));
        var runs = new List<Dictionary<string, SampledSensitivityTests.Piece>>();
        foreach (int seed in SeedTree.Derive(347, Seeds))
        {
            var cells = new Dictionary<string, SampledSensitivityTests.Piece>();
            for (int i = 0; i < cases.Length; i++)
            {
                EffectAnalysis sampled = BassSensitivity.Effects(cases[i].Fit, cases[i].Data,
                    cases[i].Spec with { Route = SensitivityRoute.Sampled, Draws = [.. Enumerable.Range(0, 20)], Points = N, Seed = seed }, grid);
                Assert.Equal(SensitivityRoute.Sampled, sampled.Route);
                Assert.Equal(exact[i].At.Count, sampled.At.Count);
                for (int c = 0; c < exact[i].At.Count; c++) AddEffects(cells, $"form {i} surface {c}", exact[i].At[c], sampled.At[c]);
            }
            runs.Add(cells);
        }
        SampledSensitivityTests.Tally t = SampledSensitivityTests.Summarize(runs, _ => scale);
        Report("functional effects", t);
        Gate(t);
    }
}
