using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Intrinsics;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Bgp;
using Maths.Regression.Changepoint;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Tests;

/// <summary>
/// The engines' per-iteration split (the brief's <c>engine:*</c> blocks): each engine's whole step as production runs it,
/// and its parts as the public pieces it calls or as verbatim test-side copies of its private loops, each with the counted
/// units the code executes.
/// </summary>
internal static class CostEngines
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double Sink { get => CostKernels.Sink; set => CostKernels.Sink = value; }
    private static object? ObjectSink { set => CostKernels.ObjectSink = value; }

    // ================================================================ BARS (and ABMRS)

    internal sealed class BarsFixture
    {
        public readonly double[] X, Y;
        public readonly SplineBasis Basis = new(3);
        public readonly WeightedNormalModel Model = new();
        public readonly IComplexityPrior Prior;
        public readonly SplineTarget Target;
        public readonly LocalBetaKernel Kernel = new(50.0);

        public BarsFixture(int n, IComplexityPrior? prior = null)
        {
            var rng = new Xoshiro256PlusPlus(12 + n);
            X = Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();
            Y = X.Select(t => Math.Sin(2.0 * Math.PI * t) + 0.15 * Gaussian(rng)).ToArray();
            Prior = prior ?? new PoissonPrior(5.0);
            Target = new SplineTarget(Basis, Model, Prior, X, Y);
        }

        public IRjMove<KnotConfig>[] Moves()
        {
            var density = new KnotBirthDensity(Kernel, 1.0);
            return [new KnotBirthMove(density), new KnotDeathMove(density), new KnotRelocateMove(Kernel)];
        }

        public static KnotConfig Equispaced(int k) => new(Enumerable.Range(1, k).Select(i => i / (k + 1.0)).ToArray());

        public ReversibleJumpChain<KnotConfig> Chain(int seed) => new(Moves(), Target, Equispaced(5), new Xoshiro256PlusPlus(seed));
    }

    // The support gate, the centring and the median as the target and the model compute them (private there).
    private static bool HasSupport(double[,] design, double minimum)
    {
        int m = design.GetLength(0), nu = design.GetLength(1);
        for (int j = 0; j < nu; j++)
        {
            double mass = 0.0;
            for (int i = 0; i < m; i++) mass += design[i, j];
            if (mass < minimum) return false;
        }
        return true;
    }

    private static double Median(double[] y)
    {
        var sorted = (double[])y.Clone();
        Array.Sort(sorted);
        int half = sorted.Length / 2;
        return sorted.Length % 2 == 1 ? sorted[half] : 0.5 * (sorted[half - 1] + sorted[half]);
    }

    private static List<CostProfile.Part> BarsEvaluationParts(BarsFixture fx, KnotConfig cfg)
    {
        double[] x = fx.X, y = fx.Y;
        int n = x.Length;
        double[,] design = fx.Basis.Design(cfg, x);
        int nu = design.GetLength(1);
        var bd = new BandedDesign(design);
        int bw = bd.Bandwidth;
        var band = new double[bw + 1, nu];
        var b = new double[nu];
        bd.Accumulate(null, band, y, b);
        var gram = new BandCholesky(nu, bw, BandFactorization.Ldlt, strict: true);
        gram.DecomposeBanded(band);
        double level = Median(y);
        var unit = new double[nu];
        Array.Fill(unit, 1.0);
        var rowSums = new double[n];
        var response = new double[n];
        for (int i = 0; i < n; i++) response[i] = y[i] - level;
        var rhs = new double[nu];
        bd.AccumulateRhs(null, response, rhs);
        double[] beta = gram.Solve(rhs);
        var fit = new double[n];
        bd.Apply(beta, fit);
        double rss = 0.0, fitted = 0.0;
        for (int i = 0; i < n; i++) { double r = response[i] - fit[i]; rss += r * r; fitted += fit[i] * fit[i]; }
        var priorMean = new double[nu];
        Array.Fill(priorMean, level);
        (InverseGamma? sigma2, _) = fx.Model.ResolvedPriors(y);
        double tau = 1.0 / n;
        double[] knots = BSpline.MakeClampedKnots(cfg.InteriorKnots, 3);
        long basisTerms = 4L * n, accumulateTerms = (long)n * (bw + 1) * (bw + 2) / 2;
        long factorTerms = CostKernels.BandFactorTerms(nu, bw);
        var checkBand = new double[bw + 1, nu];
        var checkB = new double[nu];

        return
        [
            new("evaluate_whole", c => { for (int i = 0; i < c; i++) ObjectSink = fx.Target.Evaluate(cfg); }, 1, "evaluation"),
            new("design_build", c => { for (int i = 0; i < c; i++) ObjectSink = fx.Basis.Design(cfg, x); }, basisTerms, "basis value"),
            new("evaluate_basis_rows_only", c => { for (int i = 0; i < c; i++) for (int r = 0; r < n; r++) ObjectSink = BSpline.EvaluateBasis(x[r], knots, 3); }, basisTerms, "basis value"),
            new("support_gate", c => { for (int i = 0; i < c; i++) Sink = HasSupport(design, 1.0) ? 1 : 0; }, (double)n * nu, "element"),
            new("banded_design_spans", c => { for (int i = 0; i < c; i++) ObjectSink = new BandedDesign(design); }, (double)n * nu, "element"),
            new("response_median_resolve", c => { for (int i = 0; i < c; i++) { double lv = Median(y); double s = 0; foreach (double v in y) s += (v - lv) * (v - lv); Sink = s; } }, n, "case (sort n log n)"),
            new("accumulate_gram_and_rhs", c => { for (int i = 0; i < c; i++) { Array.Clear(checkBand); Array.Clear(checkB); bd.Accumulate(null, checkBand, y, checkB); } }, accumulateTerms, "madd"),
            new("factor_strict_ldlt", c => { for (int i = 0; i < c; i++) gram.DecomposeBanded(band); }, factorTerms, "term"),
            new("admission_condition_estimate", c => { for (int i = 0; i < c; i++) Sink = DesignAdmission.Admits(gram) ? 1 : 0; }, nu, "pivot"),
            new("centre_apply_and_rowsum_scan", c => { for (int i = 0; i < c; i++) { bd.Apply(unit, rowSums); bool ok = true; for (int r = 0; r < n && ok; r++) ok = Math.Abs(rowSums[r] - 1.0) <= 1e-12; var res = new double[n]; for (int r = 0; r < n; r++) res[r] = y[r] - level; ObjectSink = res; Sink = ok ? 1 : 0; } }, (double)n * (bw + 1) + 2.0 * n, "madd"),
            new("rhs_and_solve", c => { for (int i = 0; i < c; i++) { Array.Clear(rhs); bd.AccumulateRhs(null, response, rhs); Sink = gram.Solve(rhs)[0]; } }, (double)n * (bw + 1) + CostKernels.BandSolveTerms(nu, bw), "madd"),
            new("fit_and_residual", c => { for (int i = 0; i < c; i++) { bd.Apply(beta, fit); double s = 0, f2 = 0; for (int r = 0; r < n; r++) { double e = response[r] - fit[r]; s += e * e; f2 += fit[r] * fit[r]; } Sink = s + f2; } }, (double)n * (bw + 1) + 4.0 * n, "madd"),
            new("collapsed_g_prior_block", c => { for (int i = 0; i < c; i++) { CollapsedGPrior blk = CollapsedGPrior.FromLeastSquares(gram, beta, rss, fitted, n, priorMean); GPriorMarginal mg = blk.MarginalAt(tau); Sink = mg.LogMarginal(sigma2); } }, 1, "block"),
        ];
    }

    private static List<CostProfile.Part> BarsStepParts(BarsFixture fx, ReversibleJumpChain<KnotConfig> chain, KnotConfig cfg)
    {
        IRjMove<KnotConfig>[] moves = fx.Moves();
        var byKey = moves.ToDictionary(m => m.Key);
        var rng = new Xoshiro256PlusPlus(77);
        double total = moves.Sum(m => m.Weight(cfg));
        return
        [
            new("step_whole", c => { for (int i = 0; i < c; i++) chain.Step(); }, 1, "step"),
            new("total_weight_once", c => { for (int i = 0; i < c; i++) { double t = 0; for (int m = 0; m < moves.Length; m++) { double w = moves[m].Weight(cfg); if (w > 0) t += w; } Sink = t; } }, moves.Length, "weight call"),
            new("select_move", c => { for (int i = 0; i < c; i++) { double u = rng.NextDouble() * total, cum = 0; int pick = moves.Length - 1; for (int m = 0; m < moves.Length; m++) { double w = moves[m].Weight(cfg); if (w <= 0) continue; cum += w; if (u < cum) { pick = m; break; } } Sink = pick; } }, moves.Length, "weight call"),
            new("reverse_key_lookup", c => { for (int i = 0; i < c; i++) Sink = byKey[moves[i % moves.Length].ReverseKey].Weight(cfg); }, 1, "lookup"),
            new("propose_birth", c => { for (int i = 0; i < c; i++) ObjectSink = moves[0].Propose(cfg, rng); }, 1, "proposal"),
            new("propose_death", c => { for (int i = 0; i < c; i++) ObjectSink = moves[1].Propose(cfg, rng); }, 1, "proposal"),
            new("propose_relocate", c => { for (int i = 0; i < c; i++) ObjectSink = moves[2].Propose(cfg, rng); }, 1, "proposal"),
            new("acceptance_arithmetic", c => { for (int i = 0; i < c; i++) { double la = Math.Log(0.3 / 0.7) - Math.Log(0.4 / 0.8) + i * 1e-9; Sink = rng.NextDouble() < Math.Exp(la) ? 1 : 0; } }, 1, "step"),
        ];
    }

    /// <summary>The evaluation's parts at forced configurations, and the step's parts at the chain's realized state.</summary>
    internal static void BarsBlock()
    {
        const string block = "engine:bars";
        foreach (int n in new[] { 100, 500, 2000, 10_000 })
        {
            var fx = new BarsFixture(n);
            foreach (int k in new[] { 3, 8, 16, 32 })
            {
                KnotConfig cfg = BarsFixture.Equispaced(k);
                if (fx.Target.Evaluate(cfg).LogPrior == double.NegativeInfinity) continue;
                CostProfile.Point(block, $"evaluation n={n} k={k}", BarsEvaluationParts(fx, cfg));
            }

            var chain = fx.Chain(21);
            for (int i = 0; i < 2000; i++) chain.Step();
            long attempts0 = chain.Attempts;
            double kSum = 0;
            for (int i = 0; i < 500; i++) { chain.Step(); kSum += chain.Current.InteriorKnots.Length; }
            double kReal = kSum / 500;
            CostProfile.Row(block, $"step n={n}", "(count)", "realized_k", kReal);
            CostProfile.Row(block, $"step n={n}", "(count)", "attempts_per_step", (chain.Attempts - attempts0) / 500.0);
            CostProfile.Row(block, $"step n={n}", "(count)", "accepted_per_step", chain.Accepted / (double)Math.Max(1, chain.Attempts));
            CostProfile.Point(block, $"step n={n}", BarsStepParts(fx, chain, chain.Current));
        }

        // The run: one chain, burn-in 500 and 1,500 samples, at n = 500.
        {
            var fx = new BarsFixture(500);
            double[] grid = Enumerable.Range(0, 101).Select(i => i / 100.0).ToArray();
            var cfgRun = new BarsConfig { Chains = 1, BurnIn = 500, MaxSamples = 1500, MasterSeed = 5 };
            const int iterations = 2000;
            CostProfile.Point(block, "run n=500 poisson", new List<CostProfile.Part>
            {
                new("bars_run", c => { for (int i = 0; i < c; i++) ObjectSink = Bars.Run(cfgRun, fx.X, fx.Y, grid); }, iterations, "iteration"),
            });
        }

        // ABMRS: BARS with its EBIC count law, the whole step and the prior's call.
        foreach (int candidates in new[] { 100, 1000 })
        {
            var prior = new EbicPrior(1.0, candidates);
            var fx = new BarsFixture(500, prior);
            var chain = fx.Chain(31);
            for (int i = 0; i < 2000; i++) chain.Step();
            long a0 = chain.Attempts;
            CostProfile.Point(block, $"abmrs n=500 C={candidates}", new List<CostProfile.Part>
            {
                new("step_whole", c => { for (int i = 0; i < c; i++) chain.Step(); }, 1, "step"),
                new("ebic_log_prior", c => { for (int i = 0; i < c; i++) Sink = prior.LogPrior(1 + i % 10); }, 1, "call"),
            });
            CostProfile.Row(block, $"abmrs n=500 C={candidates}", "(count)", "attempts_per_step", (chain.Attempts - a0) / (double)Math.Max(1, chain.Attempts));
        }
    }

    // ================================================================ BAPS

    private sealed class BapsFixture
    {
        public double[,] Z = null!;
        public double[] Y = null!;
        public IBandPenalty Penalty = null!;
        public int N, Nu, Bandwidth;
        public TensorDesign? Tensor;
    }

    private static BapsFixture BapsOneDimensional(int n, int nu)
    {
        var rng = new Xoshiro256PlusPlus(91 + n);
        var x = Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();
        var y = x.Select(t => Math.Sin(2.0 * Math.PI * t) + 0.15 * Gaussian(rng)).ToArray();
        var knots = BarsFixture.Equispaced(nu - 4);
        double[,] z = new SplineBasis(3).Design(knots, x);
        var penalty = new DifferencePenalty(2);
        return new BapsFixture { Z = z, Y = y, Penalty = penalty, N = n, Nu = z.GetLength(1), Bandwidth = Math.Max(new BandedDesign(z).Bandwidth, penalty.Bandwidth) };
    }

    private static BapsFixture BapsTensor(int knots)
    {
        int g = 3 * (knots + 4) / 2;
        var rng = new Xoshiro256PlusPlus(131 + knots);
        var xs = new double[g * g];
        var ys = new double[g * g];
        var y = new double[g * g];
        for (int a = 0; a < g; a++)
            for (int b = 0; b < g; b++)
            {
                int i = a * g + b;
                xs[i] = (a + 0.5) / g; ys[i] = (b + 0.5) / g;
                y[i] = Math.Sin(1.5 * Math.PI * xs[i]) * Math.Sin(5.0 * Math.PI * ys[i]) + 0.1 * Gaussian(rng);
            }
        var k = BarsFixture.Equispaced(knots);
        TensorDesign td = TensorProductDesign.Build(new SplineBasis(3), k, xs, new SplineBasis(3), k, ys);
        var penalty = new TensorPenalty(td.NuX, td.NuY);
        return new BapsFixture
        {
            Z = td.Design, Y = y, Penalty = penalty, N = y.Length, Nu = td.NuX * td.NuY, Tensor = td,
            Bandwidth = Math.Max(new BandedDesign(td.Design).Bandwidth, penalty.Bandwidth),
        };
    }

    private static List<CostProfile.Part> BapsParts(BapsFixture fx, Baps baps)
    {
        double[,] z = fx.Z;
        double[] y = fx.Y;
        int n = fx.N, nu = fx.Nu, bw = fx.Bandwidth;
        var bd = new BandedDesign(z);
        var gram0 = new double[bw + 1, nu];
        var zty = new double[nu];
        bd.Accumulate(null, gram0, y, zty);
        var band = new double[bw + 1, nu];
        var factor = new BandCholesky(nu, bw, BandFactorization.Ldlt);
        var gaussian = new GaussianConditional(nu);
        var rng = new Xoshiro256PlusPlus(5);
        var beta = new double[nu];
        var fit = new double[n];
        double lambda = 3.0;
        Array.Copy(gram0, band, gram0.Length);
        fx.Penalty.AccumulateInto(band, nu, lambda);
        factor.DecomposeBanded(band);
        gaussian.Draw(factor, zty, rng, beta, 0.3);
        long factorTerms = CostKernels.BandFactorTerms(nu, bw);
        var parts = new List<CostProfile.Part>
        {
            new("gram_copy", c => { for (int i = 0; i < c; i++) Array.Copy(gram0, band, gram0.Length); }, gram0.Length, "double"),
            new("penalty_accumulate", c => { for (int i = 0; i < c; i++) fx.Penalty.AccumulateInto(band, nu, lambda); }, nu * 3, "element"),
            new("band_factor_floored_ldlt", c => { for (int i = 0; i < c; i++) factor.DecomposeBanded(band); }, factorTerms, "term"),
            new("gaussian_draw", c => { for (int i = 0; i < c; i++) gaussian.Draw(factor, zty, rng, beta, 0.3); }, CostKernels.BandSolveTerms(nu, bw) + (double)nu * (bw + 1), "term"),
            new("residual_ss_dense_loop", c => { for (int i = 0; i < c; i++) { double ss = 0; for (int r = 0; r < n; r++) { double f = 0; for (int j = 0; j < nu; j++) f += z[r, j] * beta[j]; double d = y[r] - f; ss += d * d; } Sink = ss; } }, (double)n * nu, "madd"),
            new("residual_ss_banded_apply", c => { for (int i = 0; i < c; i++) { bd.Apply(beta, fit); double ss = 0; for (int r = 0; r < n; r++) { double d = y[r] - fit[r]; ss += d * d; } Sink = ss; } }, (double)n * (bw + 1), "madd"),
            new("roughness", c => { for (int i = 0; i < c; i++) Sink = fx.Penalty.Roughness(beta); }, nu * 3, "element"),
            new("conjugate_variance_draw", c => { for (int i = 0; i < c; i++) Sink = ConjugateUpdates.NormalVariance(1e-3, 1e-3, n, 3.0).Sample(rng); }, 1, "draw"),
            new("reml_evidence_per_lambda_proposal", c => { var ps = new PenalizedSpline(z, y, fx.Penalty); for (int i = 0; i < c; i++) Sink = ps.RemlLogEvidence(lambda * (1 + 1e-3 * (i % 7))); }, factorTerms, "term"),
        };
        if (fx.Penalty is ILocallyAdaptivePenalty local)
        {
            var weights = new double[local.DifferenceCount(nu)];
            Array.Fill(weights, lambda);
            var sq = new double[weights.Length];
            var eta = new double[weights.Length];
            parts.Add(new("local_penalty_accumulate", c => { for (int i = 0; i < c; i++) local.AccumulateInto(band, nu, weights); }, weights.Length * 3, "element"));
            parts.Add(new("squared_differences", c => { for (int i = 0; i < c; i++) local.SquaredDifferencesInto(beta, sq); }, weights.Length * 3, "element"));
            parts.Add(new("local_multiplier_draws", c => { for (int i = 0; i < c; i++) for (int j = 0; j < eta.Length; j++) eta[j] = ConjugateUpdates.NormalPrecision(1.0, 1.0, 1.0, 0.2 + sq[j]).Sample(rng); }, weights.Length, "draw"));
        }
        if (fx.Penalty is IAnisotropicPenalty aniso)
        {
            var lambdas = new double[aniso.Components];
            Array.Fill(lambdas, lambda);
            var logTau2 = new double[aniso.Components];
            var smoothing = new ComponentSmoothing(aniso, baps.ComponentPrior!, SmoothingVarianceUpdate.Slice);
            parts.Add(new("component_penalty_accumulate", c => { for (int i = 0; i < c; i++) aniso.AccumulateComponents(band, nu, lambdas); }, nu * 6, "element"));
            parts.Add(new("component_smoothing_update", c => { for (int i = 0; i < c; i++) smoothing.Update(logTau2, beta, rng); }, 1, "update"));
        }
        return parts;
    }

    internal static void BapsBlock()
    {
        const string block = "engine:baps";
        var points = new List<(string Name, BapsFixture Fx, BapsScheme[] Schemes)>
        {
            ("1-D n=120 nu=29", BapsOneDimensional(120, 29), [BapsScheme.Collapsed, BapsScheme.Gibbs, BapsScheme.LocallyAdaptiveGibbs]),
            ("1-D n=1000 nu=104", BapsOneDimensional(1000, 104), [BapsScheme.Collapsed, BapsScheme.Gibbs, BapsScheme.LocallyAdaptiveGibbs]),
            ("1-D n=5000 nu=504", BapsOneDimensional(5000, 504), [BapsScheme.Collapsed, BapsScheme.Gibbs, BapsScheme.LocallyAdaptiveGibbs]),
            ("tensor 6 knots", BapsTensor(6), [BapsScheme.AnisotropicGibbs, BapsScheme.Gibbs]),
            ("tensor 10 knots", BapsTensor(10), [BapsScheme.AnisotropicGibbs, BapsScheme.Gibbs]),
        };
        foreach ((string name, BapsFixture fx, BapsScheme[] schemes) in points)
        {
            var baps = new Baps(fx.Z, fx.Y, fx.Penalty);
            string label = $"{name} (n={fx.N} nu={fx.Nu} bw={fx.Bandwidth})";
            CostProfile.Row(block, label, "(count)", "bandwidth", fx.Bandwidth);
            CostProfile.Row(block, label, "(count)", "design_bytes", 8.0 * fx.N * fx.Nu);
            const int burn = 200, samples = 600;
            var runs = schemes.Select(s => new CostProfile.Part($"run_{s}", c => { for (int i = 0; i < c; i++) ObjectSink = baps.Run(s, chains: 1, masterSeed: 7, burn: burn, samples: samples); },
                                                                 burn + samples, "iteration")).ToList();
            CostProfile.Point(block, label + " run", runs);
            CostProfile.Point(block, label + " parts", BapsParts(fx, baps));
        }

        // The hierarchical extension: its whole run, by scheme, on six replicates of the 1-D fixture.
        {
            BapsFixture one = BapsOneDimensional(40, 29);
            var rng = new Xoshiro256PlusPlus(41);
            var replicates = Enumerable.Range(0, 6).Select(_ => { double shift = 0.3 * (rng.NextDouble() - 0.5); return one.Y.Select((v, i) => v + shift + 0.05 * Gaussian(rng)).ToArray(); }).ToArray();
            var h = new HierarchicalPSpline(one.Z, replicates, new DifferencePenalty(2));
            const int burn = 200, samples = 600;
            var runs = new[] { HierarchicalScheme.Centred, HierarchicalScheme.CollapsedPopulation }
                .Select(s => new CostProfile.Part($"hierarchical_run_{s}", c => { for (int i = 0; i < c; i++) ObjectSink = h.Run(chains: 1, masterSeed: 3, burn: burn, samples: samples, scheme: s); }, burn + samples, "iteration")).ToList();
            CostProfile.Point(block, "hierarchical n=40 nu=29 R=6", runs);
        }
    }

    // ================================================================ BGP

    private static (double[,] X, double[] Y) SwissRoll(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n, 3];
        var f = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = rng.NextDouble(), v = rng.NextDouble();
            double phi = 1.5 * Math.PI * (1.0 + 2.0 * u);
            x[i, 0] = phi * Math.Cos(phi); x[i, 1] = phi * Math.Sin(phi); x[i, 2] = 12.0 * v;
            f[i] = Math.Sin(2.0 * Math.PI * u) + 0.5 * Math.Cos(2.0 * Math.PI * v) + 0.1 * Gaussian(rng);
        }
        return (x, f);
    }

    internal static void BgpBlock()
    {
        const string block = "engine:bgp";
        foreach (int n in new[] { 50, 100, 200, 400, 800, 1600 })
        {
            (double[,] x, double[] y) = SwissRoll(n, 5 + n);
            var kernel = new SquaredExponentialKernel();
            var gp = new GpRegression(x, y, 0.01, kernel);
            const double t = 0.05;

            // The kernel matrix as the fit builds it, on squared distances of the unit-scaled inputs.
            var sq = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < i; j++)
                {
                    double s = 0;
                    for (int c = 0; c < 3; c++) { double d = (x[i, c] - x[j, c]) / 12.0; s += d * d; }
                    sq[i, j] = s;
                }
            var a = new double[n, n];
            var l = new double[n, n];
            var lInv = new double[n, n];
            var chol = new CholeskyDecomposition(n);
            var strict = new DenseCholesky(n);
            double[] row = new double[n], rowOut = new double[n];
            long fill = (long)n * (n + 1) / 2;
            void Fill()
            {
                for (int i = 0; i < n; i++)
                {
                    a[i, i] = kernel.Evaluate(0.0, t) + 0.01;
                    for (int j = 0; j < i; j++) { double v = kernel.Evaluate(sq[i, j], t); a[i, j] = v; a[j, i] = v; }
                }
            }
            Fill();
            chol.Decompose(a);
            CostKernels.CholeskyStep1(a, l, n);
            var parts = new List<CostProfile.Part>
            {
                new("fit_whole", c => { for (int i = 0; i < c; i++) Sink = gp.Fit(t).LogMarginal; }, 1, "fit"),
                new("kernel_fill_interface_call", c => { for (int i = 0; i < c; i++) Fill(); }, fill, "pair"),
                new("kernel_fill_inline_exp", c => { for (int r = 0; r < c; r++) for (int i = 0; i < n; i++) { a[i, i] = 1.01; for (int j = 0; j < i; j++) { double v = Math.Exp(-sq[i, j] / (2.0 * t)); a[i, j] = v; a[j, i] = v; } } }, fill, "pair"),
                new("kernel_fill_vector_exp_lower_triangle", c => { for (int r = 0; r < c; r++) for (int i = 0; i < n; i++) { for (int j = 0; j < i; j++) row[j] = -sq[i, j] / (2.0 * t); VectorExp(row, rowOut, i); for (int j = 0; j < i; j++) { a[i, j] = rowOut[j]; a[j, i] = rowOut[j]; } a[i, i] = 1.01; } }, fill, "pair"),
                new("cholesky_decompose_production", c => { for (int i = 0; i < c; i++) chol.Decompose(a); }, CostKernels.CholeskyStep1Terms(n) + CostKernels.CholeskyStep2Terms(n), "term"),
                new("cholesky_step1_factor", c => { for (int i = 0; i < c; i++) CostKernels.CholeskyStep1(a, l, n); }, CostKernels.CholeskyStep1Terms(n), "term"),
                new("cholesky_step2_inverse_of_L", c => { for (int i = 0; i < c; i++) CostKernels.CholeskyStep2(l, lInv, n); }, CostKernels.CholeskyStep2Terms(n), "term"),
                new("dense_cholesky_strict_factor", c => { for (int i = 0; i < c; i++) strict.Decompose(a); }, CostKernels.DenseFactorTerms(n), "term"),
                new("two_solves", c => { for (int i = 0; i < c; i++) { Sink = chol.Solve(y)[0]; Sink = chol.Solve(y)[1]; } }, 2.0 * n * (n + 1), "term"),
                new("three_n_by_n_allocations", c => { for (int i = 0; i < c; i++) { ObjectSink = new double[n, n]; ObjectSink = new double[n, n]; ObjectSink = new double[n, n]; } }, 3, "array"),
            };
            CostProfile.Point(block, $"fit n={n}", parts);
        }

        // The sampler: one chain, burn-in 100 and 400 draws, a conditional variance at 100 test points per draw.
        foreach (int n in new[] { 100, 400 })
        {
            (double[,] x, double[] y) = SwissRoll(n, 5 + n);
            var gp = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
            var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp));
            (double[,] xt, _) = SwissRoll(100, 999);
            CostProfile.Point(block, $"sampler run n={n}", new List<CostProfile.Part>
            {
                new("bgp_run", c => { for (int i = 0; i < c; i++) ObjectSink = sampler.Run(xt, draws: 400, burn: 100, chains: 1, seed: 3); }, 500, "iteration"),
            });
        }
    }

    private static void VectorExp(double[] x, double[] y, int count)
    {
        int i = 0;
        for (; i + 3 < count; i += 4)
            Vector256.Exp(Vector256.LoadUnsafe(ref x[0], (nuint)i)).StoreUnsafe(ref y[0], (nuint)i);
        for (; i < count; i++) y[i] = Math.Exp(x[i]);
    }

    // ================================================================ Changepoint

    private static double LogSumExp(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double max = Math.Max(a, b);
        return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
    }

    private static double LogBinom(int n, int k) => SpecialFunctions.LogGamma(n + 1) - SpecialFunctions.LogGamma(k + 1) - SpecialFunctions.LogGamma(n - k + 1);

    private sealed class ChangepointTables
    {
        public double[,] LogD = null!, Mean = null!, Var = null!, Fwd = null!, Bwd = null!;
        public int N, K;
    }

    // ExactChangepoint.Infer's loops, verbatim; the counter tallies the LogSumExp cells each executes.
    private static void Forward(ChangepointTables t, ref long cells)
    {
        int n = t.N, k = t.K;
        for (int s = 0; s <= k; s++) for (int b = 0; b <= n; b++) t.Fwd[s, b] = double.NegativeInfinity;
        t.Fwd[0, 0] = 0.0;
        for (int b = 1; b <= n; b++) t.Fwd[1, b] = t.LogD[1, b];
        for (int s = 2; s <= k; s++)
            for (int b = s; b <= n; b++)
            {
                double acc = double.NegativeInfinity;
                for (int a = s - 1; a < b; a++) { acc = LogSumExp(acc, t.Fwd[s - 1, a] + t.LogD[a + 1, b]); cells++; }
                t.Fwd[s, b] = acc;
            }
    }

    private static void Backward(ChangepointTables t, ref long cells)
    {
        int n = t.N, k = t.K;
        for (int s = 0; s <= k; s++) for (int b = 0; b <= n; b++) t.Bwd[s, b] = double.NegativeInfinity;
        t.Bwd[0, n] = 0.0;
        for (int a = 0; a < n; a++) t.Bwd[1, a] = t.LogD[a + 1, n];
        for (int s = 2; s <= k; s++)
            for (int a = 0; a <= n - s; a++)
            {
                double acc = double.NegativeInfinity;
                for (int c = a + 1; c <= n - s + 1; c++) { acc = LogSumExp(acc, t.LogD[a + 1, c] + t.Bwd[s - 1, c]); cells++; }
                t.Bwd[s, a] = acc;
            }
    }

    private static void Breaks(ChangepointTables t, double logEvidence, double[] breaks, ref long cells)
    {
        int n = t.N, k = t.K;
        double logPk = -Math.Log(k);
        for (int i = 1; i <= n - 1; i++)
        {
            double acc = double.NegativeInfinity;
            for (int s = 2; s <= k; s++)
            {
                double pre = logPk - LogBinom(n - 1, s - 1);
                int jMax = Math.Min(i, s - 1);
                for (int j = 1; j <= jMax; j++)
                {
                    int after = s - j;
                    if (n - i < after) continue;
                    acc = LogSumExp(acc, pre + t.Fwd[j, i] + t.Bwd[after, i]);
                    cells++;
                }
            }
            breaks[i - 1] = Math.Exp(acc - logEvidence);
        }
    }

    private static double CurveLoop(ChangepointTables t, double logEvidence, double[] diffMean, double[] diffM2, ref long cells)
    {
        int n = t.N, k = t.K;
        double logPk = -Math.Log(k), total = 0.0;
        for (int a = 1; a <= n; a++)
            for (int b = a; b <= n; b++)
            {
                double logPhi = double.NegativeInfinity;
                for (int j = 0; j < k; j++)
                {
                    double fj = t.Fwd[j, a - 1];
                    if (double.IsNegativeInfinity(fj)) continue;
                    for (int m = 0; j + 1 + m <= k; m++)
                    {
                        double bm = t.Bwd[m, b];
                        if (double.IsNegativeInfinity(bm)) continue;
                        logPhi = LogSumExp(logPhi, logPk - LogBinom(n - 1, j + m) + fj + bm);
                        cells++;
                    }
                }
                if (double.IsNegativeInfinity(logPhi)) continue;
                double w = Math.Exp(t.LogD[a, b] + logPhi - logEvidence);
                if (w <= 0.0) continue;
                double mab = t.Mean[a, b];
                diffMean[a] += w * mab; diffMean[b + 1] -= w * mab;
                double m2 = w * (t.Var[a, b] + mab * mab);
                diffM2[a] += m2; diffM2[b + 1] -= m2;
                total += w;
            }
        return total;
    }

    internal static void ChangepointBlock()
    {
        const string block = "engine:changepoint";
        foreach ((int n, int k) in new[] { (100, 5), (100, 20), (300, 5), (300, 10), (300, 20), (1000, 5), (1000, 10), (1000, 20) })
        {
            var rng = new Xoshiro256PlusPlus(3 + n);
            double[] levels = [0.0, 2.0, -1.0, 1.0];
            var y = new double[n];
            for (int i = 0; i < n; i++) y[i] = levels[i * 4 / n] + 0.5 * Gaussian(rng);
            (double sigma2, double nu, double rho2) = ExactChangepoint.DefaultHyperparameters(y);
            var model = new GaussianSegmentModel(y, sigma2, nu, rho2);
            var t = new ChangepointTables
            {
                N = n, K = k, LogD = new double[n + 1, n + 1], Mean = new double[n + 1, n + 1], Var = new double[n + 1, n + 1],
                Fwd = new double[k + 1, n + 1], Bwd = new double[k + 1, n + 1],
            };
            model.Fill(t.LogD, t.Mean, t.Var);
            long fwdCells = 0, bwdCells = 0, breakCells = 0, curveCells = 0;
            Forward(t, ref fwdCells);
            Backward(t, ref bwdCells);
            double logPk = -Math.Log(k), logEvidence = double.NegativeInfinity;
            for (int s = 1; s <= k; s++) logEvidence = LogSumExp(logEvidence, logPk - LogBinom(n - 1, s - 1) + t.Fwd[s, n]);
            var breaks = new double[n - 1];
            var diffMean = new double[n + 2];
            var diffM2 = new double[n + 2];
            Breaks(t, logEvidence, breaks, ref breakCells);
            CurveLoop(t, logEvidence, diffMean, diffM2, ref curveCells);
            ChangepointResult result = ExactChangepoint.Infer(y, model, k);
            if (Math.Abs(result.LogEvidence - logEvidence) > 1e-9 * Math.Max(1.0, Math.Abs(logEvidence)))
                throw new InvalidOperationException($"The replica's evidence {logEvidence:R} differs from Infer's {result.LogEvidence:R} at n = {n}, k = {k}.");
            string name = $"n={n} k={k}";
            CostProfile.Row(block, name, "(count)", "forward_cells", fwdCells);
            CostProfile.Row(block, name, "(count)", "backward_cells", bwdCells);
            CostProfile.Row(block, name, "(count)", "break_cells", breakCells);
            CostProfile.Row(block, name, "(count)", "curve_cells", curveCells);
            CostProfile.Row(block, name, "(count)", "table_bytes", 8.0 * 4 * (n + 1) * (n + 1));
            long sink = 0;
            var parts = new List<CostProfile.Part>
            {
                new("infer_whole", c => { for (int i = 0; i < c; i++) ObjectSink = ExactChangepoint.Infer(y, model, k); }, 1, "inference"),
                new("model_fill", c => { for (int i = 0; i < c; i++) model.Fill(t.LogD, t.Mean, t.Var); }, 0.5 * n * n, "segment"),
                new("forward_dp", c => { for (int i = 0; i < c; i++) Forward(t, ref sink); }, fwdCells, "LogSumExp cell"),
                new("backward_dp", c => { for (int i = 0; i < c; i++) Backward(t, ref sink); }, bwdCells, "LogSumExp cell"),
                new("break_probabilities", c => { for (int i = 0; i < c; i++) Breaks(t, logEvidence, breaks, ref sink); }, Math.Max(1, breakCells), "LogSumExp cell"),
                new("regression_curve_loop", c => { for (int i = 0; i < c; i++) { Array.Clear(diffMean); Array.Clear(diffM2); Sink = CurveLoop(t, logEvidence, diffMean, diffM2, ref sink); } }, Math.Max(1, curveCells), "LogSumExp cell"),
            };
            CostProfile.Point(block, name, parts);
            ChangepointPosterior posterior = result.Posterior;
            CostProfile.Point(block, name + " draws", new List<CostProfile.Part>
            {
                new("sample_1000_draws", c => { for (int i = 0; i < c; i++) ObjectSink = posterior.Sample(1000, 7); }, 1000, "draw"),
            });
        }
    }

    // ================================================================ The collapsed g-prior block

    internal static void GPriorBlock()
    {
        foreach (int m in new[] { 5, 20, 50, 100, 200 })
        {
            double[,] g = CostKernels.DenseSpd(m, 800 + m);
            var factor = new DenseCholesky(m);
            factor.Decompose(g);
            double[] beta = CostKernels.Random(m, 21);
            const int rows = 1000;
            CollapsedGPrior block = CollapsedGPrior.FromLeastSquares(factor, beta, 5.0, 40.0, rows, null);
            InverseGamma prior = Sigma2Prior.ScaledToResponse().For(1.0) ?? throw new InvalidOperationException("No scale.");
            double tau = 1.0 / rows;
            GPriorConditional conditional = block.At(tau);
            var draw = new double[m];
            var rng = new Xoshiro256PlusPlus(8);
            CostProfile.Point("engine:gprior", $"M={m}", new List<CostProfile.Part>
            {
                new("from_least_squares", c => { for (int i = 0; i < c; i++) ObjectSink = CollapsedGPrior.FromLeastSquares(factor, beta, 5.0, 40.0, rows, null); }, 1, "block"),
                new("marginal_at", c => { for (int i = 0; i < c; i++) Sink = block.MarginalAt(tau).S; }, 1, "call"),
                new("marginal_at_log_marginal", c => { for (int i = 0; i < c; i++) Sink = block.MarginalAt(tau).LogMarginal(prior); }, 1, "call"),
                new("conditional_at", c => { for (int i = 0; i < c; i++) ObjectSink = block.At(tau); }, m, "coefficient"),
                new("conditional_sample", c => { for (int i = 0; i < c; i++) Sink = conditional.Sample(prior, rng, draw); }, (double)m * m, "term"),
            });
        }
    }
}
