using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Intrinsics;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Bgp;
using Maths.Regression.Changepoint;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;

namespace Maths.Regression.Tests;

/// <summary>
/// The levers priced by test-side replicas (the brief's <c>engine:levers</c> block): each engine's step with the levers the parts
/// show applied together, checked against production in the same process (bit identity where the arithmetic and its order are
/// unchanged, the relative difference where a reduction or a vector <c>exp</c> changes them), and timed beside it. Nothing here is
/// production code; the replicas are what a change under <c>src/</c> would have to match.
/// </summary>
internal static class CostLevers
{
    private static double Sink { get => CostKernels.Sink; set => CostKernels.Sink = value; }
    private static object? ObjectSink { set => CostKernels.ObjectSink = value; }

    // ================================================================ BARS: the evaluation

    /// <summary>
    /// <see cref="SplineTarget.Evaluate"/> for the unweighted Normal model, with: the basis row computed without allocation straight
    /// into a compact n × (degree + 1) design; the support gate, the Gram band and the right-hand side read from it (n·(p + 1), not
    /// n·ν); the response's centre and σ² scale taken once; the partition-of-unity scan of the centring left out (a B-spline basis
    /// has it by construction); scratch arrays reused. The factor, its admission, the solve and the collapsed block are production's.
    /// </summary>
    internal sealed class BarsReplica
    {
        private const int Degree = 3, Width = Degree + 1;
        private readonly CostEngines.BarsFixture _fx;
        private readonly double[] _x, _y;
        private readonly int _n;
        private readonly double _level;
        private readonly InverseGamma? _sigma2;
        private readonly double[] _z, _response, _fit;
        private readonly int[] _lo;

        public BarsReplica(CostEngines.BarsFixture fx)
        {
            _fx = fx;
            _x = fx.X; _y = fx.Y; _n = _x.Length;
            (_sigma2, _level) = fx.Model.ResolvedPriors(_y);
            _z = new double[_n * Width];
            _lo = new int[_n];
            _response = new double[_n];
            _fit = new double[_n];
            for (int i = 0; i < _n; i++) _response[i] = _y[i] - _level;
        }

        // Cox–de Boor / Piegl & Tiller A2.2, as BSpline.EvaluateBasis, on the caller's scratch.
        private static int Row(double t, double[] knots, int n, Span<double> local, Span<double> left, Span<double> right)
        {
            int span = Degree;
            while (span < n - 1 && t >= knots[span + 1]) span++;
            local[0] = 1.0;
            for (int j = 1; j <= Degree; j++)
            {
                left[j] = t - knots[span + 1 - j];
                right[j] = knots[span + j] - t;
                double saved = 0.0;
                for (int r = 0; r < j; r++)
                {
                    double denom = right[r + 1] + left[j - r];
                    double temp = denom > 0.0 ? local[r] / denom : 0.0;
                    local[r] = saved + right[r + 1] * temp;
                    saved = left[j - r] * temp;
                }
                local[j] = saved;
            }
            return span - Degree;
        }

        /// <summary>ln p(y | knots) at b = 1, or −∞ where production refuses the design.</summary>
        public double LogMarginal(KnotConfig config)
        {
            double[] knots = BSpline.MakeClampedKnots(config.InteriorKnots, Degree);
            int nu = knots.Length - Degree - 1;
            if (_n <= nu) return double.NegativeInfinity;
            Span<double> local = stackalloc double[Width], left = stackalloc double[Width], right = stackalloc double[Width];
            for (int i = 0; i < _n; i++)
            {
                local.Clear();
                int lo = Row(_x[i], knots, nu, local, left, right);
                _lo[i] = lo;
                local.CopyTo(_z.AsSpan(i * Width, Width));
            }
            Span<double> mass = nu <= 256 ? stackalloc double[nu] : new double[nu];
            mass.Clear();
            for (int i = 0; i < _n; i++)
                for (int r = 0; r < Width; r++) mass[_lo[i] + r] += _z[i * Width + r];
            for (int j = 0; j < nu; j++) if (mass[j] < 1.0) return double.NegativeInfinity;

            int bw = Degree;
            var band = new double[bw + 1, nu];
            var b = new double[nu];
            for (int i = 0; i < _n; i++)
            {
                int lo = _lo[i];
                double ri = _y[i];
                for (int p = 0; p < Width; p++)
                {
                    double wzip = 1.0 * _z[i * Width + p];
                    b[lo + p] += wzip * ri;
                    for (int q = 0; q <= p; q++) band[p - q, lo + q] += wzip * _z[i * Width + q];
                }
            }
            var gram = new BandCholesky(nu, bw, BandFactorization.Ldlt, strict: true);
            gram.DecomposeBanded(band);
            if (!DesignAdmission.Admits(gram)) return double.NegativeInfinity;

            var rhs = new double[nu];
            for (int i = 0; i < _n; i++)
            {
                double wr = 1.0 * _response[i];
                int lo = _lo[i];
                for (int p = 0; p < Width; p++) rhs[lo + p] += _z[i * Width + p] * wr;
            }
            double[] beta = gram.Solve(rhs);
            double rss = 0.0, fitted = 0.0;
            for (int i = 0; i < _n; i++)
            {
                double s = 0.0;
                int lo = _lo[i];
                for (int p = 0; p < Width; p++) s += _z[i * Width + p] * beta[lo + p];
                _fit[i] = s;
                double r = _response[i] - s;
                rss += r * r;
                fitted += s * s;
            }
            var priorMean = new double[nu];
            Array.Fill(priorMean, _level);
            CollapsedGPrior block = CollapsedGPrior.FromLeastSquares(gram, beta, rss, fitted, _n, priorMean);
            return block.MarginalAt(1.0 / _n).LogMarginal(_sigma2);
        }
    }

    // ================================================================ BGP: the fit

    internal sealed class GpReplica
    {
        private readonly int _n;
        private readonly double[,] _sq;
        private readonly double[] _y, _ones;
        private readonly double _nugget;
        private readonly SquaredExponentialKernel _kernel = new();
        private readonly double[,] _a;
        private readonly DenseCholesky _factor;
        private readonly double[] _row, _rowOut, _kInvY, _kInvOne;

        public GpReplica(double[,] x, double[] y, double nugget)
        {
            _n = x.GetLength(0);
            _y = y; _nugget = nugget;
            _sq = new double[_n, _n];
            for (int i = 0; i < _n; i++)
                for (int j = 0; j < i; j++)
                {
                    double s = 0;
                    for (int c = 0; c < x.GetLength(1); c++) { double d = (x[i, c] - x[j, c]) / 12.0; s += d * d; }
                    _sq[i, j] = s;
                }
            _ones = Enumerable.Repeat(1.0, _n).ToArray();
            _a = new double[_n, _n];
            _factor = new DenseCholesky(_n);
            _row = new double[_n]; _rowOut = new double[_n]; _kInvY = new double[_n]; _kInvOne = new double[_n];
        }

        private double Reduce(double logDet, double[] kInvY, double[] kInvOne)
        {
            double oneKOne = 0.0, oneKY = 0.0;
            for (int i = 0; i < _n; i++) { oneKOne += kInvOne[i]; oneKY += kInvY[i]; }
            double mHat = oneKY / oneKOne, quad = 0.0;
            for (int i = 0; i < _n; i++) quad += (_y[i] - mHat) * (kInvY[i] - mHat * kInvOne[i]);
            return -0.5 * quad - 0.5 * logDet - 0.5 * Math.Log(oneKOne) - 0.5 * (_n - 1) * Math.Log(2.0 * Math.PI);
        }

        /// <summary>The fit as production writes it: fresh arrays, the interface kernel, <see cref="CholeskyDecomposition"/> (with L⁻¹), two solves.</summary>
        public double Verbatim(double t)
        {
            var a = new double[_n, _n];
            for (int i = 0; i < _n; i++)
            {
                a[i, i] = _kernel.Evaluate(0.0, t) + _nugget;
                for (int j = 0; j < i; j++) { double v = _kernel.Evaluate(_sq[i, j], t); a[i, j] = v; a[j, i] = v; }
            }
            var chol = new CholeskyDecomposition(_n);
            chol.Decompose(a);
            return Reduce(chol.LogDet, chol.Solve(_y), chol.Solve(_ones));
        }

        /// <summary>Scratch reused, the kernel inline, the strict factor without L⁻¹ (<see cref="DenseCholesky"/>), the solves in place.</summary>
        public double Lean(double t, bool vectorExp)
        {
            double[,] a = _a;
            double inv = 1.0 / (2.0 * t);
            for (int i = 0; i < _n; i++)
            {
                a[i, i] = 1.0 + _nugget;
                if (vectorExp)
                {
                    for (int j = 0; j < i; j++) _row[j] = -_sq[i, j] / (2.0 * t);
                    int c = 0;
                    for (; c + 3 < i; c += 4) Vector256.Exp(Vector256.LoadUnsafe(ref _row[0], (nuint)c)).StoreUnsafe(ref _rowOut[0], (nuint)c);
                    for (; c < i; c++) _rowOut[c] = Math.Exp(_row[c]);
                    for (int j = 0; j < i; j++) { a[i, j] = _rowOut[j]; a[j, i] = _rowOut[j]; }
                }
                else
                    for (int j = 0; j < i; j++) { double v = Math.Exp(-_sq[i, j] / (2.0 * t)); a[i, j] = v; a[j, i] = v; }
            }
            _ = inv;
            _factor.Decompose(a);
            Array.Copy(_y, _kInvY, _n);
            _factor.SolveInPlace(_kInvY);
            Array.Copy(_ones, _kInvOne, _n);
            _factor.SolveInPlace(_kInvOne);
            return Reduce(_factor.LogDet, _kInvY, _kInvOne);
        }
    }

    // ================================================================ Changepoint: the curve loop

    internal static double CurveHoisted(CostEngines.ChangepointTables t, double logEvidence, double[] diffMean, double[] diffM2, double[] logBinom, bool twoPass, bool vector, double[] buffer)
    {
        int n = t.N, k = t.K;
        double logPk = -Math.Log(k), total = 0.0;
        for (int a = 1; a <= n; a++)
            for (int b = a; b <= n; b++)
            {
                double logPhi = double.NegativeInfinity;
                if (!twoPass)
                {
                    for (int j = 0; j < k; j++)
                    {
                        double fj = t.Fwd[j, a - 1];
                        if (double.IsNegativeInfinity(fj)) continue;
                        for (int m = 0; j + 1 + m <= k; m++)
                        {
                            double bm = t.Bwd[m, b];
                            if (double.IsNegativeInfinity(bm)) continue;
                            logPhi = CostEngines.LogSumExp(logPhi, logPk - logBinom[j + m] + fj + bm);
                        }
                    }
                }
                else
                {
                    int count = 0;
                    double max = double.NegativeInfinity;
                    for (int j = 0; j < k; j++)
                    {
                        double fj = t.Fwd[j, a - 1];
                        if (double.IsNegativeInfinity(fj)) continue;
                        for (int m = 0; j + 1 + m <= k; m++)
                        {
                            double bm = t.Bwd[m, b];
                            if (double.IsNegativeInfinity(bm)) continue;
                            double term = logPk - logBinom[j + m] + fj + bm;
                            buffer[count++] = term;
                            if (term > max) max = term;
                        }
                    }
                    if (count > 0)
                    {
                        for (int i = 0; i < count; i++) buffer[i] -= max;
                        double sum = 0.0;
                        if (vector)
                        {
                            int i = 0;
                            for (; i + 3 < count; i += 4)
                                sum += Vector256.Sum(Vector256.Exp(Vector256.LoadUnsafe(ref buffer[0], (nuint)i)));
                            for (; i < count; i++) sum += Math.Exp(buffer[i]);
                        }
                        else
                            for (int i = 0; i < count; i++) sum += Math.Exp(buffer[i]);
                        logPhi = max + Math.Log(sum);
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

    // ================================================================ The block

    internal static void LeversBlock()
    {
        const string block = "engine:levers";

        foreach (int n in new[] { 100, 500, 2000 })
            foreach (int k in new[] { 3, 8, 16 })
            {
                var fx = new CostEngines.BarsFixture(n);
                var replica = new BarsReplica(fx);
                KnotConfig cfg = CostEngines.BarsFixture.Equispaced(k);
                double production = fx.Target.Evaluate(cfg).LogLikelihood(1.0);
                double mine = replica.LogMarginal(cfg);
                string name = $"BARS evaluation n={n} k={k}";
                CostProfile.Row(block, name, "(check)", "log_marginal_bits_identical", BitConverter.DoubleToInt64Bits(production) == BitConverter.DoubleToInt64Bits(mine) ? 1 : 0);
                CostProfile.Row(block, name, "(check)", "log_marginal_relative_difference", Math.Abs(production - mine) / Math.Max(1e-300, Math.Abs(production)));
                CostProfile.Point(block, name, new List<CostProfile.Part>
                {
                    new("production_evaluate", c => { for (int i = 0; i < c; i++) ObjectSink = fx.Target.Evaluate(cfg); }, 1, "evaluation"),
                    new("replica_levers_applied", c => { for (int i = 0; i < c; i++) Sink = replica.LogMarginal(cfg); }, 1, "evaluation"),
                });
            }

        foreach (int n in new[] { 100, 400, 800 })
        {
            (double[,] x, double[] y) = CostEngines.SwissRoll(n, 5 + n);
            var gp = new GpReplica(x, y, 0.01);
            const double t = 0.05;
            double v0 = gp.Verbatim(t), v1 = gp.Lean(t, vectorExp: false), v2 = gp.Lean(t, vectorExp: true);
            string name = $"BGP fit n={n}";
            CostProfile.Row(block, name, "(check)", "verbatim_value", v0);
            CostProfile.Row(block, name, "(check)", "lean_value", v1);
            CostProfile.Row(block, name, "(check)", "lean_relative_difference", Math.Abs(v1 - v0) / Math.Abs(v0));
            CostProfile.Row(block, name, "(check)", "lean_bits_identical", BitConverter.DoubleToInt64Bits(v1) == BitConverter.DoubleToInt64Bits(v0) ? 1 : 0);
            CostProfile.Row(block, name, "(check)", "lean_vector_exp_relative_difference", Math.Abs(v2 - v0) / Math.Abs(v0));
            CostProfile.Point(block, name, new List<CostProfile.Part>
            {
                new("replica_verbatim", c => { for (int i = 0; i < c; i++) Sink = gp.Verbatim(t + 1e-6 * i); }, 1, "fit"),
                new("replica_lean_scalar_exp", c => { for (int i = 0; i < c; i++) Sink = gp.Lean(t + 1e-6 * i, vectorExp: false); }, 1, "fit"),
                new("replica_lean_vector_exp", c => { for (int i = 0; i < c; i++) Sink = gp.Lean(t + 1e-6 * i, vectorExp: true); }, 1, "fit"),
            });
        }

        foreach ((int n, int k) in new[] { (300, 10), (1000, 10) })
        {
            var rng = new Xoshiro256PlusPlus(3 + n);
            double[] levels = [0.0, 2.0, -1.0, 1.0];
            var y = new double[n];
            for (int i = 0; i < n; i++) y[i] = levels[i * 4 / n] + 0.5 * CostEngines.Gaussian(rng);
            (double sigma2, double nu, double rho2) = ExactChangepoint.DefaultHyperparameters(y);
            var model = new GaussianSegmentModel(y, sigma2, nu, rho2);
            var t = new CostEngines.ChangepointTables
            {
                N = n, K = k, LogD = new double[n + 1, n + 1], Mean = new double[n + 1, n + 1], Var = new double[n + 1, n + 1],
                Fwd = new double[k + 1, n + 1], Bwd = new double[k + 1, n + 1],
            };
            model.Fill(t.LogD, t.Mean, t.Var);
            long cells = 0;
            CostEngines.Forward(t, ref cells);
            CostEngines.Backward(t, ref cells);
            double logPk = -Math.Log(k), logEvidence = double.NegativeInfinity;
            for (int s = 1; s <= k; s++) logEvidence = CostEngines.LogSumExp(logEvidence, logPk - CostEngines.LogBinom(n - 1, s - 1) + t.Fwd[s, n]);
            double[] logBinom = Enumerable.Range(0, k + 1).Select(s => CostEngines.LogBinom(n - 1, s)).ToArray();
            var buffer = new double[(k + 1) * (k + 2) / 2 + 2];
            var diffMean = new double[n + 2];
            var diffM2 = new double[n + 2];
            long dummy = 0;
            double verbatim = CostEngines.CurveLoop(t, logEvidence, diffMean, diffM2, ref dummy);
            Array.Clear(diffMean); Array.Clear(diffM2);
            double hoisted = CurveHoisted(t, logEvidence, diffMean, diffM2, logBinom, false, false, buffer);
            double[] meanHoisted = (double[])diffMean.Clone();
            Array.Clear(diffMean); Array.Clear(diffM2);
            double twoPass = CurveHoisted(t, logEvidence, diffMean, diffM2, logBinom, true, false, buffer);
            Array.Clear(diffMean); Array.Clear(diffM2);
            double vec = CurveHoisted(t, logEvidence, diffMean, diffM2, logBinom, true, true, buffer);
            string name = $"changepoint curve loop n={n} k={k}";
            CostProfile.Row(block, name, "(check)", "verbatim_total", verbatim);
            CostProfile.Row(block, name, "(check)", "hoisted_total", hoisted);
            CostProfile.Row(block, name, "(check)", "two_pass_total", twoPass);
            CostProfile.Row(block, name, "(check)", "hoisted_total_bits_identical", BitConverter.DoubleToInt64Bits(verbatim) == BitConverter.DoubleToInt64Bits(hoisted) ? 1 : 0);
            CostProfile.Row(block, name, "(check)", "two_pass_relative_difference", Math.Abs(twoPass - verbatim) / Math.Abs(verbatim));
            CostProfile.Row(block, name, "(check)", "two_pass_vector_relative_difference", Math.Abs(vec - verbatim) / Math.Abs(verbatim));
            _ = meanHoisted;
            long sink = 0;
            CostProfile.Point(block, name, new List<CostProfile.Part>
            {
                new("verbatim_log_binom_inside", c => { for (int i = 0; i < c; i++) { Array.Clear(diffMean); Array.Clear(diffM2); Sink = CostEngines.CurveLoop(t, logEvidence, diffMean, diffM2, ref sink); } }, 1, "loop"),
                new("log_binom_hoisted", c => { for (int i = 0; i < c; i++) { Array.Clear(diffMean); Array.Clear(diffM2); Sink = CurveHoisted(t, logEvidence, diffMean, diffM2, logBinom, false, false, buffer); } }, 1, "loop"),
                new("hoisted_two_pass_max", c => { for (int i = 0; i < c; i++) { Array.Clear(diffMean); Array.Clear(diffM2); Sink = CurveHoisted(t, logEvidence, diffMean, diffM2, logBinom, true, false, buffer); } }, 1, "loop"),
                new("hoisted_two_pass_vector_exp", c => { for (int i = 0; i < c; i++) { Array.Clear(diffMean); Array.Clear(diffM2); Sink = CurveHoisted(t, logEvidence, diffMean, diffM2, logBinom, true, true, buffer); } }, 1, "loop"),
            });
        }
    }
}
