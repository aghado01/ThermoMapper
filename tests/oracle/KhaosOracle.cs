using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Rng;

namespace Maths.Oracle.Tests;

/// <summary>
/// What the KHAOS parity classes share (E19's build brief, slice 0): the references computed on this side, each sharing nothing with
/// khaos's code — the standardized shifted Legendre recurrence, compositions, random_partition's law by counting, the Poisson-binomial
/// by enumeration, |𝒜| by enumeration and by its sum — and the collapsed marginal through the engine's <see cref="PriorPencil"/>, the
/// prior-precision route E21's builder takes. The modes run against CRAN's khaos 2.1.1, pinned in <c>r/renv.lock</c>, through
/// <see cref="GbassOracle"/>'s runner; a class returns only when R itself is absent.
/// </summary>
internal static class KhaosOracle
{
    /// <summary>Rows padded with NaN to the longest, the modes' ragged-case form.</summary>
    public static double[][] Padded(IEnumerable<double[]> rows)
    {
        double[][] all = [.. rows];
        int width = all.Max(r => r.Length);
        return [.. all.Select(r => r.Concat(Enumerable.Repeat(double.NaN, width - r.Length)).ToArray())];
    }

    /// <summary>ψⱼ(x) = √(2j + 1)·Pⱼ(2x − 1), the standardized shifted Legendre polynomial, by Bonnet's recurrence.</summary>
    public static double Psi(int j, double x)
    {
        double t = 2.0 * x - 1.0, previous = 1.0, current = t;
        if (j == 0) return 1.0;
        for (int k = 1; k < j; k++)
        {
            double next = ((2 * k + 1) * t * current - k * previous) / (k + 1);
            previous = current;
            current = next;
        }
        return Math.Sqrt(2 * j + 1) * current;
    }

    /// <summary>Every composition of <paramref name="d"/> into <paramref name="q"/> positive parts, lexicographic.</summary>
    public static List<int[]> Compositions(int d, int q)
    {
        var result = new List<int[]>();
        if (q == 1) { result.Add([d]); return result; }
        for (int a = 1; a <= d - q + 1; a++)
            foreach (int[] rest in Compositions(d - a, q - 1)) result.Add([a, .. rest]);
        return result;
    }

    /// <summary>
    /// The probability of composition <paramref name="c"/> under q − 1 cuts drawn uniformly with replacement from {0, …, d − q} and
    /// sorted, the scheme khaos's random_partition writes: the cuts are the partial sums of c − 1, and a multiset of cuts with
    /// multiplicities m arises from (q − 1)!/Πm! of the (d − q + 1)^(q − 1) ordered draws.
    /// </summary>
    public static double SortedCutsProbability(int[] c)
    {
        int q = c.Length, r = c.Sum() - q;
        var cuts = new List<int>();
        int sum = 0;
        for (int k = 0; k < q - 1; k++) { sum += c[k] - 1; cuts.Add(sum); }
        double log = LogFactorial(q - 1) - (q - 1) * Math.Log(r + 1);
        foreach (var group in cuts.GroupBy(v => v)) log -= LogFactorial(group.Count());
        return Math.Exp(log);
    }

    public static double LogFactorial(int k) => SpecialFunctions.LogGamma(k + 1.0);

    public static double Choose(int n, int k) => k < 0 || k > n ? 0.0 : Math.Round(Math.Exp(LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k)));

    /// <summary>|𝒜_{p,d,q}| by its sum Σ_{i≤q} Σ_{j≤d} C(p, i)·C(j − 1, i − 1), in exact integers.</summary>
    public static long CountBySum(int p, int d, int q)
    {
        long total = 0;
        for (int i = 1; i <= q; i++)
            for (int j = 1; j <= d; j++) total += Binomial(p, i) * Binomial(j - 1, i - 1);
        return total;
    }

    private static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long b = 1;
        for (int i = 1; i <= k; i++) b = b * (n - k + i) / i;
        return b;
    }

    /// <summary>|𝒜_{p,d,q}| by walking the multi-indices: 1 ≤ Σα ≤ d, at most q nonzero.</summary>
    public static long CountByEnumeration(int p, int d, int q)
    {
        long Walk(int index, int sum, int nonzero)
        {
            if (index == p) return sum >= 1 ? 1 : 0;
            long total = 0;
            for (int v = 0; sum + v <= d; v++)
            {
                int nz = nonzero + (v > 0 ? 1 : 0);
                if (nz > q) break;
                total += Walk(index + 1, sum + v, nz);
            }
            return total;
        }
        return Walk(0, 0, 0);
    }

    /// <summary>P(χ | w) under independent coins.</summary>
    public static double CoinProbability(IReadOnlyList<int> chi, IReadOnlyList<double> w)
    {
        double p = 1.0;
        for (int j = 0; j < w.Count; j++) p *= chi[j] == 1 ? w[j] : 1.0 - w[j];
        return p;
    }

    /// <summary>Every inclusion vector over p coins.</summary>
    public static IEnumerable<int[]> AllInclusions(int p)
    {
        for (int mask = 0; mask < 1 << p; mask++) yield return [.. Enumerable.Range(0, p).Select(j => (mask >> j) & 1)];
    }

    /// <summary>The Poisson-binomial pmf of Σχ by enumeration over the 2^p vectors.</summary>
    public static double[] PoissonBinomial(IReadOnlyList<double> w)
    {
        var pmf = new double[w.Count + 1];
        foreach (int[] chi in AllInclusions(w.Count)) pmf[chi.Sum()] += CoinProbability(chi, w);
        return pmf;
    }

    /// <summary>
    /// ln p(y | B, S₀) under β | σ² ~ N(0, σ²S₀), σ² ~ IG(a, b), through <see cref="PriorPencil"/> with likelihood Gram A = BᵀB, prior
    /// precision D = S₀⁻¹/τ and τ = <paramref name="tau"/>: the pencil's terms give ½ln|τD| − ½ln|A + τD| and the shrinkage, the least
    /// squares fit the rest of yᵀy − cᵀ(A + τD)⁻¹c.
    /// </summary>
    public static double PencilMarginal(double[][] basis, double[] y, double[,] priorPrecision, double tau, double a, double b)
    {
        int n = y.Length, p = basis[0].Length;
        var gram = new double[p, p];
        var cross = new double[p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
            {
                cross[j] += basis[i][j] * y[i];
                for (int k = 0; k < p; k++) gram[j, k] += basis[i][j] * basis[i][k];
            }
        PriorPencil pencil = PriorPencil.Of(gram, priorPrecision, cross, p - 1, interceptInside: true)
                             ?? throw new InvalidOperationException("The pencil is not positive definite.");
        RidgeTerms terms = pencil.At(tau);
        var chol = new DenseCholesky(p);
        chol.Decompose(gram);
        double rss = y.Sum(v => v * v) - chol.InverseQuadraticForm(cross);
        return SpecialFunctions.LogGamma(a + n / 2.0) - SpecialFunctions.LogGamma(a) + a * Math.Log(b) - n / 2.0 * Math.Log(2 * Math.PI)
               + terms.LogDeterminant - (a + n / 2.0) * Math.Log(b + 0.5 * (rss + terms.Shrinkage));
    }

    /// <summary>D_g⁻¹·BᵀB·D_g⁻¹, the modified g-prior's precision per unit 1/g₀² (Rumsey, Francom et al., arXiv 2510.25036).</summary>
    public static double[,] GPriorPrecision(double[][] basis, double[] g)
    {
        int p = g.Length;
        var d = new double[p, p];
        foreach (double[] row in basis)
            for (int j = 0; j < p; j++)
                for (int k = 0; k < p; k++) d[j, k] += row[j] * row[k] / (g[j] * g[k]);
        return d;
    }

    /// <summary>The identity, the raw ridge's precision per unit 1/τ².</summary>
    public static double[,] Identity(int p)
    {
        var d = new double[p, p];
        for (int j = 0; j < p; j++) d[j, j] = 1.0;
        return d;
    }

    /// <summary>The first <paramref name="columns"/> columns of each row.</summary>
    public static double[][] Columns(double[][] basis, int columns) => [.. basis.Select(r => r[..columns])];

    /// <summary>ln IG(x; shape, rate)'s kernel, −(shape + 1)ln x − rate/x.</summary>
    public static double LogInverseGammaKernel(double x, double shape, double rate) => -(shape + 1.0) * Math.Log(x) - rate / x;

    /// <summary>The spread of a − b over paired values: zero exactly when they differ by a constant.</summary>
    public static double Spread(IEnumerable<double> differences)
    {
        double[] d = [.. differences];
        return d.Max() - d.Min();
    }

    // ---- the conditionals mode's fixture, shared by the conditional and defect classes ----------------------------------------

    /// <summary>The frozen state's terms: each column's variables and degrees (inputs 0-based), after the intercept.</summary>
    public static readonly (int[] Vars, int[] Degrees)[] StateTerms = [([0], [1]), ([1], [2]), ([0, 2], [1, 2]), ([1], [3])];

    public const int StateRows = 60;
    public static readonly double[] StateBeta = [1.0, 0.8, -0.5, 0.3, 0.05];
    public static readonly double[] Theta = [0.5, 2, 10, 40, 200, 1000, 5000];
    public static readonly double[] Sigma2Grid = [0.03, 0.06, 0.09, 0.15, 0.3];
    public const double Zeta = 1, ASigma = 2, BSigma = 0.5, AG = 2, BG = 50, G0Sq = 40, Tau2 = 100, S2 = 0.09;

    /// <summary>The frozen state: sixty cases on [0, 1]³, the basis (intercept first) and y = f + 0.3·N(0, 1).</summary>
    public static (double[][] Basis, double[] Y) State()
    {
        var rng = new Xoshiro256PlusPlus(20261107);
        var basis = new double[StateRows][];
        var y = new double[StateRows];
        for (int i = 0; i < StateRows; i++)
        {
            double[] x = [rng.NextDouble(), rng.NextDouble(), rng.NextDouble()];
            basis[i] = [1.0, .. StateTerms.Select(t => t.Vars.Select((v, k) => Psi(t.Degrees[k], x[v])).Aggregate(1.0, (u, w) => u * w))];
            y[i] = 1.0 + 0.8 * basis[i][1] - 0.5 * basis[i][2] + 0.3 * basis[i][3] + 0.3 * Normal.SampleStandard(rng);
        }
        return (basis, y);
    }

    /// <summary>Each column's g_m = (1 + q(d + q − 2))^(−ζ/2), the intercept's 1.</summary>
    public static double[] StateG() =>
        [1.0, .. StateTerms.Select(t => { int q = t.Vars.Length, d = t.Degrees.Sum(); return Math.Pow(1.0 + q * (d + q - 2.0), -Zeta / 2.0); })];

    /// <summary>The proposal cases: births from two states (order below p, and order = p), two degree mutations, a variable swap,
    /// and the prior ratios with the move probabilities at M = 0.</summary>
    public static readonly double[][] ProposalCases =
    [
        [1, 3, 6, 3, 0, 1, 0],
        [1, 3, 5, 0, 2, 1],
        [2, 6, 0, 2, 3, 5],
        [2, 6, 1, 3, 4, 6],
        [3, 4, 2, 1, 3, 2, 1, 1, 0],
        [4, 4, 2.5, 30],
    ];

    public static readonly Lazy<JsonElement?> Conditionals = new(() =>
    {
        if (!ROracle.IsAvailable) return null;
        (double[][] basis, double[] y) = State();
        string proposals = GbassOracle.WriteCsv(Padded(ProposalCases));
        static string Join(params IEnumerable<double> v) => string.Join(";", v.Select(x => x.ToString("R", System.Globalization.CultureInfo.InvariantCulture)));
        try
        {
            return GbassOracle.Run("khaos_conditionals_oracle.R", basis.Select((r, i) => new[] { y[i] }.Concat(r).ToArray()), proposals,
                $"q={Join(StateTerms.Select(t => (double)t.Vars.Length))}", $"d={Join(StateTerms.Select(t => (double)t.Degrees.Sum()))}",
                $"zeta={Join(Zeta)}", $"a_sigma={Join(ASigma)}", $"b_sigma={Join(BSigma)}", $"a_g={Join(AG)}", $"b_g={Join(BG)}",
                $"g0sq={Join(G0Sq)}", $"tau2={Join(Tau2)}", $"beta={Join(StateBeta)}", $"s2={Join(S2)}", $"theta={Join(Theta)}",
                $"s2grid={Join(Sigma2Grid)}", "draws=20000", "seed=7");
        }
        finally { System.IO.File.Delete(proposals); }
    });
}
