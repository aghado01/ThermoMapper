using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Regression.Conformal;

/// <summary>The settings of <see cref="ScoreMixture.Fit"/>; each default is Zhang and Candès's (arXiv 2409.19712, §6 and §6.1).</summary>
public sealed record ScoreMixtureOptions
{
    /// <summary>s, the levels ξ_t at which the score's conditional distribution function is matched: the training scores' empirical
    /// quantiles at t/(s + 1), the deciles at the default 9. Levels that repeat, or that every score reaches, are dropped.</summary>
    public int Levels { get; init; } = 9;

    /// <summary>J fixed; null chooses it by their rule: J grows from 1 until R² gains less than <see cref="ClusterGain"/>.</summary>
    public int? Clusters { get; init; }

    /// <summary>The least gain in R² from J to J + 1 that keeps J growing.</summary>
    public double ClusterGain { get; init; } = 0.05;

    /// <summary>The largest J the rule reaches.</summary>
    public int MaxClusters { get; init; } = 10;

    /// <summary>m fixed; null chooses it by their rule: the largest m in [<see cref="MinPrecision"/>, <see cref="MaxPrecision"/>], found by
    /// bisection, at which the training split's randomized weights have mean effective sample size above
    /// <see cref="TargetEffectiveSampleSize"/> and mean self-weight at most <see cref="MaxSelfWeight"/>.</summary>
    public int? Precision { get; init; }

    /// <summary>The least m the rule considers.</summary>
    public int MinPrecision { get; init; } = 5;

    /// <summary>The largest m the rule considers.</summary>
    public int MaxPrecision { get; init; } = 500;

    /// <summary>The mean effective sample size the rule's m must exceed, read on the training split: where the calibration split is smaller,
    /// the calibration's effective sample size is smaller in proportion, so scale the target by the training split's size over the
    /// calibration's (their training and calibration splits are the same size).</summary>
    public double TargetEffectiveSampleSize { get; init; } = 100.0;

    /// <summary>The mean self-weight the rule's m must not exceed.</summary>
    public double MaxSelfWeight { get; init; } = 1.0 / 30.0;

    /// <summary>The logistic fits' ridge per pooled point on the standardized features (<see cref="DensityRatio.Logistic(double[,], double[,], double, int, double)"/>).</summary>
    public double Ridge { get; init; } = 1e-4;

    /// <summary>The factorization's k-means++ starts per J; the fit with the least residual is kept.</summary>
    public int Restarts { get; init; } = 5;

    /// <summary>The seed of the k-means++ starts and of the precision rule's draws.</summary>
    public int Seed { get; init; }
}

/// <summary>What a fit cost (<see cref="ScoreMixture.Cost"/>): seconds by phase, the factorization's alternating passes, and the simplex
/// step's solves, active-set iterations and solves stopped at the iteration cap.</summary>
internal sealed class ScoreMixtureCost
{
    public double LevelSeconds, FactorSeconds, PrecisionSeconds;
    public long Passes, Solves, SolverIterations, CappedSolves;

    public override string ToString()
        => $"levels {LevelSeconds:F2} s, factorization {FactorSeconds:F2} s ({Passes} passes, {Solves} solves, {SolverIterations} active-set iterations, "
           + $"{CappedSolves} capped), precision {PrecisionSeconds:F2} s";
}

/// <summary>
/// A mixture model of a conformal score, R | X ~ Σ_k π_k(X)·F_k, fitted on a training split by Zhang and Candès's matching algorithm
/// (arXiv 2409.19712, §6; the conformal layer's slice 15b), for posterior conformal prediction's memberships π (<see cref="PosteriorConformal"/>).
/// <list type="bullet">
/// <item><b>Levels.</b> s levels ξ_t, the training scores' quantiles at t/(s + 1).</item>
/// <item><b>Fitting.</b> At each level, τ_t(x) = P(R ≤ ξ_t | x) by Bayes' rule from the density ratio g_t of the inputs at or below the level
/// against those above it: τ_t = n_t·g_t/(n − n_t + n_t·g_t), their (6.4), with g_t estimated by the logistic two-sample fit
/// (<see cref="DensityRatio"/>) on the caller's features. They fit g_t by linear least squares (Kanamori et al. 2009) because their π is
/// refitted with the test point for every candidate response, which the Sherman–Morrison–Woodbury identity makes cheap; fitted on a separate
/// training split, π is fitted once and the logistic fit, whose τ lies in (0, 1), costs nothing extra.</item>
/// <item><b>Clustering.</b> The profiles τ(xᵢ) ∈ ℝ^s are factored as τ(xᵢ) ≈ Σ_k π_ik·γ_k with each πᵢ on the simplex: γ initialized at the
/// k-means++ centres, then alternating least squares for γ given π with the exact simplex-constrained least squares for each πᵢ given γ
/// (an active-set solve, where they take mirror-descent steps). γ_k estimates cluster k's distribution function at the levels.</item>
/// <item><b>Hyperparameters</b> (their §6.1). J grows until R² gains less than 0.05; m is the largest in [5, 500] whose randomized weights over
/// the training split, each point weighed against all of them, have mean effective sample size above 100 and mean self-weight at most 1/30.</item>
/// </list>
/// The memberships at new inputs (<see cref="Memberships"/>) are the same projection: τ at the inputs, then the simplex-constrained least
/// squares against the fitted γ. Fitted on a training split apart from the calibration units and the test cases, π̂ is a fixed function of
/// the inputs, so their symmetry condition holds and the randomized weights' coverage is exact whatever its quality; the quality decides how
/// close to nominal each cluster's coverage comes. The training scores should be distributed as the calibration scores are (held out, as
/// cross-fitted residuals are).
/// </summary>
public sealed class ScoreMixture
{
    private readonly DensityRatio[] _fits;
    private readonly double[] _prior;      // ln((n − n_t)/n_t) at each level
    private readonly double[,] _profiles;  // γ, J × s
    private readonly double[,] _gram;      // γγᵀ with its ridge, J × J
    private readonly int _width;

    private ScoreMixture(DensityRatio[] fits, double[] prior, double[] levels, double[,] profiles, int width, double[] rSquared)
    {
        _fits = fits;
        _prior = prior;
        _profiles = profiles;
        _width = width;
        _gram = Gram(profiles);
        Levels = levels;
        RSquared = rSquared;
    }

    /// <summary>The levels ξ_t, ascending.</summary>
    public IReadOnlyList<double> Levels { get; }

    /// <summary>J.</summary>
    public int Clusters => _profiles.GetLength(0);

    /// <summary>γ, a row per cluster: its distribution function's estimate at each level.</summary>
    public double[,] Profiles => (double[,])_profiles.Clone();

    /// <summary>The factorization's R² on the training split at J = 1, 2, … as far as the rule went (only the fixed J's where J is fixed).</summary>
    public IReadOnlyList<double> RSquared { get; }

    /// <summary>m, fixed or chosen by the rule.</summary>
    public int Precision { get; private set; }

    /// <summary>The training split's mean effective sample size at <see cref="Precision"/>.</summary>
    public double MeanEffectiveSampleSize { get; private set; }

    /// <summary>The training split's mean self-weight at <see cref="Precision"/>.</summary>
    public double MeanSelfWeight { get; private set; }

    /// <summary>Whether <see cref="Precision"/> meets the rule's two conditions; false where even the least m does not, which then is m.</summary>
    public bool PrecisionRuleMet { get; private set; }

    /// <summary>What the fit cost.</summary>
    internal ScoreMixtureCost Cost { get; private set; } = new();

    /// <summary>
    /// Fits the mixture from the training split's <paramref name="features"/> (a row per unit, the inputs' features for the logistic fits)
    /// and <paramref name="scores"/> (its conformal scores, at least two distinct, all finite).
    /// </summary>
    public static ScoreMixture Fit(double[,] features, IReadOnlyList<double> scores, ScoreMixtureOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(features);
        ArgumentNullException.ThrowIfNull(scores);
        ScoreMixtureOptions o = Checked(options ?? new ScoreMixtureOptions());
        int n = features.GetLength(0), p = features.GetLength(1);
        if (scores.Count != n) throw new ArgumentException($"One score per row of features is required: {scores.Count} scores for {n} rows.", nameof(scores));
        for (int i = 0; i < n; i++)
            if (!double.IsFinite(scores[i])) throw new ArgumentException($"Score {i} is {scores[i]}: each is finite.", nameof(scores));
        var cost = new ScoreMixtureCost();
        var clock = System.Diagnostics.Stopwatch.StartNew();
        double[] sorted = [.. scores];
        Array.Sort(sorted);

        // The levels: the quantiles at t/(s + 1), each kept once and only where some score lies above it.
        var levels = new List<double>();
        for (int t = 1; t <= o.Levels; t++)
        {
            int at = Math.Clamp((int)Math.Ceiling(t * (double)n / (o.Levels + 1)) - 1, 0, n - 1);
            double level = sorted[at];
            if (level < sorted[n - 1] && (levels.Count == 0 || level > levels[^1])) levels.Add(level);
        }
        if (levels.Count == 0) throw new ArgumentException("The scores take one value, so no level splits them.", nameof(scores));

        int s = levels.Count;
        var fits = new DensityRatio[s];
        var prior = new double[s];
        for (int t = 0; t < s; t++)
        {
            int below = 0;
            for (int i = 0; i < n; i++) if (scores[i] <= levels[t]) below++;
            var source = new double[n - below, p];
            var target = new double[below, p];
            for (int i = 0, a = 0, b = 0; i < n; i++)
            {
                bool under = scores[i] <= levels[t];
                for (int j = 0; j < p; j++)
                {
                    if (under) target[b, j] = features[i, j];
                    else source[a, j] = features[i, j];
                }
                if (under) b++; else a++;
            }
            fits[t] = DensityRatio.Logistic(source, target, o.Ridge);
            prior[t] = Math.Log((double)(n - below) / below);
        }

        var tau = new double[n][];
        var row = new double[p];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < p; j++) row[j] = features[i, j];
            tau[i] = Tau(fits, prior, row);
        }

        cost.LevelSeconds = clock.Elapsed.TotalSeconds;
        clock.Restart();
        var rng = new Xoshiro256PlusPlus(o.Seed);
        var rSquared = new List<double>();
        double[,] chosen;
        if (o.Clusters is int fixedJ)
        {
            (double[,] gamma, double r2) = Factor(tau, fixedJ, o.Restarts, rng, cost)
                                           ?? throw new ArgumentException($"The profiles have fewer than {fixedJ} distinct points.", nameof(options));
            rSquared.Add(r2);
            chosen = gamma;
        }
        else
        {
            (chosen, double previous) = Factor(tau, 1, 1, rng, cost)!.Value;
            rSquared.Add(previous);
            for (int j = 2; j <= Math.Min(o.MaxClusters, n); j++)
            {
                (double[,] Gamma, double RSquared)? fit = Factor(tau, j, o.Restarts, rng, cost);
                if (fit is null) break;
                (double[,] gamma, double r2) = fit.Value;
                rSquared.Add(r2);
                if (r2 - previous < o.ClusterGain) break;
                chosen = gamma;
                previous = r2;
            }
        }

        cost.FactorSeconds = clock.Elapsed.TotalSeconds;
        clock.Restart();
        var mixture = new ScoreMixture(fits, prior, [.. levels], chosen, p, [.. rSquared]) { Cost = cost };
        double[][] memberships = mixture.Rows(features);
        if (o.Precision is int m)
        {
            mixture.Precision = m;
            (mixture.MeanEffectiveSampleSize, mixture.MeanSelfWeight) = Reading(memberships, m, o.Seed);
            mixture.PrecisionRuleMet = mixture.MeanEffectiveSampleSize > o.TargetEffectiveSampleSize && mixture.MeanSelfWeight <= o.MaxSelfWeight;
        }
        else
        {
            var readings = new Dictionary<int, (double Ess, double Self)>();
            bool Meets(int precision)
            {
                if (!readings.TryGetValue(precision, out var r)) readings[precision] = r = Reading(memberships, precision, o.Seed);
                return r.Ess > o.TargetEffectiveSampleSize && r.Self <= o.MaxSelfWeight;
            }
            int lo = o.MinPrecision, hi = o.MaxPrecision;
            bool met = Meets(lo);
            if (met && Meets(hi)) lo = hi;
            else if (met)
                while (hi - lo > 1)
                {
                    int mid = lo + (hi - lo) / 2;
                    if (Meets(mid)) lo = mid; else hi = mid;
                }
            mixture.Precision = lo;
            mixture.PrecisionRuleMet = met;
            (mixture.MeanEffectiveSampleSize, mixture.MeanSelfWeight) = readings[lo];
        }
        cost.PrecisionSeconds = clock.Elapsed.TotalSeconds;
        return mixture;
    }

    /// <summary>π at each row of <paramref name="features"/>: τ at the row, projected onto the simplex against γ. A row per input, a column
    /// per cluster.</summary>
    public double[,] Memberships(double[,] features)
    {
        double[][] rows = Rows(features);
        var result = new double[rows.Length, Clusters];
        for (int i = 0; i < rows.Length; i++)
            for (int k = 0; k < Clusters; k++) result[i, k] = rows[i][k];
        return result;
    }

    /// <summary>Posterior conformal prediction's weights (<see cref="PosteriorConformal.Weights"/>) at precision <see cref="Precision"/>, from
    /// the memberships at the calibration units' features <paramref name="units"/> and the test cases' <paramref name="cases"/>.</summary>
    /// <param name="units">The calibration (or training) units' features, a row each.</param>
    /// <param name="cases">The test cases' features, a row each.</param>
    /// <param name="seed">The draws' seed under the randomized weights.</param>
    /// <param name="grids">Each case's grid, for the density forms; null for the closed forms.</param>
    /// <param name="weighting">The randomized weights, the default, or the KL weights.</param>
    public PosteriorWeights Weights(double[,] units, double[,] cases, int seed = 0, IReadOnlyList<double[]>? grids = null,
                                    PosteriorWeighting weighting = PosteriorWeighting.Randomized)
        => PosteriorConformal.Weights(Memberships(units), Memberships(cases), Precision, seed, grids, weighting);

    private double[][] Rows(double[,] features)
    {
        ArgumentNullException.ThrowIfNull(features);
        if (features.GetLength(1) != _width) throw new ArgumentException($"{_width} features are required.", nameof(features));
        int rows = features.GetLength(0), k = Clusters;
        var result = new double[rows][];
        var row = new double[_width];
        var solver = new SimplexSolver(k);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < _width; j++) row[j] = features[i, j];
            double[] tau = Tau(_fits, _prior, row);
            var linear = new double[k];
            for (int c = 0; c < k; c++)
                for (int t = 0; t < tau.Length; t++) linear[c] += _profiles[c, t] * tau[t];
            var pi = new double[k];
            Array.Fill(pi, 1.0 / k);
            solver.Solve(_gram, linear, pi);
            result[i] = pi;
        }
        return result;
    }

    // τ_t = n_t·g_t/(n − n_t + n_t·g_t) = σ(ln g_t − ln((n − n_t)/n_t)).
    private static double[] Tau(DensityRatio[] fits, double[] prior, ReadOnlySpan<double> features)
    {
        var tau = new double[fits.Length];
        for (int t = 0; t < fits.Length; t++)
        {
            double eta = fits[t].LogRatio(features) - prior[t];
            tau[t] = eta >= 0.0 ? 1.0 / (1.0 + Math.Exp(-eta)) : Math.Exp(eta) / (1.0 + Math.Exp(eta));
        }
        return tau;
    }

    // ---- The factorization τᵢ ≈ Σ_k π_ik·γ_k ----

    // The best of the restarts at J clusters: γ and R²; null where the profiles have fewer than J distinct points.
    internal static (double[,] Gamma, double RSquared)? Factor(double[][] tau, int clusters, int restarts, Xoshiro256PlusPlus rng,
                                                               ScoreMixtureCost? cost = null)
    {
        int n = tau.Length, s = tau[0].Length;
        var mean = new double[s];
        foreach (double[] row in tau)
            for (int t = 0; t < s; t++) mean[t] += row[t] / n;
        double total = 0.0;
        foreach (double[] row in tau)
            for (int t = 0; t < s; t++) total += (row[t] - mean[t]) * (row[t] - mean[t]);
        double[,]? best = null;
        double bestResidual = double.PositiveInfinity;
        for (int r = 0; r < restarts; r++)
        {
            double[,]? gamma = KMeansPlusPlus(tau, clusters, rng);
            if (gamma is null) return null;
            double residual = Alternate(tau, gamma, cost);
            if (residual < bestResidual)
            {
                best = gamma;
                bestResidual = residual;
            }
        }
        return (best!, total > 0.0 ? 1.0 - bestResidual / total : 1.0);
    }

    // k-means++ seeding (Arthur and Vassilvitskii 2007), then Lloyd's iterations; null where fewer than J points are distinct.
    private static double[,]? KMeansPlusPlus(double[][] points, int clusters, Xoshiro256PlusPlus rng)
    {
        int n = points.Length, s = points[0].Length;
        var centres = new double[clusters, s];
        var nearest = new double[n];
        int first = rng.NextInt(n);
        for (int t = 0; t < s; t++) centres[0, t] = points[first][t];
        for (int i = 0; i < n; i++) nearest[i] = Distance(points[i], centres, 0);
        for (int c = 1; c < clusters; c++)
        {
            double sum = 0.0;
            foreach (double d in nearest) sum += d;
            if (!(sum > 0.0)) return null;
            double u = rng.NextDouble() * sum, cumulative = 0.0;
            int pick = n - 1;
            for (int i = 0; i < n; i++)
            {
                cumulative += nearest[i];
                if (u < cumulative && nearest[i] > 0.0) { pick = i; break; }
            }
            while (nearest[pick] == 0.0) pick--;
            for (int t = 0; t < s; t++) centres[c, t] = points[pick][t];
            for (int i = 0; i < n; i++) nearest[i] = Math.Min(nearest[i], Distance(points[i], centres, c));
        }
        var label = new int[n];
        Array.Fill(label, -1);
        for (int iteration = 0; iteration < 100; iteration++)
        {
            bool moved = false;
            for (int i = 0; i < n; i++)
            {
                int to = 0;
                double closest = double.PositiveInfinity;
                for (int c = 0; c < clusters; c++)
                {
                    double d = Distance(points[i], centres, c);
                    if (d < closest) { closest = d; to = c; }
                }
                if (to != label[i]) { label[i] = to; moved = true; }
            }
            if (!moved) break;
            var sums = new double[clusters, s];
            var counts = new int[clusters];
            for (int i = 0; i < n; i++)
            {
                counts[label[i]]++;
                for (int t = 0; t < s; t++) sums[label[i], t] += points[i][t];
            }
            for (int c = 0; c < clusters; c++)
                if (counts[c] > 0)
                    for (int t = 0; t < s; t++) centres[c, t] = sums[c, t] / counts[c];
        }
        return centres;
    }

    private static double Distance(double[] point, double[,] centres, int c)
    {
        double d = 0.0;
        for (int t = 0; t < point.Length; t++) d += (point[t] - centres[c, t]) * (point[t] - centres[c, t]);
        return d;
    }

    // Alternating minimization from γ, updated in place: π given γ by the exact simplex-constrained least squares, γ given π by least
    // squares, until a pass lowers the residual by less than 10⁻⁷ of itself (300 passes at most: where J exceeds the clusters the
    // factorization has flat directions the passes crawl along, and R² and π are settled long before). Returns Σᵢ ‖τᵢ − Σ_k π_ik γ_k‖².
    private static double Alternate(double[][] tau, double[,] gamma, ScoreMixtureCost? cost)
    {
        int n = tau.Length, s = tau[0].Length, k = gamma.GetLength(0);
        var pi = new double[n][];
        for (int i = 0; i < n; i++) { pi[i] = new double[k]; Array.Fill(pi[i], 1.0 / k); }
        var factor = new DenseCholesky(k);
        var solver = new SimplexSolver(k);
        double previous = double.PositiveInfinity, residual = ProjectRows(tau, gamma, pi, solver);
        int passes = 0;
        for (int iteration = 0; iteration < 300 && previous - residual > 1e-7 * residual; iteration++)
        {
            passes++;
            // γ given π: (ΠᵀΠ + εI)γ = ΠᵀT, a ridge of 10⁻¹² of the mean diagonal; a cluster that holds no weight keeps its profile.
            var normal = new double[k, k];
            var right = new double[k, s];
            for (int i = 0; i < n; i++)
                for (int a = 0; a < k; a++)
                {
                    if (pi[i][a] == 0.0) continue;
                    for (int b = 0; b < k; b++) normal[a, b] += pi[i][a] * pi[i][b];
                    for (int t = 0; t < s; t++) right[a, t] += pi[i][a] * tau[i][t];
                }
            double trace = 0.0;
            for (int a = 0; a < k; a++)
            {
                if (normal[a, a] == 0.0)
                {
                    normal[a, a] = 1.0;
                    for (int t = 0; t < s; t++) right[a, t] = gamma[a, t];
                }
                trace += normal[a, a];
            }
            for (int a = 0; a < k; a++) normal[a, a] += 1e-12 * trace / k;
            factor.Decompose(normal);
            if (factor.FailedAt >= 0) break;
            var column = new double[k];
            for (int t = 0; t < s; t++)
            {
                for (int a = 0; a < k; a++) column[a] = right[a, t];
                double[] solved = factor.Solve(column);
                for (int a = 0; a < k; a++) gamma[a, t] = solved[a];
            }
            previous = residual;
            residual = ProjectRows(tau, gamma, pi, solver);
        }
        if (cost is not null)
        {
            cost.Passes += passes;
            cost.Solves += solver.Solves;
            cost.SolverIterations += solver.Total;
            cost.CappedSolves += solver.Capped;
        }
        return residual;
    }

    // Each πᵢ given γ, warm-started from its current value; the residual.
    private static double ProjectRows(double[][] tau, double[,] gamma, double[][] pi, SimplexSolver solver)
    {
        int k = gamma.GetLength(0), s = gamma.GetLength(1);
        double[,] gram = Gram(gamma);
        var linear = new double[k];
        double residual = 0.0;
        for (int i = 0; i < tau.Length; i++)
        {
            for (int c = 0; c < k; c++)
            {
                linear[c] = 0.0;
                for (int t = 0; t < s; t++) linear[c] += gamma[c, t] * tau[i][t];
            }
            solver.Solve(gram, linear, pi[i]);
            for (int t = 0; t < s; t++)
            {
                double fit = 0.0;
                for (int c = 0; c < k; c++) fit += pi[i][c] * gamma[c, t];
                residual += (tau[i][t] - fit) * (tau[i][t] - fit);
            }
        }
        return residual;
    }

    // γγᵀ with a ridge of 10⁻¹² of its mean diagonal, so the simplex problem is strictly convex where two profiles coincide.
    private static double[,] Gram(double[,] gamma)
    {
        int k = gamma.GetLength(0), s = gamma.GetLength(1);
        var gram = new double[k, k];
        double trace = 0.0;
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
            {
                for (int t = 0; t < s; t++) gram[a, b] += gamma[a, t] * gamma[b, t];
                if (a == b) trace += gram[a, a];
            }
        for (int a = 0; a < k; a++) gram[a, a] += 1e-12 * (trace / k + 1e-300);
        return gram;
    }

    /// <summary>The simplex-constrained least squares at <paramref name="pi"/>, overwritten (<see cref="SimplexSolver"/>).</summary>
    internal static void Project(double[,] h, double[] c, double[] pi) => new SimplexSolver(c.Length).Solve(h, c, pi);

    /// <summary>
    /// The simplex-constrained least squares min ½πᵀHπ − cᵀπ over π ≥ 0, Σπ = 1, H = γγᵀ and c = γτ, by a primal active-set method from a
    /// feasible π: on the free coordinates, the equality-constrained step; a step that would leave the simplex stops at its boundary and
    /// fixes the coordinate; where the step's predicted decrease is negligible the point is stationary, and the fixed coordinate with the
    /// most negative multiplier is freed, until none is. Its buffers are reused across calls.
    /// </summary>
    private sealed class SimplexSolver(int k)
    {
        private readonly bool[] _free = new bool[k];
        private readonly int[] _on = new int[k];
        private readonly double[] _g = new double[k], _step = new double[k], _right = new double[k + 1];
        private readonly double[,] _system = new double[k + 1, k + 1];

        /// <summary>The solves, their active-set iterations, and the solves stopped at the cap, over the solver's life.</summary>
        public long Solves, Total, Capped;

        public void Solve(double[,] h, double[] c, double[] pi)
        {
            for (int a = 0; a < k; a++) _free[a] = pi[a] > 0.0;
            double scale = 0.0;
            for (int a = 0; a < k; a++) scale = Math.Max(scale, Math.Max(Math.Abs(h[a, a]), Math.Abs(c[a])));
            double tolerance = 1e-10 * (1.0 + scale);
            int released = -1;
            Solves++;
            for (int iteration = 0; iteration < 50 * (k + 1); iteration++)
            {
                Total++;
                int f = 0;
                for (int a = 0; a < k; a++)
                {
                    _g[a] = -c[a];
                    for (int b = 0; b < k; b++) _g[a] += h[a, b] * pi[b];
                    if (_free[a]) _on[f++] = a;
                }
                // [H_FF −1; 1ᵀ 0][p_F; ν] = [−g_F; 0].
                for (int a = 0; a <= f; a++)
                {
                    for (int b = 0; b <= f; b++)
                        _system[a, b] = a < f && b < f ? h[_on[a], _on[b]] : a == f && b == f ? 0.0 : a == f ? 1.0 : -1.0;
                    _right[a] = a < f ? -_g[_on[a]] : 0.0;
                }
                SolveSmall(_system, _right, f + 1);
                double nu = _right[f], decrease = 0.0;
                Array.Clear(_step);
                for (int a = 0; a < f; a++) _step[_on[a]] = _right[a];
                for (int a = 0; a < k; a++)
                {
                    if (_step[a] == 0.0) continue;
                    double hp = 0.0;
                    for (int b = 0; b < k; b++) hp += h[a, b] * _step[b];
                    decrease -= _step[a] * (_g[a] + 0.5 * hp);
                }
                if (decrease <= 1e-15 * (1.0 + scale))
                {
                    // Stationary on the free set: the fixed coordinates' multipliers g_a − ν, all nonnegative at the optimum.
                    int release = -1;
                    double worst = -tolerance;
                    for (int a = 0; a < k; a++)
                        if (!_free[a] && _g[a] - nu < worst) { worst = _g[a] - nu; release = a; }
                    if (release < 0) return;
                    _free[release] = true;
                    released = release;
                    continue;
                }
                double length = 1.0;
                int block = -1;
                for (int a = 0; a < k; a++)
                    if (_free[a] && _step[a] < 0.0 && -pi[a] / _step[a] < length) { length = -pi[a] / _step[a]; block = a; }
                // The coordinate just freed would be fixed again at once: its multiplier's sign was rounding, and the point is optimal.
                if (block >= 0 && block == released && length == 0.0) return;
                released = -1;
                double sum = 0.0;
                for (int a = 0; a < k; a++)
                {
                    pi[a] = a == block ? 0.0 : Math.Max(0.0, pi[a] + length * _step[a]);
                    sum += pi[a];
                }
                for (int a = 0; a < k; a++) pi[a] /= sum;
                if (block >= 0) _free[block] = false;
            }
            Capped++;
        }
    }

    // Gaussian elimination with partial pivoting on the leading n × n block, the solution left in b.
    private static void SolveSmall(double[,] a, double[] b, int n)
    {
        for (int col = 0; col < n; col++)
        {
            int pivot = col;
            for (int r = col + 1; r < n; r++) if (Math.Abs(a[r, col]) > Math.Abs(a[pivot, col])) pivot = r;
            if (pivot != col)
            {
                for (int j = 0; j < n; j++) (a[col, j], a[pivot, j]) = (a[pivot, j], a[col, j]);
                (b[col], b[pivot]) = (b[pivot], b[col]);
            }
            double d = a[col, col];
            if (d == 0.0) throw new InvalidOperationException("The simplex step's system is singular.");
            for (int r = col + 1; r < n; r++)
            {
                double factor = a[r, col] / d;
                if (factor == 0.0) continue;
                for (int j = col; j < n; j++) a[r, j] -= factor * a[col, j];
                b[r] -= factor * b[col];
            }
        }
        for (int r = n - 1; r >= 0; r--)
        {
            double v = b[r];
            for (int j = r + 1; j < n; j++) v -= a[r, j] * b[j];
            b[r] = v / a[r, r];
        }
    }

    // ---- The precision rule ----

    // The training split's mean effective sample size and mean self-weight at m: each point draws L ~ Multinomial(m, π(xᵢ)) on its own
    // stream and weighs every training point, itself included, by Π_k π_k(xⱼ)^{L_k}.
    internal static (double Ess, double Self) Reading(double[][] memberships, int precision, int seed)
    {
        int n = memberships.Length, k = memberships.Length > 0 ? memberships[0].Length : 0;
        double ess = 0.0, self = 0.0;
        var log = new double[n];
        var logs = new double[n][];
        for (int j = 0; j < n; j++) logs[j] = [.. Array.ConvertAll(memberships[j], v => Math.Log(v))];
        for (int i = 0; i < n; i++)
        {
            int[] draw = PosteriorConformal.Draw(memberships[i], precision, PosteriorConformal.Stream(seed, i));
            double top = double.NegativeInfinity;
            for (int j = 0; j < n; j++)
            {
                // ln Π_k π_jk^{L_k}, a cluster not drawn contributing nothing even where π_jk is 0.
                double l = 0.0;
                for (int c = 0; c < k; c++) if (draw[c] > 0) l += draw[c] * logs[j][c];
                log[j] = l;
                top = Math.Max(top, l);
            }
            double total = 0.0, squares = 0.0;
            for (int j = 0; j < n; j++)
            {
                double w = Math.Exp(log[j] - top);
                total += w;
                squares += w * w;
            }
            ess += total * total / squares;
            self += Math.Exp(log[i] - top) / total;
        }
        return (ess / n, self / n);
    }

    private static ScoreMixtureOptions Checked(ScoreMixtureOptions o)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(o.Levels, 1, nameof(o.Levels));
        if (o.Clusters is int j) ArgumentOutOfRangeException.ThrowIfLessThan(j, 1, nameof(o.Clusters));
        ArgumentOutOfRangeException.ThrowIfLessThan(o.MaxClusters, 1, nameof(o.MaxClusters));
        if (!(o.ClusterGain >= 0.0 && double.IsFinite(o.ClusterGain))) throw new ArgumentOutOfRangeException(nameof(o.ClusterGain), "The gain is finite and nonnegative.");
        if (o.Precision is int m) ArgumentOutOfRangeException.ThrowIfLessThan(m, 1, nameof(o.Precision));
        ArgumentOutOfRangeException.ThrowIfLessThan(o.MinPrecision, 1, nameof(o.MinPrecision));
        if (o.MaxPrecision < o.MinPrecision) throw new ArgumentOutOfRangeException(nameof(o.MaxPrecision), "The precision range is empty.");
        if (!(o.TargetEffectiveSampleSize >= 0.0 && double.IsFinite(o.TargetEffectiveSampleSize))) throw new ArgumentOutOfRangeException(nameof(o.TargetEffectiveSampleSize));
        if (!(o.MaxSelfWeight > 0.0 && o.MaxSelfWeight <= 1.0)) throw new ArgumentOutOfRangeException(nameof(o.MaxSelfWeight), "The self-weight bound is in (0, 1].");
        ArgumentOutOfRangeException.ThrowIfLessThan(o.Restarts, 1, nameof(o.Restarts));
        return o;
    }
}
