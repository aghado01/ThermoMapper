using System;
using System.Collections.Generic;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// A <see cref="SamcLearning"/>'s state, as values: the weights, the iterations taken, the gain's scale t₀ and the
/// iteration its count starts from.
/// </summary>
public sealed record SamcLearningState(MulticanonicalWeights Weights, long Iterations, double GainScale, long GainStart);

/// <summary>
/// Stochastic approximation Monte Carlo (Liang, Liu and Carroll 2007) for the <see cref="MulticanonicalWeights"/>: each
/// iteration's visits move the nodes' log weights toward equal occupancy, with the gain γ_t = t₀/max(t₀, t), the
/// Robbins–Monro sequence their convergence theorem takes (Σγ = ∞, Σγ^ζ &lt; ∞ for ζ in (1, 2)); Belardinelli and
/// Pereyra's 1/t Wang–Landau is the same schedule with t₀ = K. Several chains share the one weight function, each
/// iteration's update averaging their visits (Bornn, Jacob, Del Moral and Doucet 2013, PAWL), summed in chain order, so
/// the learned weights do not depend on how the chains were scheduled. It accumulates the chains' visits into the
/// learned field ln w ≈ −ln g; it runs during burn-in only, and the weights then freeze.
/// <para><b>The update with interpolated nodes.</b> With θ = −ln w, a visit at ℓ inside the window spreads SAMC's
/// indicator over the two nodes of its segment by their hat functions φ_j(ℓ) = 1 − u, φ_{j+1}(ℓ) = u, and the update is
/// <c>θ_k ← θ_k + γ_t·(1/C)·Σ_c 1_in(ℓ_c)·(φ_k(ℓ_c) − π_k)</c>, π the flat occupancy
/// (<see cref="MulticanonicalWeights.FlatOccupancy"/>). A visit outside the window does not update: its term is zero.
/// Since Σ_k φ_k = Σ_k π_k = 1, an update conserves Σ_k θ_k. The mean field is
/// <c>P_in(θ)·(E_θ[φ | in] − π)</c>, which is −P_in times the gradient of the convex
/// <c>G(θ) = ln ∫_in prior·e^{−θ·φ(ℓ)} + π·θ</c>, so its root, the limit, is the unique (up to a constant) θ* with
/// <c>E_θ*[φ_k | in] = π_k</c>: a function of g over the window alone, whatever the extensions. When ln g is piecewise
/// linear on the nodes, θ* = ln g up to a constant. A node whose hat sees no mass falls without bound, and the others'
/// limit takes π_k + d, d their equal share of the empty nodes' π (Liang et al.'s theorem); <see cref="FixedPoint"/>
/// computes θ* for a known density of states.</para>
/// </summary>
public sealed class SamcLearning
{
    private readonly double[] _logWeights;
    private readonly double[] _pi;
    private readonly double[] _visits;   // per update: Σ_c 1_in·φ_k

    /// <param name="start">The weights to learn from: flat, or a warm start.</param>
    /// <param name="gainScale">t₀ ≥ 1: the gain is 1 for the first t₀ iterations, then decays as t₀/t.</param>
    /// <param name="gainStart">The iteration the gain's count starts from, ≥ 0: 0 from flat; a warm start that already
    /// holds t_w iterations' information continues the sequence from t_w.</param>
    public SamcLearning(MulticanonicalWeights start, double gainScale, long gainStart = 0)
        : this(new SamcLearningState(start, 0, gainScale, gainStart))
    {
    }

    /// <summary>Continue a learning from its saved state.</summary>
    public SamcLearning(SamcLearningState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(state.Weights);
        if (!(state.GainScale >= 1.0) || double.IsInfinity(state.GainScale))
            throw new ArgumentOutOfRangeException(nameof(state), "The gain's scale t₀ must be finite and at least 1.");
        if (state.Iterations < 0) throw new ArgumentOutOfRangeException(nameof(state), "The iterations taken must be non-negative.");
        if (state.GainStart < 0) throw new ArgumentOutOfRangeException(nameof(state), "The gain's start must be non-negative.");
        Weights = state.Weights;
        Iterations = state.Iterations;
        GainScale = state.GainScale;
        GainStart = state.GainStart;
        _logWeights = [.. state.Weights.LogWeights];
        _pi = state.Weights.FlatOccupancy();
        _visits = new double[_pi.Length];
    }

    /// <summary>The current weights.</summary>
    public MulticanonicalWeights Weights { get; private set; }

    /// <summary>Iterations taken: updates applied.</summary>
    public long Iterations { get; private set; }

    /// <summary>t₀.</summary>
    public double GainScale { get; }

    /// <summary>The iteration the gain's count starts from.</summary>
    public long GainStart { get; }

    /// <summary>The gain the next update applies: γ at t = <see cref="GainStart"/> + <see cref="Iterations"/> + 1.</summary>
    public double NextGain => Gain(GainScale, GainStart + Iterations + 1);

    /// <summary>γ_t = t₀/max(t₀, t), Liang et al.'s (2007) sequence, for t ≥ 1.</summary>
    public static double Gain(double gainScale, long iteration)
    {
        if (iteration < 1) throw new ArgumentOutOfRangeException(nameof(iteration), "The gain is indexed from t = 1.");
        return gainScale / Math.Max(gainScale, iteration);
    }

    /// <summary>This learning's state, as values.</summary>
    public SamcLearningState Save() => new(Weights, Iterations, GainScale, GainStart);

    /// <summary>
    /// One iteration: every chain's visit, its untempered log likelihood ℓ₁ in chain order, moves ln w by
    /// <c>−γ_t·(1/C)·Σ_c 1_in(ℓ_c)·(φ(ℓ_c) − π)</c>. Returns the new weights.
    /// </summary>
    public MulticanonicalWeights Update(ReadOnlySpan<double> levels)
    {
        if (levels.Length == 0) throw new ArgumentException("An update needs at least one chain's visit.", nameof(levels));
        Array.Clear(_visits);
        int inside = 0;
        foreach (double level in levels)
        {
            if (!Weights.TryLocate(level, out int segment, out double fraction)) continue;
            _visits[segment] += 1.0 - fraction;
            _visits[segment + 1] += fraction;
            inside++;
        }
        double step = NextGain / levels.Length;
        Iterations++;
        if (inside == 0) return Weights;
        for (int k = 0; k < _logWeights.Length; k++) _logWeights[k] -= step * (_visits[k] - inside * _pi[k]);
        Weights = Weights.WithLogWeights(_logWeights);
        return Weights;
    }

    /// <summary>
    /// SAMC's limit θ* for a known density of states, given as a discrete measure: masses <c>e^{logMasses[j]}</c> at
    /// <c>levels[j]</c> (quadrature cells of the prior, an enumeration's energy levels, a WHAM estimate's draws). It
    /// returns the weights over <paramref name="shape"/>'s window and slope whose hat occupancy, conditional on the window,
    /// is flat: <c>E[φ_k | in] = π_k</c> under <c>mass·w</c>, the minimizer of the convex G by Newton's method. Nodes whose
    /// hat sees no mass are empty: the others take π_k + d (<see cref="SamcLearning"/>), and an empty node's ln w is read
    /// off the line through its nearest non-empty neighbours, extrapolated beyond the last. The constant is fixed by the
    /// first non-empty node's ln w = 0. Levels outside the window, and masses of −∞, play no part. A solve that cannot
    /// bring every node's occupancy within 1e−3 of flat is refused (<see cref="InvalidOperationException"/>).
    /// </summary>
    /// <param name="shape">The window, slope and node count; its values are not read.</param>
    /// <param name="levels">The measure's levels ℓ₁.</param>
    /// <param name="logMasses">Each level's log mass.</param>
    /// <param name="empty">Per node, whether its hat saw no mass.</param>
    public static MulticanonicalWeights FixedPoint(MulticanonicalWeights shape, ReadOnlySpan<double> levels,
                                                  ReadOnlySpan<double> logMasses, out bool[] empty)
    {
        ArgumentNullException.ThrowIfNull(shape);
        if (levels.Length != logMasses.Length) throw new ArgumentException("One log mass per level is required.", nameof(logMasses));
        int k = shape.Nodes;

        // The measure inside the window, as (segment, fraction, log mass).
        var segments = new List<int>();
        var fractions = new List<double>();
        var masses = new List<double>();
        var seen = new bool[k];
        for (int j = 0; j < levels.Length; j++)
        {
            if (double.IsNaN(logMasses[j]) || double.IsPositiveInfinity(logMasses[j]))
                throw new ArgumentException($"Log mass {j} is not a finite number or −∞.", nameof(logMasses));
            if (double.IsNegativeInfinity(logMasses[j]) || !shape.TryLocate(levels[j], out int s, out double u)) continue;
            segments.Add(s);
            fractions.Add(u);
            masses.Add(logMasses[j]);
            if (u < 1.0) seen[s] = true;
            if (u > 0.0) seen[s + 1] = true;
        }
        if (masses.Count == 0) throw new ArgumentException("No mass falls inside the window.", nameof(levels));

        // The free nodes are the non-empty ones bar the first, which is pinned at θ = 0.
        empty = new bool[k];
        int full = 0;
        for (int i = 0; i < k; i++)
        {
            empty[i] = !seen[i];
            if (seen[i]) full++;
        }
        double[] pi = shape.FlatOccupancy();
        double emptyShare = 0.0;
        for (int i = 0; i < k; i++) if (empty[i]) emptyShare += pi[i];
        var target = new double[k];
        for (int i = 0; i < k; i++) target[i] = empty[i] ? 0.0 : pi[i] + emptyShare / full;

        int pin = Array.IndexOf(seen, true);
        var free = new List<int>(full);
        for (int i = 0; i < k; i++) if (seen[i] && i != pin) free.Add(i);
        var index = new int[k];
        Array.Fill(index, -1);
        for (int f = 0; f < free.Count; f++) index[free[f]] = f;

        // Newton starts from ln g at each node, fitted to the measure near it: the node's cell (half a spacing either
        // side, inside the window) is cut into sub-bins, their log masses summed in logs about each one's largest term (a
        // density of states can span more than a double's 745 nats), and a line through the occupied ones is read at the
        // node, less ln(sinh(βh/2)/(βh/2)) for its slope β over the sub-bins' width h, the bias of an exponential density
        // averaged over a bin. So the start keeps ln g's slope however steep it is against the spacing, where hat or cell
        // masses lean on their low or high edges. A cell holding no level takes the line through its nearest held
        // neighbours. From a start far down a steep ln g, Newton would stall in damped steps.
        const int subBins = 8;
        double sub = shape.Spacing / subBins;
        var subMax = new double[k * subBins];
        var subSum = new double[k * subBins];
        Array.Fill(subMax, double.NegativeInfinity);
        for (int j = 0; j < masses.Count; j++)
        {
            int node = fractions[j] < 0.5 ? segments[j] : segments[j] + 1;
            double offsetInCell = (segments[j] + fractions[j] - node) * shape.Spacing;   // in [−Δ/2, Δ/2]
            int b = Math.Clamp((int)Math.Floor((offsetInCell + 0.5 * shape.Spacing) / sub), 0, subBins - 1);
            int slot = node * subBins + b;
            double logTerm = masses[j];
            if (logTerm > subMax[slot])
            {
                subSum[slot] = subSum[slot] * Math.Exp(subMax[slot] - logTerm) + 1.0;
                subMax[slot] = logTerm;
            }
            else subSum[slot] += Math.Exp(logTerm - subMax[slot]);
        }
        var theta = new double[k];
        var cellEmpty = new bool[k];
        for (int i = 0; i < k; i++)
        {
            var xs = new List<double>(subBins);
            var ys = new List<double>(subBins);
            for (int b = 0; b < subBins; b++)
            {
                int slot = i * subBins + b;
                if (double.IsNegativeInfinity(subMax[slot])) continue;
                xs.Add((b + 0.5) * sub - 0.5 * shape.Spacing);
                ys.Add(subMax[slot] + Math.Log(subSum[slot]) - Math.Log(sub));
            }
            cellEmpty[i] = xs.Count == 0;
            if (cellEmpty[i]) continue;
            if (xs.Count == 1)
            {
                theta[i] = ys[0];
                continue;
            }
            double mx = 0.0, my = 0.0;
            for (int t = 0; t < xs.Count; t++) { mx += xs[t]; my += ys[t]; }
            mx /= xs.Count;
            my /= xs.Count;
            double sxy = 0.0, sxx = 0.0;
            for (int t = 0; t < xs.Count; t++) { sxy += (xs[t] - mx) * (ys[t] - my); sxx += (xs[t] - mx) * (xs[t] - mx); }
            double slope = sxx > 0.0 ? sxy / sxx : 0.0;
            double half = 0.5 * Math.Abs(slope) * sub;
            double bias = half < 1e-8 ? 0.0 : half + Math.Log((1.0 - Math.Exp(-2.0 * half)) / (2.0 * half));   // ln(sinh(h)/h)
            theta[i] = my - slope * mx - bias;
        }
        FillEmpty(theta, cellEmpty);
        double offset = theta[pin];
        for (int i = 0; i < k; i++) theta[i] = seen[i] ? theta[i] - offset : 0.0;
        int n = free.Count;
        var gradient = new double[n];
        var hessian = new double[n * n];
        var mean = new double[k];
        var second = new double[3 * k];   // E[φ_iφ_i], E[φ_iφ_{i+1}] per node
        var exponents = new double[masses.Count];
        double value = Objective(theta, segments, fractions, masses, target, mean, second, exponents);
        for (int iteration = 0; iteration < 200 && n > 0; iteration++)
        {
            // G's gradient is π' − E[φ | in] and its Hessian Cov[φ | in], tridiagonal less a rank one.
            double norm = 0.0;
            for (int f = 0; f < n; f++)
            {
                gradient[f] = target[free[f]] - mean[free[f]];
                norm = Math.Max(norm, Math.Abs(gradient[f]));
            }
            if (norm < 1e-10) break;
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++) hessian[a * n + b] = -mean[free[a]] * mean[free[b]];
            for (int a = 0; a < n; a++)
            {
                int i = free[a];
                hessian[a * n + a] += second[3 * i];
                if (i + 1 < k && index[i + 1] is int up and >= 0)
                {
                    hessian[a * n + up] += second[3 * i + 1];
                    hessian[up * n + a] += second[3 * i + 1];
                }
            }
            double[] direction = SolveSpd(hessian, gradient, n);

            // Backtracking on G, which is convex: the Newton step is a descent direction. Once a step no longer lowers G
            // beyond rounding, the occupancy is as flat as doubles resolve, and the solve stops there.
            double slope = 0.0;
            for (int f = 0; f < n; f++) slope -= gradient[f] * direction[f];
            double t = 1.0;
            var trial = (double[])theta.Clone();
            double trialValue = double.NaN;
            bool descended = false;
            for (int halving = 0; halving < 40; halving++, t *= 0.5)
            {
                for (int f = 0; f < n; f++) trial[free[f]] = theta[free[f]] - t * direction[f];
                trialValue = Objective(trial, segments, fractions, masses, target, mean, second, exponents);
                if (trialValue <= value + 1e-4 * t * slope) { descended = true; break; }
            }
            if (!descended)
            {
                Objective(theta, segments, fractions, masses, target, mean, second, exponents);
                break;
            }
            Array.Copy(trial, theta, k);
            // A step that no longer moves G ends the solve only once the occupancy is flat to 1e−8: rounding then
            // limits it. Short of that, a small step is a damped step on the way, not convergence.
            bool stalled = value - trialValue <= 1e-14 * Math.Max(1.0, Math.Abs(value)) && norm < 1e-8;
            value = trialValue;
            if (stalled) break;
        }

        // A solve that stopped short of flat is refused rather than returned: its weights would not be SAMC's limit.
        double worst = 0.0;
        for (int i = 0; i < k; i++) if (seen[i]) worst = Math.Max(worst, Math.Abs(mean[i] / target[i] - 1.0));
        if (!(worst <= 1e-3))
            throw new InvalidOperationException($"SAMC's limit was not reached: a node's occupancy departs from flat by {worst:G3}.");

        var logWeights = new double[k];
        for (int i = 0; i < k; i++) logWeights[i] = -theta[i];
        FillEmpty(logWeights, empty);
        return shape.WithLogWeights(logWeights);
    }

    // G(θ) = ln Σ_j m_j e^{−θ·φ(ℓ_j)} + π'·θ, with E[φ | in] and the hats' second moments under the tilted measure.
    private static double Objective(double[] theta, List<int> segments, List<double> fractions, List<double> masses,
                                    double[] target, double[] mean, double[] second, double[] exponents)
    {
        int count = masses.Count;
        double max = double.NegativeInfinity;
        for (int j = 0; j < count; j++)
        {
            int s = segments[j];
            double u = fractions[j];
            exponents[j] = masses[j] - ((1.0 - u) * theta[s] + u * theta[s + 1]);
            if (exponents[j] > max) max = exponents[j];
        }
        Array.Clear(mean);
        Array.Clear(second);
        double total = 0.0;
        for (int j = 0; j < count; j++)
        {
            double p = Math.Exp(exponents[j] - max);
            total += p;
            int s = segments[j];
            double u = fractions[j];
            mean[s] += p * (1.0 - u);
            mean[s + 1] += p * u;
            second[3 * s] += p * (1.0 - u) * (1.0 - u);
            second[3 * (s + 1)] += p * u * u;
            second[3 * s + 1] += p * (1.0 - u) * u;
        }
        for (int i = 0; i < mean.Length; i++)
        {
            mean[i] /= total;
            second[3 * i] /= total;
            second[3 * i + 1] /= total;
        }
        double linear = 0.0;
        for (int i = 0; i < theta.Length; i++) linear += target[i] * theta[i];
        return max + Math.Log(total) + linear;
    }

    // An empty node's ln w on the line through its nearest non-empty neighbours, or beyond the last two.
    private static void FillEmpty(double[] logWeights, bool[] empty)
    {
        int k = logWeights.Length;
        var full = new List<int>();
        for (int i = 0; i < k; i++) if (!empty[i]) full.Add(i);
        if (full.Count == k) return;
        for (int i = 0; i < k; i++)
        {
            if (!empty[i]) continue;
            int lo = -1, hi = -1;
            foreach (int f in full)
            {
                if (f < i) lo = f;
                else if (hi < 0) hi = f;
            }
            if (lo < 0 && full.Count >= 2) (lo, hi) = (full[0], full[1]);
            else if (hi < 0 && full.Count >= 2) (lo, hi) = (full[^2], full[^1]);
            if (lo < 0 || hi < 0)
            {
                logWeights[i] = logWeights[lo < 0 ? hi : lo];
                continue;
            }
            logWeights[i] = logWeights[lo] + (logWeights[hi] - logWeights[lo]) * (i - lo) / (hi - lo);
        }
    }

    // Solve H x = g for a symmetric positive definite H (n × n, row-major) by Cholesky, with a ridge where it is singular
    // to rounding. WHAM's Newton steps share it.
    internal static double[] SolveSpd(double[] h, double[] g, int n)
    {
        var l = new double[n * n];
        double scale = 0.0;
        for (int i = 0; i < n; i++) scale = Math.Max(scale, h[i * n + i]);
        double ridge = 1e-14 * Math.Max(scale, 1e-300);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double s = h[i * n + j];
                for (int m = 0; m < j; m++) s -= l[i * n + m] * l[j * n + m];
                if (i == j) l[i * n + i] = Math.Sqrt(Math.Max(s, 0.0) + ridge);
                else l[i * n + j] = s / l[j * n + j];
            }
        }
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = g[i];
            for (int m = 0; m < i; m++) s -= l[i * n + m] * y[m];
            y[i] = s / l[i * n + i];
        }
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double s = y[i];
            for (int m = i + 1; m < n; m++) s -= l[m * n + i] * x[m];
            x[i] = s / l[i * n + i];
        }
        return x;
    }
}
