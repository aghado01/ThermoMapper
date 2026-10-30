using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;

namespace Maths.Samplers.Tests;

/// <summary>
/// The mean-field (Curie–Weiss) q-state Potts model on N sites, a toy whose first-order transition sharpens with N and
/// whose exact density of states stays cheap far past enumeration. The prior is uniform over the qᴺ configurations and
/// ℓ₁ = (J/N)·Σ_{i&lt;j} δ(σ_i, σ_j) = (J/(2N))·(Σ_k n_k² − N), a function of the occupation counts n alone.
/// <para><b>The state is n.</b> Single-site Metropolis on configurations (a site uniformly, then one of the other q − 1
/// colours uniformly) lumps exactly onto n: from n it moves a site of colour k to l with probability (n_k/N)/(q − 1)
/// whatever the configuration. On n the prior is the multinomial N!/Π n_k!·q^{−N}, the move's reverse has probability
/// (n_l + 1)/(N(q − 1)), and the prior ratio n_k/(n_l + 1) cancels the proposal ratio (n_l + 1)/n_k, so the Hastings
/// ratio on n is the likelihood ratio, as on configurations. A step costs O(q).</para>
/// <para><b>Exact references.</b> ln g sums the multinomial over the n of one Σ n_k² (the compositions of N into q
/// parts, about N^{q−1}/(q − 1)! of them); the exact law at any b and ln Z(b) follow. Levels sit J/N apart, and the
/// count of compositions per level is irregular, so peaks and dips are read on a binned law.</para>
/// <para><b>The transition</b>, in the limit, with x the colour fractions: the free energy per site
/// φ(x) = (βJ/2)·Σ x_k² − Σ x_k ln x_k, on the ordered branch x_1 = (1 + (q − 1)s)/q, x_{k&gt;1} = (1 − s)/q,
/// gains Δφ(s) = (βJ(q − 1)/(2q))·s² − (1/q)·[(1 + (q − 1)s)·ln(1 + (q − 1)s) + (q − 1)(1 − s)·ln(1 − s)] over the
/// disordered point. Δφ′ = ((q − 1)/q)·[βJ·s − ln((1 + (q − 1)s)/(1 − s))] = 0 and Δφ = 0 together give
/// s_c = (q − 2)/(q − 1) and β_cJ = 2(q − 1)·ln(q − 1)/(q − 2) (4 ln 2 at q = 3), with the jump in Σx² of
/// (q − 1)·s_c²/q, so a latent ℓ₁ of J·N·(q − 1)·s_c²/(2q).</para>
/// </summary>
internal sealed class MeanFieldPotts
{
    private readonly double[] _logFactorial;

    public MeanFieldPotts(int sites, int q, double coupling)
    {
        if (sites < 2) throw new ArgumentOutOfRangeException(nameof(sites));
        if (q < 3) throw new ArgumentOutOfRangeException(nameof(q), "The transition is first order from q = 3.");
        Sites = sites;
        Q = q;
        Coupling = coupling;
        _logFactorial = new double[sites + 1];
        for (int k = 2; k <= sites; k++) _logFactorial[k] = _logFactorial[k - 1] + Math.Log(k);
    }

    public int Sites { get; }
    public int Q { get; }
    public double Coupling { get; }

    /// <summary>β_cJ = 2(q − 1)·ln(q − 1)/(q − 2), the limit's transition, derived above.</summary>
    public static double CriticalCoupling(int q) => 2.0 * (q - 1) * Math.Log(q - 1) / (q - 2);

    /// <summary>ℓ₁ at Σ n_k² = <paramref name="squares"/>.</summary>
    public double LevelOf(long squares) => Coupling / (2.0 * Sites) * (squares - Sites);

    public double Level(int[] n)
    {
        long s = 0;
        foreach (int c in n) s += (long)c * c;
        return LevelOf(s);
    }

    /// <summary>ln of the multinomial N!/Π n_k!·q^{−N}: the prior's mass on n.</summary>
    public double LogPrior(int[] n)
    {
        double v = _logFactorial[Sites] - Sites * Math.Log(Q);
        foreach (int c in n)
        {
            if (c < 0) return double.NegativeInfinity;
            v -= _logFactorial[c];
        }
        return v;
    }

    public ISamplingTarget<int[]> Target => SamplingTarget.Separable<int[]>(LogPrior, Level);

    /// <summary>The most even split: the disordered phase's centre.</summary>
    public int[] Disordered()
    {
        var n = new int[Q];
        for (int k = 0; k < Q; k++) n[k] = Sites / Q + (k < Sites % Q ? 1 : 0);
        return n;
    }

    /// <summary>
    /// The exact density of states: each level ℓ₁ with g's log mass there (the prior's mass of the n with that Σ n_k²,
    /// summing to 1), ascending. Every composition of N into q parts is counted, through the partitions
    /// n_1 ≥ … ≥ n_q, each weighed by its q!/Π_v m_v! orderings (m_v the parts equal to v): q! times fewer visits.
    /// </summary>
    public (double[] Levels, double[] LogMasses) DensityOfStates()
    {
        long minSquares = (long)Sites * Sites / Q, maxSquares = (long)Sites * Sites;
        long first = minSquares - ((minSquares - Sites) % 2 + 2) % 2;   // Σ n² ≡ N (mod 2)
        int slots = checked((int)((maxSquares - first) / 2 + 1));
        var max = new double[slots];
        var sum = new double[slots];
        Array.Fill(max, double.NegativeInfinity);
        double constant = _logFactorial[Sites] - Sites * Math.Log(Q) + _logFactorial[Q];
        // Parts in non-increasing order: part k at most the last part and at most what is left; the last part takes the
        // rest when it fits. run is the length of the current run of equal parts, for the orderings' count.
        void Visit(int k, int left, int bound, int previous, int run, double logMass, long squares)
        {
            if (k == Q - 1)
            {
                if (left > bound) return;
                int r = left == previous ? run + 1 : 1;
                double v = logMass - _logFactorial[left] - (left == previous ? Math.Log(r) : 0.0);
                int j = (int)((squares + (long)left * left - first) / 2);
                if (v > max[j])
                {
                    sum[j] = sum[j] * Math.Exp(max[j] - v) + 1.0;
                    max[j] = v;
                }
                else sum[j] += Math.Exp(v - max[j]);
                return;
            }
            int lowest = (left + (Q - k) - 1) / (Q - k);   // the largest of the remaining parts is at least their mean
            for (int c = Math.Min(bound, left); c >= lowest; c--)
            {
                int r = c == previous ? run + 1 : 1;
                // Dividing by r at each step builds Π_v m_v! one factor at a time.
                Visit(k + 1, left - c, c, c, r, logMass - _logFactorial[c] - (c == previous ? Math.Log(r) : 0.0), squares + (long)c * c);
            }
        }
        Visit(0, Sites, Sites, -1, 0, constant, 0);
        var levels = new List<double>();
        var masses = new List<double>();
        for (int j = 0; j < slots; j++)
        {
            if (double.IsNegativeInfinity(max[j])) continue;
            levels.Add(LevelOf(first + 2L * j));
            masses.Add(max[j] + Math.Log(sum[j]));
        }
        return ([.. levels], [.. masses]);
    }
}

/// <summary>The exact density of states of a mean-field Potts model and what it gives: laws, ln Z(b), the transition.</summary>
internal sealed class MeanFieldExact
{
    public MeanFieldExact(MeanFieldPotts model)
    {
        Model = model;
        (Levels, LogMasses) = model.DensityOfStates();
        // Bins for reading peaks and dips: fifty across the limit's latent ℓ₁, so each holds many levels at large N.
        double sc = (model.Q - 2.0) / (model.Q - 1.0);
        double latent = model.Coupling * model.Sites * (model.Q - 1) * sc * sc / (2.0 * model.Q);
        BinWidth = Math.Max(latent / 50.0, model.Coupling / model.Sites);
    }

    public MeanFieldPotts Model { get; }
    public double[] Levels { get; }
    public double[] LogMasses { get; }
    public double BinWidth { get; }

    public double LogEvidence(double b)
    {
        var terms = new double[Levels.Length];
        for (int j = 0; j < terms.Length; j++) terms[j] = LogMasses[j] + b * Levels[j];
        return FlatHistogramToys.LogSumExp(terms);
    }

    /// <summary>The exact law of ℓ₁ at b, per level.</summary>
    public double[] Law(double b)
    {
        double z = LogEvidence(b);
        var p = new double[Levels.Length];
        for (int j = 0; j < p.Length; j++) p[j] = Math.Exp(LogMasses[j] + b * Levels[j] - z);
        return p;
    }

    /// <summary>E_b[f(ℓ₁)] exactly.</summary>
    public double Expectation(double b, Func<double, double> f)
    {
        double[] p = Law(b);
        double s = 0.0;
        for (int j = 0; j < p.Length; j++) s += p[j] * f(Levels[j]);
        return s;
    }

    /// <summary>The law at b binned from the lowest level by <see cref="BinWidth"/>, as log masses: bin i covers
    /// [ℓ_min + i·h, ℓ_min + (i + 1)·h).</summary>
    public double[] BinnedLogLaw(double b)
    {
        double z = LogEvidence(b), lo = Levels[0];
        int bins = (int)((Levels[^1] - lo) / BinWidth) + 1;
        var max = new double[bins];
        var sum = new double[bins];
        Array.Fill(max, double.NegativeInfinity);
        for (int j = 0; j < Levels.Length; j++)
        {
            int i = Math.Min(bins - 1, (int)((Levels[j] - lo) / BinWidth));
            double v = LogMasses[j] + b * Levels[j] - z;
            if (v > max[i])
            {
                sum[i] = sum[i] * Math.Exp(max[i] - v) + 1.0;
                max[i] = v;
            }
            else sum[i] += Math.Exp(v - max[i]);
        }
        var log = new double[bins];
        for (int i = 0; i < bins; i++) log[i] = double.IsNegativeInfinity(max[i]) ? double.NegativeInfinity : max[i] + Math.Log(sum[i]);
        return log;
    }

    public double BinCentre(int bin) => Levels[0] + (bin + 0.5) * BinWidth;

    /// <summary>The binned law smoothed over five bins (in probability), as log masses: the level counts' irregularity
    /// averaged out, so peaks and dips are the phases', not the lattice's.</summary>
    public double[] SmoothedLogLaw(double b)
    {
        double[] log = BinnedLogLaw(b);
        double top = double.NegativeInfinity;
        foreach (double v in log) top = Math.Max(top, v);
        var smooth = new double[log.Length];
        for (int i = 0; i < log.Length; i++)
        {
            double sum = 0.0;
            int count = 0;
            for (int d = -2; d <= 2; d++)
            {
                int k = i + d;
                if (k < 0 || k >= log.Length) continue;
                sum += Math.Exp(log[k] - top);
                count++;
            }
            smooth[i] = top + Math.Log(sum / count);
        }
        return smooth;
    }

    /// <summary>
    /// The smoothed law's two phases at b: the disordered peak (the most likely bin below <paramref name="split"/>), the
    /// ordered peak (above it) and the dip, the least likely bin between them; as bins.
    /// </summary>
    public (int Disordered, int Dip, int Ordered) Phases(double b, double split)
    {
        double[] log = SmoothedLogLaw(b);
        int s = Math.Clamp((int)((split - Levels[0]) / BinWidth), 1, log.Length - 2);
        int dis = 0, ord = s;
        for (int i = 0; i < s; i++) if (log[i] > log[dis]) dis = i;
        for (int i = s; i < log.Length; i++) if (log[i] > log[ord]) ord = i;
        int dip = dis;
        for (int i = dis; i <= ord; i++) if (log[i] < log[dip]) dip = i;
        return (dis, dip, ord);
    }

    /// <summary>The smoothed law's log mass in a bin at b.</summary>
    public double SmoothedLogMass(double b, int bin) => SmoothedLogLaw(b)[bin];

    /// <summary>The level that splits the phases: the least likely bin of the smoothed law at the limit's b_c = β_cJ/J,
    /// between the limit's disordered centre (Σx² = 1/q) and its ordered one at s_c (Σx² = (1 + (q − 1)s_c²)/q).</summary>
    public double Split()
    {
        int q = Model.Q, n = Model.Sites;
        double sc = (q - 2.0) / (q - 1.0);
        double disordered = Model.LevelOf((long)Math.Round((double)n * n / q));
        double ordered = Model.LevelOf((long)Math.Round((double)n * n * (1.0 + (q - 1) * sc * sc) / q));
        double[] log = SmoothedLogLaw(MeanFieldPotts.CriticalCoupling(q) / Model.Coupling);
        int lo = Math.Max(0, (int)((disordered - Levels[0]) / BinWidth) + 1), hi = Math.Min(log.Length - 1, (int)((ordered - Levels[0]) / BinWidth) - 1);
        if (hi <= lo) throw new InvalidOperationException("The phases' centres share a bin at this size.");
        int dip = lo;
        for (int i = lo; i <= hi; i++) if (log[i] < log[dip]) dip = i;
        return BinCentre(dip);
    }

    /// <summary>The ordered phase's mass at b: the law's mass above <paramref name="split"/>.</summary>
    public double OrderedMass(double b, double split)
    {
        double[] p = Law(b);
        double m = 0.0;
        for (int j = 0; j < p.Length; j++) if (Levels[j] > split) m += p[j];
        return m;
    }

    /// <summary>b_t, the equal-weight point: where the phases split at <paramref name="split"/> hold equal mass, by
    /// bisection on (0, 1].</summary>
    public double EqualWeight(double split)
    {
        double lo = 0.0, hi = 1.0;
        for (int i = 0; i < 60; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (OrderedMass(mid, split) < 0.5) lo = mid;
            else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    /// <summary>ℓ₁'s q-quantile under the exact law at b.</summary>
    public double Quantile(double b, double q)
    {
        double[] p = Law(b);
        double c = 0.0;
        for (int j = 0; j < p.Length; j++)
        {
            c += p[j];
            if (c >= q) return Levels[j];
        }
        return Levels[^1];
    }

    public double StandardDeviation(double b)
    {
        double mean = Expectation(b, l => l);
        return Math.Sqrt(Expectation(b, l => (l - mean) * (l - mean)));
    }
}

/// <summary>Single-site Metropolis lumped onto the occupation counts: a site of colour k with probability n_k/N, then
/// one of the other q − 1 colours uniformly; the proposal ratio is (n_l + 1)/n_k.</summary>
internal sealed class OccupationMove(MeanFieldPotts model) : IRjMove<int[]>
{
    public string Key => "site";
    public string ReverseKey => "site";
    public double Weight(int[] state) => 1.0;

    public Proposal<int[]>? Propose(int[] current, Xoshiro256PlusPlus rng)
    {
        int site = rng.NextInt(model.Sites), k = 0;
        for (int seen = current[0]; seen <= site; seen += current[++k]) { }
        int l = (k + 1 + rng.NextInt(model.Q - 1)) % model.Q;
        var next = (int[])current.Clone();
        next[k]--;
        next[l]++;
        return new Proposal<int[]>(next, Math.Log((current[l] + 1.0) / current[k]), 0.0);
    }
}

/// <summary>The mean-field chains: lumped single-site Metropolis on n from the disordered centre; the functional is ℓ₁.</summary>
internal sealed class MeanFieldModel(MeanFieldPotts model) : IMulticanonicalModel<int[]>
{
    public int FunctionalDim => 1;
    public int MonitoredDim => 1;
    public IReadOnlyList<string> FunctionalNames => ["l1"];
    public string Fingerprint => new RunFingerprint().Add(model.Sites).Add(model.Q).Add(model.Coupling).ToString();

    public IMulticanonicalChain<int[]> StartChain(MulticanonicalTarget<int[]> target, int chainSeed, int readoutSeed)
        => new MulticanonicalChain<int[]>(new ReversibleJumpChain<int[]>([new OccupationMove(model)], target, model.Disordered(),
                                                                         new Xoshiro256PlusPlus(chainSeed)), target, Reduce);

    private void Reduce(in int[] n, Span<double> functionals) => functionals[0] = model.Level(n);
}
