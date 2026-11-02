using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Regression.Tests;

/// <summary>
/// A finite reversible chain as the spectral-gap brief reads it: π in logs and the off-diagonal K in logs, by state.
/// Built from <see cref="BassKernelEnumeration"/>'s flows, or by hand for the solver's closed-form checks.
/// </summary>
internal sealed class ReversibleChain
{
    private ReversibleChain(double[] logPi, int[] starts, int[] to, double[] logRate)
    {
        LogPi = logPi;
        Starts = starts;
        To = to;
        LogRate = logRate;
        Rate = [.. logRate.Select(Math.Exp)];
        Hold = new double[logPi.Length];
        for (int i = 0; i < Size; i++)
        {
            double leaving = 0.0;
            for (int e = starts[i]; e < starts[i + 1]; e++) leaving += Rate[e];
            Hold[i] = Math.Max(0.0, 1.0 - leaving);
        }
    }

    public int Size => LogPi.Length;

    /// <summary>ln π, normalized.</summary>
    public double[] LogPi { get; }

    /// <summary>Each state's first edge in <see cref="To"/>, and the end.</summary>
    public int[] Starts { get; }

    public int[] To { get; }

    /// <summary>ln K(S, S′) by edge.</summary>
    public double[] LogRate { get; }

    public double[] Rate { get; }

    /// <summary>K(S, S).</summary>
    public double[] Hold { get; }

    public int Transitions => To.Length;

    /// <summary>The chain with π ∝ exp(<paramref name="logPi"/>) and the transitions of <paramref name="flows"/> summed.</summary>
    public static ReversibleChain Of(double[] logPi, IEnumerable<IReadOnlyDictionary<(int From, int To), double>> flows)
    {
        double total = LogMath.Sum(logPi);
        var rows = new SortedDictionary<int, double>[logPi.Length];
        for (int i = 0; i < rows.Length; i++) rows[i] = [];
        foreach (IReadOnlyDictionary<(int From, int To), double> group in flows)
            foreach (((int from, int to), double logK) in group)
                rows[from][to] = rows[from].TryGetValue(to, out double held) ? LogMath.Add(held, logK) : logK;
        var starts = new int[logPi.Length + 1];
        for (int i = 0; i < rows.Length; i++) starts[i + 1] = starts[i] + rows[i].Count;
        var targets = new int[starts[^1]];
        var rate = new double[starts[^1]];
        for (int i = 0; i < rows.Length; i++)
        {
            int e = starts[i];
            foreach ((int j, double logK) in rows[i])
            {
                targets[e] = j;
                rate[e++] = logK;
            }
        }
        return new ReversibleChain([.. logPi.Select(v => v - total)], starts, targets, rate);
    }

    /// <summary>The communicating classes.</summary>
    public int[][] Classes()
    {
        var parent = new int[Size];
        for (int i = 0; i < Size; i++) parent[i] = i;
        int Find(int i)
        {
            while (parent[i] != i) i = parent[i] = parent[parent[i]];
            return i;
        }
        for (int i = 0; i < Size; i++)
            for (int e = Starts[i]; e < Starts[i + 1]; e++) parent[Find(To[e])] = Find(i);
        return [.. Enumerable.Range(0, Size).GroupBy(Find).Select(g => g.ToArray())];
    }

    /// <summary>The chain on a communicating class, π renormalized.</summary>
    public ReversibleChain On(int[] states)
    {
        var index = new Dictionary<int, int>(states.Length);
        for (int k = 0; k < states.Length; k++) index.Add(states[k], k);
        var flows = new Dictionary<(int From, int To), double>();
        foreach (int i in states)
            for (int e = Starts[i]; e < Starts[i + 1]; e++) flows[(index[i], index[To[e]])] = LogRate[e];
        return Of([.. states.Select(i => LogPi[i])], [flows]);
    }

    /// <summary>D^(1/2) K D^(−1/2), each off-diagonal √(K(S, S′)K(S′, S)).</summary>
    public double[,] Symmetrized()
    {
        var s = new double[Size, Size];
        var log = new double[Size, Size];
        for (int i = 0; i < Size; i++)
        {
            s[i, i] = Hold[i];
            for (int e = Starts[i]; e < Starts[i + 1]; e++) log[i, To[e]] = LogRate[e];
        }
        for (int i = 0; i < Size; i++)
            for (int e = Starts[i]; e < Starts[i + 1]; e++)
            {
                int j = To[e];
                if (j > i) s[i, j] = s[j, i] = Math.Exp(0.5 * (log[i, j] + log[j, i]));
            }
        return s;
    }

    /// <summary>μK, into <paramref name="next"/>.</summary>
    public void Step(double[] mu, double[] next)
    {
        for (int i = 0; i < Size; i++) next[i] = mu[i] * Hold[i];
        for (int i = 0; i < Size; i++)
        {
            double m = mu[i];
            if (m == 0.0) continue;
            for (int e = Starts[i]; e < Starts[i + 1]; e++) next[To[e]] += m * Rate[e];
        }
    }

    /// <summary>Kf, into <paramref name="next"/>.</summary>
    public void Apply(double[] f, double[] next)
    {
        for (int i = 0; i < Size; i++)
        {
            double v = Hold[i] * f[i];
            for (int e = Starts[i]; e < Starts[i + 1]; e++) v += Rate[e] * f[To[e]];
            next[i] = v;
        }
    }

    /// <summary>π normalized.</summary>
    public double[] Pi() => [.. LogPi.Select(Math.Exp)];

    /// <summary>Each pair's conductance π(S)K(S, S′), by S &lt; S′, descending.</summary>
    public List<(double C, int I, int J)> EdgesByConductance()
    {
        var edges = new List<(double C, int I, int J)>(Transitions / 2);
        for (int i = 0; i < Size; i++)
            for (int e = Starts[i]; e < Starts[i + 1]; e++)
                if (To[e] > i) edges.Add((LogPi[i] + LogRate[e], i, To[e]));
        edges.Sort((a, b) => b.C.CompareTo(a.C));
        return edges;
    }
}

/// <summary>Sums and adds in logs.</summary>
internal static class LogMath
{
    public static double Add(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double top = Math.Max(a, b);
        return top + Math.Log(Math.Exp(a - top) + Math.Exp(b - top));
    }

    public static double Sum(ReadOnlySpan<double> values)
    {
        double top = double.NegativeInfinity;
        foreach (double v in values) top = Math.Max(top, v);
        if (double.IsNegativeInfinity(top)) return top;
        double sum = 0.0;
        foreach (double v in values) sum += Math.Exp(v - top);
        return top + Math.Log(sum);
    }
}

/// <summary>
/// I − K eliminated with no subtraction (Grassmann, Taksar and Heyman): states go ascending π, each pivot the sum of its
/// row's remaining rates, so the factors keep K's off-diagonals' relative accuracy however small γ is. Solves the
/// Poisson equation in the symmetrized scale, and reads the slowest modes by Lanczos on the inverse orthogonal to √π.
/// </summary>
internal sealed class EliminatedChain
{
    // The smallest pivot read, about e⁻⁶⁶⁷: past it the solves can overflow, so the rung says t_rel is beyond range.
    private const double PivotFloor = 1e-290;

    private readonly int[] _order;
    private readonly double[] _sqrtPi;     // √π, ascending π
    private readonly double[][] _lower;    // by pivot k, the multipliers for rows k + 1 … n − 2
    private readonly double[][] _upper;    // by pivot k, the rates to columns k + 1 … n − 2
    private readonly double[] _pivot;

    private EliminatedChain(int[] order, double[] sqrtPi, double[][] lower, double[][] upper, double[] pivot)
    {
        _order = order;
        _sqrtPi = sqrtPi;
        _lower = lower;
        _upper = upper;
        _pivot = pivot;
    }

    public int Size => _order.Length;

    /// <summary>The elimination of <paramref name="chain"/>, which must be irreducible; null when a pivot falls below
    /// the range doubles can solve with (t_rel beyond about e⁶⁶⁷).</summary>
    public static EliminatedChain? Of(ReversibleChain chain, int threads = 1)
    {
        int n = chain.Size;
        int[] order = [.. Enumerable.Range(0, n).OrderBy(i => chain.LogPi[i]).ThenBy(i => i)];
        var position = new int[n];
        for (int k = 0; k < n; k++) position[order[k]] = k;
        var rows = new double[n][];
        for (int k = 0; k < n; k++)
        {
            var row = new double[n];
            int i = order[k];
            for (int e = chain.Starts[i]; e < chain.Starts[i + 1]; e++) row[position[chain.To[e]]] = chain.Rate[e];
            rows[k] = row;
        }
        double[] log = [.. order.Select(i => chain.LogPi[i])];
        var pivot = new double[n - 1];
        var options = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, threads) };
        for (int k = 0; k < n - 1; k++)
        {
            double[] rk = rows[k];
            double d = 0.0;
            for (int j = k + 1; j < n; j++) d += rk[j];
            if (!(d > PivotFloor)) return null;
            pivot[k] = d;
            int pivotRow = k;
            void Eliminate(int i)
            {
                double[] ri = rows[i];
                double f = ri[pivotRow] / d;
                if (f == 0.0) return;
                int j = pivotRow + 1, width = Vector<double>.Count;
                var scale = new Vector<double>(f);
                for (; j <= n - width; j += width)
                    (new Vector<double>(ri, j) + scale * new Vector<double>(rk, j)).CopyTo(ri, j);
                for (; j < n; j++) ri[j] += f * rk[j];
            }
            // The root's row is never used, and K's row i at column i is a self-loop no pivot reads.
            if (threads > 1 && n - k > 128) Parallel.For(k + 1, n - 1, options, Eliminate);
            else for (int i = k + 1; i < n - 1; i++) Eliminate(i);
        }
        var lower = new double[n - 1][];
        var upper = new double[n - 1][];
        for (int k = 0; k < n - 1; k++)
        {
            int width = Math.Max(0, n - 2 - k);
            var l = new double[width];
            var u = new double[width];
            double logD = Math.Log(pivot[k]);
            for (int m = 0; m < width; m++)
            {
                int i = k + 1 + m;
                double down = rows[i][k], up = rows[k][i];
                if (down > 0.0) l[m] = Math.Exp(Math.Log(down) - logD + 0.5 * (log[i] - log[k]));
                if (up > 0.0) u[m] = Math.Exp(Math.Log(up) + 0.5 * (log[k] - log[i]));
            }
            lower[k] = l;
            upper[k] = u;
        }
        double[] sqrtPi = [.. log.Select(v => Math.Exp(0.5 * v))];
        double norm = Math.Sqrt(sqrtPi.Sum(v => v * v));
        for (int k = 0; k < n; k++) sqrtPi[k] /= norm;
        return new EliminatedChain(order, sqrtPi, lower, upper, pivot);
    }

    // (I − S)x = y orthogonal to √π, in the ascending order; the root's x is 0 until the projection.
    private double[] Solve(double[] y)
    {
        int n = Size;
        double[] x = [.. y];
        for (int k = 0; k < n - 1; k++)
        {
            double xk = x[k];
            if (xk == 0.0) continue;
            double[] l = _lower[k];
            for (int m = 0; m < l.Length; m++) x[k + 1 + m] += l[m] * xk;
        }
        x[n - 1] = 0.0;
        for (int k = n - 2; k >= 0; k--)
        {
            double s = x[k];
            double[] u = _upper[k];
            for (int m = 0; m < u.Length; m++) s += u[m] * x[k + 1 + m];
            x[k] = s / _pivot[k];
        }
        Remove(x, _sqrtPi);
        return x;
    }

    private static void Remove(double[] v, double[] unit)
    {
        double dot = Dot(v, unit);
        for (int i = 0; i < v.Length; i++) v[i] -= dot * unit[i];
    }

    private static double Dot(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
        return s;
    }

    // ‖v‖, scaled so a huge v's square does not overflow.
    private static double Norm(double[] v)
    {
        double top = 0.0;
        foreach (double x in v) top = Math.Max(top, Math.Abs(x));
        if (top == 0.0) return 0.0;
        double s = 0.0;
        foreach (double x in v) s += x / top * (x / top);
        return top * Math.Sqrt(s);
    }

    private static void Scale(double[] v, double by)
    {
        for (int i = 0; i < v.Length; i++) v[i] *= by;
    }

    /// <summary>The next slowest mode orthogonal to √π and <paramref name="deflated"/>, by Lanczos with full
    /// reorthogonalization: its gap, its unit vector in the symmetrized scale (by state, in the chain's order) and the
    /// relative residual.</summary>
    public (double Gap, double[] Vector, double Residual) Slowest(IReadOnlyList<double[]> deflated, int seed, int maxIterations = 90)
    {
        int n = Size;
        double[][] against = [_sqrtPi, .. deflated.Select(ToAscending)];
        void Clean(double[] v)
        {
            foreach (double[] u in against) Remove(v, u);
        }
        var rng = new Xoshiro256PlusPlus(seed);
        var start = new double[n];
        for (int i = 0; i < n; i++) start[i] = BassQualification.Gaussian(rng);
        Clean(start);
        Scale(start, 1.0 / Norm(start));
        double[] first = Solve(start);
        Clean(first);
        double unit = Norm(first);   // the operator's scale, so the tridiagonal stays near 1

        var basis = new List<double[]> { start };
        var alpha = new List<double>();
        var beta = new List<double>();
        double[] ritz = [];
        double top = 0.0;
        for (int m = 0; m < Math.Min(maxIterations, n - against.Length); m++)
        {
            double[] w = m == 0 ? first : Solve(basis[m]);
            if (m > 0) Clean(w);
            Scale(w, 1.0 / unit);
            alpha.Add(Dot(w, basis[m]));
            for (int pass = 0; pass < 2; pass++)
                foreach (double[] b in basis)
                {
                    double dot = Dot(w, b);
                    for (int i = 0; i < n; i++) w[i] -= dot * b[i];
                }
            double next = Norm(w);
            var t = new double[m + 1, m + 1];
            for (int i = 0; i <= m; i++)
            {
                t[i, i] = alpha[i];
                if (i < m) t[i, i + 1] = t[i + 1, i] = beta[i];
            }
            EigenResult tri = DenseEigen.DecomposeSymmetric(t, 256, 1e-15);
            top = tri.Eigenvalues[0];
            ritz = tri.Eigenvectors[0];
            bool done = Math.Abs(next * ritz[m]) <= 1e-10 * Math.Abs(top) || next <= 1e-13;
            if (done) break;
            beta.Add(next);
            Scale(w, 1.0 / next);
            basis.Add(w);
        }
        var y = new double[n];
        for (int k = 0; k < ritz.Length; k++)
            for (int i = 0; i < n; i++) y[i] += ritz[k] * basis[k][i];
        Clean(y);
        Scale(y, 1.0 / Norm(y));

        // One step of inverse iteration, then the Rayleigh quotient and residual on the operator itself.
        double[] refined = Solve(y);
        Clean(refined);
        Scale(refined, 1.0 / unit);
        double refinedNorm = Norm(refined);
        Scale(refined, 1.0 / refinedNorm);
        double[] again = Solve(refined);
        Clean(again);
        Scale(again, 1.0 / unit);
        double theta = Dot(again, refined);
        double residual = 0.0;
        for (int i = 0; i < n; i++)
        {
            double r = again[i] - theta * refined[i];
            residual += r * r;
        }
        double[] vector = new double[n];
        for (int k = 0; k < n; k++) vector[_order[k]] = refined[k];
        return (1.0 / (theta * unit), vector, Math.Sqrt(residual) / Math.Abs(theta));
    }

    private double[] ToAscending(double[] byState)
    {
        var v = new double[Size];
        for (int k = 0; k < Size; k++) v[k] = byState[_order[k]];
        return v;
    }

    /// <summary>f's integrated autocorrelation time under K, from the Poisson equation: 2⟨f̃, (I − K)⁻¹f̃⟩_π / Var f − 1;
    /// NaN where Var f is under 10⁻⁸, since f − E f is then mostly the mean's rounding.</summary>
    public double IntegratedTime(double[] f)
    {
        double mean = 0.0;
        for (int k = 0; k < Size; k++) mean += _sqrtPi[k] * _sqrtPi[k] * f[_order[k]];
        var y = new double[Size];
        double variance = 0.0;
        for (int k = 0; k < Size; k++)
        {
            y[k] = _sqrtPi[k] * (f[_order[k]] - mean);
            variance += y[k] * y[k];
        }
        if (!(variance >= 1e-8)) return double.NaN;
        double[] x = Solve(y);
        return 2.0 * Dot(x, y) / variance - 1.0;
    }
}

/// <summary>SC-2's reading of one set: its size, its mass and its complement's in logs, the flow out, the bound
/// L(A) = π(A)(1 − π(A)) / flow(A), and the barrier D(A) = ln π(A) − ln of the largest conductance leaving.</summary>
internal sealed record SetReading(string Family, string Label, int[] Members, double LogMass, double LogComplement, double LogFlow, double Barrier)
{
    public int Size => Members.Length;
    public double LogBound => LogMass + LogComplement - LogFlow;

    /// <summary>The barrier read from the complement's side: the same largest conductance against its mass.</summary>
    public double ComplementBarrier => Barrier - LogMass + LogComplement;

    /// <summary>Whether A and its complement both hold at least 1%, the brief's relevant sets.</summary>
    public bool Relevant => LogMass >= Math.Log(0.01) && LogComplement >= Math.Log(0.01);
}

/// <summary>SC-2 over the brief's four families: the sets tried, the bottleneck and the relevant bottleneck, each read
/// from its lighter side.</summary>
internal sealed class SetSurvey
{
    private readonly ReversibleChain _chain;
    private readonly double[] _buffer;

    public SetSurvey(ReversibleChain chain)
    {
        _chain = chain;
        _buffer = new double[chain.Size + chain.Transitions];
    }

    public int Tried { get; private set; }
    public SetReading? Bottleneck { get; private set; }
    public SetReading? RelevantBottleneck { get; private set; }

    private void Offer(SetReading reading)
    {
        Tried++;
        if (Bottleneck is null || reading.LogBound > Bottleneck.LogBound) Bottleneck = Lighter(reading);
        if (reading.Relevant && (RelevantBottleneck is null || reading.LogBound > RelevantBottleneck.LogBound)) RelevantBottleneck = Lighter(reading);
    }

    // L(A) = L(Aᶜ), so a cut is read from its lighter side: what is held there is what leaves slowly.
    private SetReading Lighter(SetReading set)
    {
        if (set.LogMass <= set.LogComplement) return set;
        var inside = new bool[_chain.Size];
        foreach (int i in set.Members) inside[i] = true;
        return set with
        {
            Members = [.. Enumerable.Range(0, _chain.Size).Where(i => !inside[i])],
            LogMass = set.LogComplement, LogComplement = set.LogMass, Barrier = set.ComplementBarrier,
        };
    }

    /// <summary>Reads the set <paramref name="members"/> (of 1 … n − 1 states).</summary>
    public void Read(string family, string label, int[] members)
    {
        if (members.Length == 0 || members.Length == _chain.Size) return;
        var inside = new bool[_chain.Size];
        foreach (int i in members) inside[i] = true;
        Offer(Measure(family, label, members, i => inside[i]));
    }

    // The set `members` is `isIn`; flow(A) = flow(Aᶜ) by reversibility, so the edges read are the smaller side's.
    private SetReading Measure(string family, string label, int[] members, Func<int, bool> isIn)
    {
        int count = 0;
        Span<double> values = _buffer;
        foreach (int i in members) values[count++] = _chain.LogPi[i];
        double logMass = LogMath.Sum(values[..count]);
        count = 0;
        for (int i = 0; i < _chain.Size; i++)
            if (!isIn(i)) values[count++] = _chain.LogPi[i];
        double logComplement = LogMath.Sum(values[..count]);
        double leaving = double.NegativeInfinity;
        count = 0;
        void Edges(int i, bool wantInside)
        {
            for (int e = _chain.Starts[i]; e < _chain.Starts[i + 1]; e++)
                if (isIn(_chain.To[e]) == wantInside)
                {
                    double c = _chain.LogPi[i] + _chain.LogRate[e];
                    _buffer[count++] = c;
                    leaving = Math.Max(leaving, c);
                }
        }
        if (members.Length <= _chain.Size / 2)
            foreach (int i in members) Edges(i, false);
        else
            for (int i = 0; i < _chain.Size; i++)
                if (!isIn(i)) Edges(i, true);
        double logFlow = LogMath.Sum(_buffer.AsSpan(0, count));
        return new SetReading(family, label, members, logMass, logComplement, logFlow, logMass - leaving);
    }

    /// <summary>The merge tree: every connected set formed while adding the edges in descending conductance.</summary>
    public void ReadMergeTree()
    {
        List<(double C, int I, int J)> edges = _chain.EdgesByConductance();
        var label = new int[_chain.Size];
        var groups = new List<int>[_chain.Size];
        for (int i = 0; i < _chain.Size; i++)
        {
            label[i] = i;
            groups[i] = [i];
        }
        foreach ((_, int i, int j) in edges)
        {
            int a = label[i], b = label[j];
            if (a == b) continue;
            if (groups[a].Count < groups[b].Count) (a, b) = (b, a);
            foreach (int k in groups[b]) label[k] = a;
            groups[a].AddRange(groups[b]);
            groups[b] = [];
            if (groups[a].Count == _chain.Size) break;
            int root = a;
            Offer(Measure("merge tree", "", [.. groups[a]], k => label[k] == root));
        }
    }

    /// <summary>Thresholds of the slow eigenvector <paramref name="vector"/> (symmetrized scale) at its eight largest
    /// jumps in f = y/√π, and at its sign.</summary>
    public void ReadEigenvectorCuts(double[] vector)
    {
        double[] f = new double[_chain.Size];
        for (int i = 0; i < f.Length; i++) f[i] = vector[i] * Math.Exp(-0.5 * _chain.LogPi[i]);
        int[] sorted = [.. Enumerable.Range(0, f.Length).OrderBy(i => f[i])];
        IEnumerable<int> jumps = Enumerable.Range(1, f.Length - 1).OrderByDescending(k => f[sorted[k]] - f[sorted[k - 1]]).Take(8);
        int sign = Array.FindIndex(sorted, i => f[i] >= 0.0);
        foreach (int k in jumps.Append(sign).Distinct())
            if (k > 0 && k < f.Length) Read("eigenvector cut", "", sorted[k..]);
    }
}

/// <summary>What the brief reads off a kernel: the gap by the elimination, the second slow mode's, and by the dense
/// solver where asked.</summary>
internal static class BassKernelSpectrum
{
    /// <summary>Above this the dense reading's absolute error near 10⁻¹⁴ leaves it six good digits: the brief's check 3.</summary>
    internal const double DenseFloor = 1e-8;

    /// <summary>The elimination's reading of an irreducible chain: γ and t_rel; the next mode's γ; the eigenvector; the
    /// relative residual. Null gap when t_rel is beyond range.</summary>
    internal sealed record Reading(double? Gap, double? NextGap, double[]? Vector, double Residual, EliminatedChain? Solver);

    internal static Reading Read(ReversibleChain chain, int threads = 1, int seed = 41)
    {
        if (chain.Size < 2) return new Reading(null, null, null, 0.0, null);
        if (EliminatedChain.Of(chain, threads) is not { } solver) return new Reading(null, null, null, double.NaN, null);
        (double gap, double[] vector, double residual) = solver.Slowest([], seed);
        double? next = chain.Size > 2 ? solver.Slowest([vector], seed).Gap : null;
        return new Reading(gap, next, vector, residual, solver);
    }

    /// <summary>The dense reading: λ₂, λ_min, and how far the top eigenpair sits from (1, √π).</summary>
    internal static (double Lambda2, double LambdaMin, double TopOff, double VectorOff) Dense(ReversibleChain chain)
    {
        EigenResult eigen = DenseEigen.DecomposeSymmetric(chain.Symmetrized(), 256, 1e-14);
        double[] v = eigen.Eigenvectors[0];
        double dot = 0.0;
        for (int i = 0; i < chain.Size; i++) dot += v[i] * Math.Exp(0.5 * chain.LogPi[i]);
        return (eigen.Eigenvalues[1], eigen.Eigenvalues[^1], Math.Abs(eigen.Eigenvalues[0] - 1.0), Math.Abs(1.0 - Math.Abs(dot)));
    }

    /// <summary>
    /// Each local mode's depth by the merge tree: when a cluster joins one with a likelier top, ln of its mass less ln of
    /// the joining conductance (0-dimensional persistence on conductance). The deepest at each of <paramref name="slack"/>'s
    /// values, the free slots M_max − M at the cluster's top, with that top; the states with no likelier neighbour; and
    /// each state's depth, NaN at the top state, which joins no likelier one.
    /// </summary>
    internal static (Dictionary<int, (double Depth, int Mode)> Deepest, int LocalModes, double[] DepthOf) ModeDepths(ReversibleChain chain, int[] slack)
    {
        int n = chain.Size;
        var parent = new int[n];
        var top = new int[n];
        double[] logMass = [.. chain.LogPi];
        var depthOf = new double[n];
        Array.Fill(depthOf, double.NaN);
        for (int i = 0; i < n; i++) parent[i] = top[i] = i;
        int Find(int i)
        {
            while (parent[i] != i) i = parent[i] = parent[parent[i]];
            return i;
        }
        var deepest = new Dictionary<int, (double Depth, int Mode)>();
        foreach ((double c, int i, int j) in chain.EdgesByConductance())
        {
            int a = Find(i), b = Find(j);
            if (a == b) continue;
            bool aHolds = chain.LogPi[top[a]] > chain.LogPi[top[b]] || (chain.LogPi[top[a]] == chain.LogPi[top[b]] && top[a] < top[b]);
            int keep = aHolds ? a : b, lose = aHolds ? b : a;
            double depth = logMass[lose] - c;
            depthOf[top[lose]] = depth;
            int s = slack[top[lose]];
            if (!deepest.TryGetValue(s, out (double Depth, int Mode) held) || depth > held.Depth) deepest[s] = (depth, top[lose]);
            parent[lose] = keep;
            logMass[keep] = LogMath.Add(logMass[keep], logMass[lose]);
        }
        int modes = 0;
        for (int i = 0; i < n; i++)
        {
            bool highest = true;
            for (int e = chain.Starts[i]; e < chain.Starts[i + 1] && highest; e++) highest = chain.LogPi[chain.To[e]] <= chain.LogPi[i];
            if (highest) modes++;
        }
        return (deepest, modes, depthOf);
    }

    /// <summary>TV(μ₀K^t, π) at each of <paramref name="times"/> (ascending), and μ_t's mass on <paramref name="set"/> at
    /// them; past a TV below 10⁻⁹ the rest read as there.</summary>
    internal static (double[] Tv, double[] Mass) Law(ReversibleChain chain, double[] start, int[] times, int[]? set)
    {
        double[] pi = chain.Pi();
        double[] mu = [.. start], next = new double[chain.Size];
        var tv = new double[times.Length];
        var mass = new double[times.Length];
        int t = 0;
        double Tv()
        {
            double s = 0.0;
            for (int i = 0; i < mu.Length; i++) s += Math.Abs(mu[i] - pi[i]);
            return 0.5 * s;
        }
        double Mass()
        {
            if (set is null) return double.NaN;
            double s = 0.0;
            foreach (int i in set) s += mu[i];
            return s;
        }
        bool settled = false;
        for (int k = 0; k < times.Length; k++)
        {
            while (t < times[k] && !settled)
            {
                chain.Step(mu, next);
                (mu, next) = (next, mu);
                t++;
                if (t % 500 == 0 && Tv() < 1e-9) settled = true;
            }
            tv[k] = Tv();
            mass[k] = Mass();
        }
        return (tv, mass);
    }

    /// <summary>The brief's scaling reading of y over doubling N: the linear and logarithmic fits, the last three
    /// increments, and "linear", "logarithmic", "flat" or "neither"; "too few" under four rungs.</summary>
    internal static (double B, double RssLinear, double C, double RssLog, double[] Increments, string Reading) Scaling(double[] n, double[] y)
    {
        (double slope, double rss) Fit(double[] x)
        {
            double mx = x.Average(), my = y.Average(), sxx = 0.0, sxy = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                sxx += (x[i] - mx) * (x[i] - mx);
                sxy += (x[i] - mx) * (y[i] - my);
            }
            double b = sxx > 0.0 ? sxy / sxx : 0.0, rss = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                double r = y[i] - my - b * (x[i] - mx);
                rss += r * r;
            }
            return (b, rss);
        }
        (double b, double rssLinear) = Fit(n);
        (double c, double rssLog) = Fit([.. n.Select(v => Math.Log(v))]);
        double[] delta = [.. Enumerable.Range(1, y.Length - 1).Select(k => y[k] - y[k - 1])];
        if (delta.Length < 3) return (b, rssLinear, c, rssLog, delta, "too few");
        double[] last = delta[^3..];
        string reading = "neither";
        if (last.All(d => Math.Abs(d) < 0.1)) reading = "flat";
        else if (last.All(d => d > 0.0) && last[1] >= 1.5 * last[0] && last[2] >= 1.5 * last[1]) reading = "linear";
        else if (last.All(d => d > 0.0) && last.Max() <= 1.25 * last.Min()) reading = "logarithmic";
        return (b, rssLinear, c, rssLog, delta, reading);
    }
}
