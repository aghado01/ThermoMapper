using System;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Birth move: draw a new knot from the birth density <see cref="KnotBirthDensity"/>, which perturbs an existing
/// knot with the proposal <see cref="IKnotKernel"/>, draws from its base with probability ε, and draws from the base
/// when k = 0. With equal birth/death weight the move-selection probabilities cancel; the Hastings correction is
/// <c>−log(k+1) − log q(cand | knots)</c>, and the Jacobian is 0. With <see cref="UniformKernel"/> and the uniform base
/// this reduces to <c>−log(k+1)</c>.
/// </summary>
public sealed class KnotBirthMove : IRjMove<KnotConfig>
{
    private readonly KnotBirthDensity _density;
    private readonly double _weight;

    /// <summary>Local births from <paramref name="kernel"/> (ε = 0), uniform when there are no knots.</summary>
    public KnotBirthMove(IKnotKernel kernel, double weight = 0.4) : this(new KnotBirthDensity(kernel), weight) { }

    /// <summary>Births from <paramref name="density"/>; its <see cref="KnotDeathMove"/> must share it.</summary>
    public KnotBirthMove(KnotBirthDensity density, double weight = 0.4)
    {
        ArgumentNullException.ThrowIfNull(density);
        if (!(weight > 0.0)) throw new ArgumentOutOfRangeException(nameof(weight));
        _density = density;
        _weight = weight;
    }

    public string Key => "birth";
    public string ReverseKey => "death";
    public double Weight(KnotConfig state) => _weight;

    public Proposal<KnotConfig>? Propose(KnotConfig current, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(current);
        ArgumentNullException.ThrowIfNull(rng);

        double[] knots = current.InteriorKnots;
        int k = knots.Length;
        double candidate = _density.Sample(knots, rng);

        var next = new double[k + 1];
        int j = 0;
        while (j < k && knots[j] < candidate) { next[j] = knots[j]; j++; }
        next[j] = candidate;
        for (int t = j; t < k; t++) next[t + 1] = knots[t];

        double logForward = _density.LogDensity(candidate, knots);
        return new Proposal<KnotConfig>(new KnotConfig(next), -Math.Log(k + 1) - logForward, 0.0);
    }
}

/// <summary>
/// Death move: delete one interior knot chosen uniformly. Exact reverse of <see cref="KnotBirthMove"/> with the same
/// <see cref="KnotBirthDensity"/> — Hastings correction <c>+log(k) + log q(dropped | remaining)</c>, Jacobian 0.
/// Unavailable (null) at k = 0.
/// </summary>
public sealed class KnotDeathMove : IRjMove<KnotConfig>
{
    private readonly KnotBirthDensity _density;
    private readonly double _weight;

    /// <summary>The reverse of local births from <paramref name="kernel"/> (ε = 0).</summary>
    public KnotDeathMove(IKnotKernel kernel, double weight = 0.4) : this(new KnotBirthDensity(kernel), weight) { }

    /// <summary>The reverse of births from <paramref name="density"/>, which its <see cref="KnotBirthMove"/> shares.</summary>
    public KnotDeathMove(KnotBirthDensity density, double weight = 0.4)
    {
        ArgumentNullException.ThrowIfNull(density);
        if (!(weight > 0.0)) throw new ArgumentOutOfRangeException(nameof(weight));
        _density = density;
        _weight = weight;
    }

    public string Key => "death";
    public string ReverseKey => "birth";
    public double Weight(KnotConfig state) => _weight;

    public Proposal<KnotConfig>? Propose(KnotConfig current, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(current);
        ArgumentNullException.ThrowIfNull(rng);

        double[] knots = current.InteriorKnots;
        int k = knots.Length;
        if (k == 0) return null;

        int drop = rng.NextInt(k);
        double dropped = knots[drop];
        var next = new double[k - 1];
        for (int t = 0, w = 0; t < k; t++)
            if (t != drop) next[w++] = knots[t];

        double logReverse = _density.LogDensity(dropped, next);
        return new Proposal<KnotConfig>(new KnotConfig(next), Math.Log(k) + logReverse, 0.0);
    }
}

/// <summary>
/// Relocate move: move one knot (chosen uniformly) by perturbing it with the proposal kernel. The Hastings
/// correction is the asymmetric <c>log h(old | new) − log h(new | old)</c> (0 for the symmetric
/// <see cref="UniformKernel"/>); Jacobian 0. Unavailable (null) at k = 0.
/// </summary>
public sealed class KnotRelocateMove : IRjMove<KnotConfig>
{
    private readonly IKnotKernel _kernel;
    private readonly double _weight;

    public KnotRelocateMove(IKnotKernel kernel, double weight = 0.2)
    {
        ArgumentNullException.ThrowIfNull(kernel);
        if (!(weight > 0.0)) throw new ArgumentOutOfRangeException(nameof(weight));
        _kernel = kernel;
        _weight = weight;
    }

    public string Key => "relocate";
    public string ReverseKey => "relocate";
    public double Weight(KnotConfig state) => _weight;

    public Proposal<KnotConfig>? Propose(KnotConfig current, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(current);
        ArgumentNullException.ThrowIfNull(rng);

        double[] knots = current.InteriorKnots;
        int k = knots.Length;
        if (k == 0) return null;

        int idx = rng.NextInt(k);
        double oldKnot = knots[idx];
        double newKnot = _kernel.Sample(oldKnot, rng);

        var next = new double[k];
        Array.Copy(knots, next, k);
        next[idx] = newKnot;
        Array.Sort(next);

        double logRatio = _kernel.LogDensity(oldKnot, newKnot) - _kernel.LogDensity(newKnot, oldKnot);
        return new Proposal<KnotConfig>(new KnotConfig(next), logRatio, 0.0);
    }
}
