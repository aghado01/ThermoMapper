using System;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Birth move: pick an existing knot uniformly as a center, perturb it with the proposal
/// <see cref="IKnotKernel"/> to get a new knot (uniform on (0,1) when k = 0). With equal birth/death weight the
/// move-selection probabilities cancel; the Hastings correction is <c>−log(k+1) − log q(cand | knots)</c> where
/// q is the mixture density, and the Jacobian is 0. With <see cref="UniformKernel"/> this reduces to
/// <c>−log(k+1)</c>.
/// </summary>
public sealed class KnotBirthMove : IRjMove<KnotConfig>
{
    private readonly IKnotKernel _kernel;
    private readonly double _weight;

    public KnotBirthMove(IKnotKernel kernel, double weight = 0.4)
    {
        ArgumentNullException.ThrowIfNull(kernel);
        if (!(weight > 0.0)) throw new ArgumentOutOfRangeException(nameof(weight));
        _kernel = kernel;
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
        double candidate = k == 0
            ? rng.NextDouble()
            : _kernel.Sample(knots[rng.NextInt(k)], rng);

        var next = new double[k + 1];
        int j = 0;
        while (j < k && knots[j] < candidate) { next[j] = knots[j]; j++; }
        next[j] = candidate;
        for (int t = j; t < k; t++) next[t + 1] = knots[t];

        double logForward = ProposalMath.LogMixtureDensity(candidate, knots, _kernel);
        return new Proposal<KnotConfig>(new KnotConfig(next), -Math.Log(k + 1) - logForward, 0.0);
    }
}

/// <summary>
/// Death move: delete one interior knot chosen uniformly. Exact reverse of <see cref="KnotBirthMove"/> —
/// Hastings correction <c>+log(k) + log q(dropped | remaining)</c>, Jacobian 0. Unavailable (null) at k = 0.
/// </summary>
public sealed class KnotDeathMove : IRjMove<KnotConfig>
{
    private readonly IKnotKernel _kernel;
    private readonly double _weight;

    public KnotDeathMove(IKnotKernel kernel, double weight = 0.4)
    {
        ArgumentNullException.ThrowIfNull(kernel);
        if (!(weight > 0.0)) throw new ArgumentOutOfRangeException(nameof(weight));
        _kernel = kernel;
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

        double logReverse = ProposalMath.LogMixtureDensity(dropped, next, _kernel);
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
