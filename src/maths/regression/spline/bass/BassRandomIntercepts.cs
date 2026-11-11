using System;
using System.Collections.Generic;
using System.Text.Json.Serialization;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>How a random-intercept block keeps clear of the level-centred intercept (E21).</summary>
public enum RandomInterceptConvention
{
    /// <summary>
    /// The group effects' observation-weighted mean is zero, Σ_g n_g u_g = 0, so their contribution Zu has mean zero over the cases and
    /// the intercept alone carries the response's level. The unweighted sum-to-zero constraint is not offered: on unbalanced groups Zu
    /// keeps a nonzero mean, which is the overlap with the level-centred intercept the convention exists to remove.
    /// </summary>
    ObservationCentred,

    /// <summary>The observation-weighted effect mean is zero separately in each parent level.</summary>
    WithinParent,
}

/// <summary>
/// s's prior when it is drawn (E21's slice 1; <c>ridge-blocks-brief.md</c>, decision 2): the half-Cauchy C⁺(0, A) on s, A =
/// <paramref name="Scale"/> in units of σ, 1 by default, a group effect about as large as the noise (Gelman 2006; Polson and Scott 2012).
/// It is proper, so a reference rung and SMC's start draw s from it.
/// </summary>
/// <param name="Scale">A, the half-Cauchy's scale, positive and finite.</param>
public sealed record BassGroupScalePrior(double Scale = 1.0)
{
    /// <summary>ln p(s) for s &gt; 0: ln(2/(πA)) − ln(1 + (s/A)²); −∞ elsewhere.</summary>
    public double LogDensity(double s)
    {
        if (!(s > 0.0) || !double.IsFinite(s)) return double.NegativeInfinity;
        double r = s <= Scale ? s / Scale : Scale / s;
        double logScale = Math.Log(Scale);
        return Math.Log(2.0 / Math.PI) + (s <= Scale ? -logScale : logScale - 2.0 * Math.Log(s)) - SpecialFunctions.Log1p(r * r);
    }

    /// <summary>A draw of s: A·tan(πu/2), u uniform on (0, 1).</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(Scale > 0.0) || !double.IsFinite(Scale)) throw new ArgumentOutOfRangeException(nameof(Scale));
        double s;
        do s = Scale * Math.Tan(0.5 * Math.PI * rng.NextDouble()); while (!(s > 0.0) || !double.IsFinite(s));
        return s;
    }

    // −∂²/∂t² ln p(eᵗ) at t = ln s: 4r/(1 + r)², r = (s/A)².
    internal double CurvatureInLog(double s)
    {
        double ratio = s <= Scale ? s / Scale : Scale / s;
        double r = ratio * ratio;
        return 4.0 * r / ((1.0 + r) * (1.0 + r));
    }
}

/// <summary>
/// Random-effect intercepts on a named data grouping (E21; <see cref="BassData.Groupings"/>): y = α + Σⱼ βⱼBⱼ(x) + Σ_k u_kg + ε, with
/// each u_k | σ² ~ N(0, σ²s_k²Π_k) on its observed levels, Π_k the projector onto its declared constraints (<see cref="Convention"/>), so
/// <paramref name="Scale"/> s is the group effects' standard deviation in units of σ: held, or drawn under <paramref name="ScalePrior"/>
/// (<see cref="BassGroupScaleBlock"/>). No term reads the grouping, and the coefficient prior is the config's, of any form. The block is
/// jointly integrated out of every structure's marginal (<see cref="RandomEffectsBlock"/>), and each retained draw carries its named
/// groups' scales and effects (<see cref="BassDraw.Groups"/>).
/// </summary>
/// <param name="Grouping">The data grouping's unique name.</param>
/// <param name="Scale">s: the group effects' standard deviation over σ, positive and finite; held, or the chain's start when it is drawn.</param>
/// <param name="Convention">Observation centring, or centring within each level of a declared parent.</param>
/// <param name="ScalePrior">s's prior, which draws it (E21's slice 1); null, the default, holds s. Left out of the serialized form when
/// null.</param>
/// <param name="ParentGrouping">Within-parent centring's parent data grouping; absent under observation centring.</param>
public sealed record BassRandomIntercept(string Grouping, double Scale, RandomInterceptConvention Convention = RandomInterceptConvention.ObservationCentred,
                                         [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] BassGroupScalePrior? ScalePrior = null,
                                         [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? ParentGrouping = null) : IBassRandomEffect;

/// <summary>
/// A random-intercept block on one data set (E21), integrated out. With Z the cases' level indicators over the L observed levels, D = ZᵀZ
/// their counts and H an orthonormal basis of {u : nᵀu = 0}, u = Hγ with γ ~ N(0, σ²s²I), so y has covariance σ²V,
/// V = I + s²(ZH)(ZH)ᵀ. Since nᵀH = 0, V1 = 1: the constant column is V's eigenvector with eigenvalue 1, so V^(−½) leaves the intercept's
/// column, every column's mean and ȳ as they are.
/// <list type="bullet">
/// <item><b>The marginal.</b> Whitening by V^(−½) turns the model into BASS's own, y_w = α + X_wβ + ε with ε ~ N(0, σ²I), at the price
/// −½ln|V| on the log density, which the working response carries as weights carry ½Σ ln wᵢ (<see cref="WorkingResponse"/>). With
/// HᵀDH = EΛEᵀ, V^(−½)x = x + Z·M·Zᵀx, M = HE·diag(((1 + s²λ)^(−½) − 1)/λ)·EᵀHᵀ: O(n + L²) a column, and ln|V| = Σ ln(1 + s²λ).
/// Every coefficient form then reads the whitened design as it reads a design.</item>
/// <item><b>The group effects.</b> Given α, β and σ², u's conditional is exact: with r = y − α − Xβ on the unwhitened columns,
/// γ ~ N((HᵀDH + I/s²)⁻¹HᵀZᵀr, σ²(HᵀDH + I/s²)⁻¹). A composed draw of (α, β, σ²) from the whitened model, u integrated, then u from
/// this conditional, is a joint draw.</item>
/// </list>
/// </summary>
public sealed class RandomInterceptBlock
{
    private readonly int[] _group;          // each case's observed-level index
    private readonly int[] _observedOf;     // each declared level's observed index, −1 for a level no case takes
    private readonly double[] _counts;      // n_g over the observed levels
    private readonly double[,] _he;         // HE, L × (L − 1): HᵀDH's eigenvectors carried back to the levels, free of s
    private readonly double[] _eigenvalues; // Λ, HᵀDH's spectrum, free of s
    private readonly double[,] _whiten;     // M, L × L
    private readonly double[,] _mean;       // HE·diag(1/(λ + 1/s²))·EᵀHᵀ, L × L
    private readonly double[,] _spread;     // HE·diag((λ + 1/s²)^(−½)), L × (L − 1)

    // The block at spec's s on the grouping's basis and spectrum, which do not depend on s.
    private RandomInterceptBlock(BassRandomIntercept spec, int[] group, int[] observedOf, double[] counts, double[,] he, double[] eigenvalues)
    {
        Spec = spec;
        LogScalePrior = spec.ScalePrior is { } prior ? prior.LogDensity(spec.Scale) + Math.Log(spec.Scale) : 0.0;
        _group = group;
        _observedOf = observedOf;
        _counts = counts;
        _he = he;
        _eigenvalues = eigenvalues;
        int levels = counts.Length, q = eigenvalues.Length;
        double s2 = spec.Scale * spec.Scale, logDeterminant = 0.0;
        var whitenScale = new double[q];
        var meanScale = new double[q];
        var spread = new double[levels, q];
        for (int k = 0; k < q; k++)
        {
            double lambda = eigenvalues[k];
            logDeterminant += double.IsPositiveInfinity(s2 * lambda)
                ? 2.0 * Math.Log(spec.Scale) + Math.Log(lambda) : Math.Log(1.0 + s2 * lambda);
            whitenScale[k] = (1.0 / Math.Sqrt(1.0 + s2 * lambda) - 1.0) / lambda;
            double precision = lambda + 1.0 / s2;
            meanScale[k] = 1.0 / precision;
            for (int g = 0; g < levels; g++) spread[g, k] = he[g, k] / Math.Sqrt(precision);
        }
        _whiten = Outer(he, whitenScale);
        _mean = Outer(he, meanScale);
        _spread = spread;
        LogDeterminant = logDeterminant;
    }

    /// <summary>The block <paramref name="spec"/> declares on <paramref name="grouping"/>'s cases.</summary>
    public static RandomInterceptBlock Of(CategoricalPredictor grouping, BassRandomIntercept spec)
    {
        ArgumentNullException.ThrowIfNull(grouping);
        ArgumentNullException.ThrowIfNull(spec);
        if (spec.Grouping != grouping.Name) throw new ArgumentException("The declaration must name this grouping.", nameof(spec));
        if (!(spec.Scale > 0.0) || !double.IsFinite(spec.Scale)) throw new ArgumentOutOfRangeException(nameof(spec), "The scale is positive and finite.");
        if (spec.Convention != RandomInterceptConvention.ObservationCentred) throw new ArgumentOutOfRangeException(nameof(spec), "The convention is not built.");
        ReadOnlySpan<int> codes = grouping.Codes.Span;
        var observedOf = new int[grouping.Levels.Count];
        Array.Fill(observedOf, -1);
        int levels = 0;
        for (int l = 0; l < observedOf.Length; l++)
            if (grouping.Observed[l]) observedOf[l] = levels++;
        var group = new int[codes.Length];
        var counts = new double[levels];
        for (int i = 0; i < codes.Length; i++)
        {
            group[i] = observedOf[codes[i]];
            counts[group[i]]++;
        }

        // H: an orthonormal basis of n's complement, by Gram–Schmidt over the unit vectors with n̂ removed, twice for its digits.
        double norm = 0.0;
        foreach (double c in counts) norm += c * c;
        norm = Math.Sqrt(norm);
        var basis = new List<double[]>(levels - 1);
        for (int k = 0; k < levels && basis.Count < levels - 1; k++)
        {
            var v = new double[levels];
            v[k] = 1.0;
            for (int pass = 0; pass < 2; pass++)
            {
                Remove(v, counts, norm);
                foreach (double[] b in basis) Remove(v, b, 1.0);
            }
            double length = Math.Sqrt(Dot(v, v));
            if (length < 1e-8) continue;
            for (int g = 0; g < levels; g++) v[g] /= length;
            basis.Add(v);
        }
        int q = basis.Count;   // L − 1

        // HᵀDH = EΛEᵀ.
        var gram = new double[q, q];
        for (int a = 0; a < q; a++)
            for (int b = 0; b <= a; b++)
            {
                double sum = 0.0;
                for (int g = 0; g < levels; g++) sum += basis[a][g] * counts[g] * basis[b][g];
                gram[a, b] = gram[b, a] = sum;
            }
        EigenResult eigen = DenseEigen.DecomposeSymmetric(gram, tol: 1e-15);
        // HE: each eigenvector carried back to the levels.
        var he = new double[levels, q];
        for (int k = 0; k < q; k++)
            for (int g = 0; g < levels; g++)
            {
                double sum = 0.0;
                for (int a = 0; a < q; a++) sum += basis[a][g] * eigen.Eigenvectors[k][a];
                he[g, k] = sum;
            }
        var eigenvalues = new double[q];
        for (int k = 0; k < q; k++)
        {
            eigenvalues[k] = eigen.Eigenvalues[k];
            if (!(eigenvalues[k] > 0.0)) throw new InvalidOperationException("The groups' Gram on n's complement is not positive definite.");
        }
        return new RandomInterceptBlock(spec, group, observedOf, counts, he, eigenvalues);
    }

    /// <summary>This block at scale <paramref name="scale"/>: the grouping's basis and spectrum shared, its whitening, its effects'
    /// conditional and ln|V| at s = <paramref name="scale"/>; this block itself at its own s. O(L²·(L − 1)).</summary>
    public RandomInterceptBlock At(double scale)
    {
        if (!(scale > 0.0) || !double.IsFinite(scale)) throw new ArgumentOutOfRangeException(nameof(scale), "The scale is positive and finite.");
        return scale == Spec.Scale ? this : new RandomInterceptBlock(Spec with { Scale = scale }, _group, _observedOf, _counts, _he, _eigenvalues);
    }

    /// <summary>The expected information about ln s at this block's scale, ½tr((V⁻¹∂V/∂ln s)²) = 2Σ_k (s²λ_k/(1 + s²λ_k))², λ the levels'
    /// spectrum on n's complement: what the scale's step reads for its proposal (<see cref="BassGroupScaleBlock"/>).</summary>
    internal double LogScaleInformation
    {
        get
        {
            double s2 = Spec.Scale * Spec.Scale, sum = 0.0;
            foreach (double lambda in _eigenvalues)
            {
                double x = s2 * lambda;
                double share = double.IsPositiveInfinity(x) ? 1.0 : x / (1.0 + x);
                sum += share * share;
            }
            return 2.0 * sum;
        }
    }

    private static void Remove(double[] v, double[] direction, double length)
    {
        double along = Dot(v, direction) / (length * length);
        for (int g = 0; g < v.Length; g++) v[g] -= along * direction[g];
    }

    private static double Dot(double[] a, double[] b)
    {
        double sum = 0.0;
        for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
        return sum;
    }

    // A·diag(d)·Aᵀ for A L × q.
    private static double[,] Outer(double[,] a, double[] d)
    {
        int l = a.GetLength(0), q = a.GetLength(1);
        var m = new double[l, l];
        for (int g = 0; g < l; g++)
            for (int h = 0; h <= g; h++)
            {
                double sum = 0.0;
                for (int k = 0; k < q; k++) sum += a[g, k] * d[k] * a[h, k];
                m[g, h] = m[h, g] = sum;
            }
        return m;
    }

    /// <summary>The declaration.</summary>
    public BassRandomIntercept Spec { get; }

    // Density in the sampled coordinate ln s, including its Jacobian. Held s contributes no prior factor.
    internal double LogScalePrior { get; }

    /// <summary>L, the observed levels.</summary>
    public int Levels => _counts.Length;

    /// <summary>Each observed level's case count, n_g.</summary>
    public ReadOnlySpan<double> Counts => _counts;

    /// <summary>ln|V|, the block's price on the log density: −½ln|V| at b = 1.</summary>
    public double LogDeterminant { get; }

    /// <summary>Case <paramref name="i"/>'s observed level.</summary>
    public int GroupOf(int i) => _group[i];

    /// <summary>The observed level of declared level <paramref name="code"/>; −1 for a level no training case takes.</summary>
    public int ObservedOf(int code) => (uint)code < (uint)_observedOf.Length ? _observedOf[code] : -1;

    /// <summary>Zᵀx: <paramref name="x"/>'s sums over each observed level's cases.</summary>
    public double[] GroupSums(ReadOnlySpan<double> x)
    {
        if (x.Length != _group.Length) throw new ArgumentException($"The column holds {x.Length} values for {_group.Length} cases.", nameof(x));
        var sums = new double[_counts.Length];
        for (int i = 0; i < x.Length; i++) sums[_group[i]] += x[i];
        return sums;
    }

    /// <summary>V^(−½)x, in place: x + Z·M·Zᵀx.</summary>
    public void Whiten(Span<double> x)
    {
        double[] sums = GroupSums(x);
        int l = _counts.Length;
        var shift = new double[l];
        for (int g = 0; g < l; g++)
        {
            double sum = 0.0;
            for (int h = 0; h < l; h++) sum += _whiten[g, h] * sums[h];
            shift[g] = sum;
        }
        for (int i = 0; i < x.Length; i++) x[i] += shift[_group[i]];
    }

    /// <summary>
    /// The group effects u over the observed levels from their conditional given α, β and σ² (<paramref name="sigma2"/>), with
    /// <paramref name="residualSums"/> = Zᵀ(y − α − Xβ) on the unwhitened columns: u = H·γ, γ's mean and spread as the block holds them.
    /// Σ_g n_g u_g = 0.
    /// </summary>
    public double[] DrawEffects(ReadOnlySpan<double> residualSums, double sigma2, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int l = _counts.Length, q = _spread.GetLength(1);
        if (residualSums.Length != l) throw new ArgumentException($"The sums hold {residualSums.Length} values for {l} levels.", nameof(residualSums));
        if (!(sigma2 > 0.0) || !double.IsFinite(sigma2)) throw new ArgumentOutOfRangeException(nameof(sigma2));
        var z = new double[q];
        for (int k = 0; k < q; k++) z[k] = Normal.SampleStandard(rng);
        double sigma = Math.Sqrt(sigma2);
        var u = new double[l];
        for (int g = 0; g < l; g++)
        {
            double mean = 0.0, noise = 0.0;
            for (int h = 0; h < l; h++) mean += _mean[g, h] * residualSums[h];
            for (int k = 0; k < q; k++) noise += _spread[g, k] * z[k];
            u[g] = mean + sigma * noise;
        }
        return u;
    }

    /// <summary>
    /// Zᵀr, r = y − α − Σⱼ βⱼBⱼ, on <paramref name="model"/>'s observations and the unwhitened case columns of <paramref name="terms"/>, with
    /// <paramref name="coefficients"/> = [α, β₁, …] in their order: the sums <see cref="DrawEffects"/> and <see cref="EffectMean"/> read.
    /// </summary>
    internal double[] ResidualSums(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients)
    {
        if (coefficients.Length != terms.Count + 1)
            throw new ArgumentException($"Give the intercept and one coefficient per term ({terms.Count + 1}).", nameof(coefficients));
        double[] sums = GroupSums(model.Data.Response.Values);
        for (int g = 0; g < sums.Length; g++) sums[g] -= coefficients[0] * _counts[g];
        var column = new double[_group.Length];
        for (int k = 0; k < terms.Count; k++)
        {
            model.Basis.CaseColumn(terms[k], model.Inputs.Cases, column);
            double[] termSums = GroupSums(column);
            for (int g = 0; g < sums.Length; g++) sums[g] -= coefficients[k + 1] * termSums[g];
        }
        return sums;
    }

    /// <summary>The group effects' conditional mean given Zᵀr = <paramref name="residualSums"/>: <see cref="DrawEffects"/>'s centre.</summary>
    public double[] EffectMean(ReadOnlySpan<double> residualSums)
    {
        int l = _counts.Length;
        var u = new double[l];
        for (int g = 0; g < l; g++)
            for (int h = 0; h < l; h++) u[g] += _mean[g, h] * residualSums[h];
        return u;
    }
}

/// <summary>
/// s's step when it is drawn (E21's slice 1; <see cref="BassRandomIntercept.ScalePrior"/>), in τ's slot after τ's: a Metropolis step on
/// t = ln s against the marginal with a, σ² and u integrated, at the chain's temperature, by a refit at the proposal
/// (<see cref="BassStructureMoves.ProposeFollow"/>) on the reading whitened at s′ = s·eᵗ.
/// The joint target includes the half-Cauchy's density and Jacobian in ln s.
/// <para>ln A = Δ ln π_b + ln q(−t | x′) − ln q(t | x),</para>
/// <para>t ~ N(0, h(x)²), h = 2.4/√C(x), C the expected information about ln s (<see cref="RandomInterceptBlock.LogScaleInformation"/>) plus
/// the prior's curvature in ln s, at least 1. It reads the state's s and the levels alone, so its Hastings ratio is in A and the step is
/// exact. Nothing reads the group effects between this step and their draw at composition, which takes the draw's s. A move whose
/// structure leaves F is declined, as π is 0 there.</para>
/// </summary>
public sealed class BassGroupScaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private const double Step = 2.4;   // the random walk's scale over the target's sd along t
    private readonly BassGroupScalePrior _prior;
    private readonly int _group;
    private readonly string _name;

    internal BassGroupScaleBlock(BassModel model, int group = 0)
    {
        ArgumentNullException.ThrowIfNull(model);
        RandomEffectsBlock blocks = model.RandomEffects
            ?? throw new ArgumentException("The group scale's step needs random intercepts.", nameof(model));
        if ((uint)group >= (uint)blocks.Count) throw new ArgumentOutOfRangeException(nameof(group));
        _group = group;
        _prior = blocks.Specifications[group].ScalePrior
                 ?? throw new ArgumentException("The group scale's step needs a drawn scale (RandomIntercepts[].ScalePrior).", nameof(model));
        _name = blocks.Count == 1 ? "group scale" : $"group scale {blocks.Specifications[group].Grouping}"
            + (blocks.Layout(group).Input is { } input ? $" on {input}" : "");
    }

    /// <inheritdoc/>
    public string Name => _name;

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        RandomEffectsBlock block = Reading(state).Whitening ?? throw new InvalidOperationException("The chain's reading carries no block.");
        double h = Scale(block), t = h * Normal.SampleStandard(rng);
        if (Propose(state, block, h, t) is not { } move) return false;
        if (!(Math.Log(1.0 - rng.NextDouble()) < move.LogAcceptance)) return false;
        state.SetState(move.Candidate, move.Evaluation);
        return true;
    }

    /// <summary>A proposed step: the candidate state on the reading at s′, its evaluation, and ln A.</summary>
    internal readonly record struct Move(BassState Candidate, ITargetEvaluation Evaluation, double LogAcceptance);

    // The reading the chain's state is on.
    internal static WorkingResponse Reading(ReversibleJumpChain<BassState> state) => BassStructureMoves.Of(state).DesignAt(state.Current).Working;

    // h = Step/√C at the block's s.
    internal double Scale(RandomEffectsBlock block)
    {
        double c = block.ScaleInformation(_group) + _prior.CurvatureInLog(block.Scales.Span[_group]);
        return Step / Math.Sqrt(double.IsFinite(c) ? Math.Max(c, 1.0) : 1.0);
    }

    // The step by t from the block's s, or null when s′ is not a scale or the structure leaves F there.
    internal Move? Propose(ReversibleJumpChain<BassState> state, RandomEffectsBlock block, double h, double t)
    {
        double s = block.Scales.Span[_group], moved = s * Math.Exp(t);
        if (!(moved > 0.0) || !double.IsFinite(moved)) return null;
        RandomEffectsBlock at = block.At(_group, moved);
        double reverse = Scale(at);
        double proposalRatio = Math.Log(h) - Math.Log(reverse) - t * t / (2.0 * reverse * reverse) + t * t / (2.0 * h * h);
        BassState candidate = BassStructureMoves.Of(state).ProposeFollow(state.Current, at.Rewhiten(Reading(state)));
        if ((state.Target as BassTarget ?? throw new InvalidOperationException("The chain's target is BASS's.")).Refusal(candidate.Fit) is not null) return null;
        ITargetEvaluation evaluation = state.Target.Evaluate(candidate);
        double target = state.Temperature.LogDensity(evaluation) - state.CurrentLogDensity;
        return new Move(candidate, evaluation, target + proposalRatio);
    }
}
