using System;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// F-X4, the proposal-kernel contract: a structure proposal over a <see cref="BassStructurePrior"/>'s term space that states
/// its forward and reverse probabilities exactly, given the current structure's <see cref="BassTermCounts"/>, so a
/// reversible-jump chain's acceptance holds the law. A birth's term with its exact log probability; a death's ratio, the
/// removed term's birth under the reduced counts against the death's uniform pick; a change of one term with its within-move
/// log ratio. The chain's sweep (<see cref="BassStructureMoves"/>) reads its moves through this contract. Members: D01's
/// hinge proposal (<see cref="BassTermProposal"/>), its draws unchanged, and KHAOS's Bernoulli-inclusion proposal over a
/// Legendre term space (<see cref="LegendreTermProposal"/>).
/// </summary>
public interface IBassProposalKernel
{
    /// <summary>The prior whose term space the proposals range over.</summary>
    BassStructurePrior Prior { get; }

    /// <summary>A birth's term into the structure <paramref name="counts"/> describes; null for a draw that proposes nothing (D01's
    /// intercept).</summary>
    BassTerm? DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts);

    /// <summary>ln q_b(term | counts), the exact log probability that a birth draws <paramref name="term"/>.</summary>
    double LogBirth(BassTerm term, BassTermCounts counts);

    /// <summary>A birth's within-move log ratio into M terms: ln[1/(M + 1)] − ln q_b(born | counts).</summary>
    double LogBirthRatio(BassTerm born, BassTermCounts counts);

    /// <summary>A death's within-move log ratio from M terms: ln q_b(dying | the counts without it) − ln(1/M).</summary>
    double LogDeathRatio(BassTerm dying, BassTermCounts counts);

    /// <summary>A change of <paramref name="term"/>, a term of the structure <paramref name="counts"/> describes, with its
    /// within-move log ratio ln q(reverse)/q(forward).</summary>
    BassTermMove ProposeChange(Xoshiro256PlusPlus rng, BassTerm term, BassTermCounts counts);
}

/// <summary>
/// KHAOS's proposal over a Legendre term space (E4, E19; the KHAOS brief's slice 1), the corrected reading of <c>khaos</c> 2.1.1
/// that slice 0 checked against the exact law (<see cref="LegendreBasis"/> for the settings):
/// <list type="bullet">
/// <item><b>The birth.</b> q₀ ∝ q₀^(−s_q) on {1, …, q_max}; η = <c>make_weights</c>(appearances, q₀); χⱼ ~ Bern(ηⱼ), redrawn while
/// q = Σχ is 0 or above q_max; d ∝ (d − q + 1)^(−s_d) on {q, …, d_max}; a uniform composition of d over the chosen variables.
/// q₀ is drawn once, before the redraws, so the inclusion's law is Σ_{q₀} J(q₀)·P(χ | q₀)/P(1 ≤ q ≤ q_max | q₀), each validity
/// the Poisson-binomial's (<see cref="PoissonBinomial"/>); with the degree's and the composition's 1/C(d − 1, q − 1) that is the
/// forward probability, exactly. <c>khaos</c> multiplies every redraw's probability in instead, and approximates the death's
/// validity by a normal; slice 0 confirmed both.</item>
/// <item><b>The death.</b> A term uniformly; its reverse is the same expression under the counts without it.</item>
/// <item><b>The change.</b> With probability <see cref="LegendreBasis.DegreeMutation"/> (1 when p ≤ 3) the degree is redrawn
/// from the birth's law at the term's order with a new uniform composition; otherwise one variable, picked uniformly, is
/// swapped for one the term does not use, ∝ w + its appearances, its degree kept. A swap onto a variable the term holds would
/// repeat a variable, outside the basis (slice 0 confirmed <c>khaos</c> draws it), so those are excluded on both sides. A term
/// that holds every variable has no swap: a null move.</item>
/// </list>
/// Immutable, with no scratch state, so one serves every chain of a fit.
/// </summary>
public sealed class LegendreTermProposal : IBassProposalKernel, IBassWeightedBirth
{
    private const int StackLimit = 256;

    private readonly double[] _logOrder;     // ln J(q₀), normalized, by q₀ (index 0 unused)
    private readonly double[][] _logDegree;  // ln D(d | q), normalized, by q then d (indices below q unused)

    /// <summary>KHAOS's proposal over <paramref name="prior"/>'s Legendre term space with <paramref name="settings"/>.</summary>
    public LegendreTermProposal(BassStructurePrior prior, LegendreBasis settings)
    {
        ArgumentNullException.ThrowIfNull(prior);
        ArgumentNullException.ThrowIfNull(settings);
        if (prior.Legendre is not { } space) throw new ArgumentException("KHAOS's proposal ranges over a Legendre term space.", nameof(prior));
        if (settings.Problem() is { } problem) throw new ArgumentOutOfRangeException(nameof(settings), $"{problem.Setting} {problem.Reason}.");
        if (settings.MaxDegree != space.MaxDegree) throw new ArgumentException("The settings' d_max is not the term space's.", nameof(settings));
        Prior = prior;
        Space = space;
        Settings = settings;
        int top = space.MaxOrder, d = space.MaxDegree;
        _logOrder = Normalized(1, top, q0 => -settings.OrderPenalty * Math.Log(q0), top + 1);
        _logDegree = new double[top + 1][];
        for (int q = 1; q <= top; q++)
        {
            int order = q;
            _logDegree[q] = Normalized(q, d, degree => -settings.DegreePenalty * Math.Log(degree - order + 1), d + 1);
        }
    }

    /// <inheritdoc/>
    public BassStructurePrior Prior { get; }

    /// <summary>The term space.</summary>
    public LegendreTermSpace Space { get; }

    /// <summary>The proposal's settings.</summary>
    public LegendreBasis Settings { get; }

    /// <summary>The probability a change mutates the degree: 1 when p ≤ 3, KHAOS's rule, the setting otherwise.</summary>
    public double DegreeMutation => Space.Variables <= 3 ? 1.0 : Settings.DegreeMutation;

    /// <summary>
    /// <c>make_weights</c> as <c>khaos</c> 2.1.1 writes it: vⱼ = (aⱼ^α + ε)/Σ(a^α + ε)·q₀/p from the appearance counts a, then
    /// <see cref="LegendreBasis.InclusionPasses"/> steps of β = mean((ln(q₀ + δ) − ln p)/ln v), δ += q₀ − Σv^β, giving ηⱼ = vⱼ^β,
    /// whose sum is about q₀. One variable is always in (η = 1); a value past 1 is held there, though none arises with p ≥ 2,
    /// where every vⱼ &lt; 1. The proposal reads the same η both ways, so their sum need not be q₀ exactly.
    /// </summary>
    public void InclusionProbabilities(ReadOnlySpan<int> appearances, int q0, Span<double> eta)
    {
        int p = appearances.Length;
        if (eta.Length != p) throw new ArgumentException("One probability per variable.", nameof(eta));
        for (int j = 0; j < p; j++) eta[j] = Math.Pow(appearances[j], Settings.InclusionExponent) + Settings.InclusionPseudoCount;
        InclusionFromRaw(q0, eta);
    }

    /// <summary><c>make_weights</c>' passes from raw weights <paramref name="raw"/> in place of (aⱼ^α + ε): E4's history birth
    /// gives max(ĥⱼ, δ). The weights must be positive.</summary>
    public void InclusionProbabilities(ReadOnlySpan<double> raw, int q0, Span<double> eta)
    {
        if (eta.Length != raw.Length) throw new ArgumentException("One probability per variable.", nameof(eta));
        raw.CopyTo(eta);
        InclusionFromRaw(q0, eta);
    }

    // make_weights from the raw weights in `eta`, in place: v = raw/Σraw·q₀/p, then the passes.
    private void InclusionFromRaw(int q0, Span<double> eta)
    {
        int p = eta.Length;
        if (q0 < 1 || q0 > p) throw new ArgumentOutOfRangeException(nameof(q0));
        if (p == 1)
        {
            eta[0] = 1.0;
            return;
        }
        double total = 0.0;
        for (int j = 0; j < p; j++) total += eta[j];
        for (int j = 0; j < p; j++) eta[j] = eta[j] / total * q0 / p;
        double delta = 0.0, beta = 1.0, logP = Math.Log(p);
        Span<double> v = p <= StackLimit ? stackalloc double[p] : new double[p];
        eta.CopyTo(v);
        for (int pass = 0; pass < Settings.InclusionPasses; pass++)
        {
            double sum = 0.0;
            double target = Math.Log(q0 + delta) - logP;
            for (int j = 0; j < p; j++) sum += target / Math.Log(v[j]);
            beta = sum / p;
            double powered = 0.0;
            for (int j = 0; j < p; j++) powered += Math.Pow(v[j], beta);
            delta += q0 - powered;
        }
        for (int j = 0; j < p; j++) eta[j] = Math.Min(1.0, Math.Pow(v[j], beta));
    }

    /// <inheritdoc/>
    public BassTerm? DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts) => Draw(rng, counts, null);

    /// <summary>The birth under <paramref name="weights"/> (E4's history birth): η from <c>make_weights</c>' passes over the weights
    /// in place of (aⱼ^α + ε); the rest as <see cref="DrawBirth(Xoshiro256PlusPlus, BassTermCounts)"/> draws it.</summary>
    BassTerm? IBassWeightedBirth.DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts, BassVariableWeights weights)
    {
        ArgumentNullException.ThrowIfNull(weights);
        return Draw(rng, counts, weights);
    }

    /// <summary>ln q_b(term) under <paramref name="weights"/>; the counts do not enter, since the weights replace them.</summary>
    double IBassWeightedBirth.LogBirth(BassTerm term, BassTermCounts counts, int shift, BassVariableWeights weights)
    {
        ArgumentNullException.ThrowIfNull(weights);
        Check(counts);
        Prior.Validate(term);
        if (counts.Terms - shift < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
        int p = Space.Variables;
        Span<bool> included = p <= StackLimit ? stackalloc bool[p] : new bool[p];
        foreach (LegendreFactor f in term.Legendre) included[f.Variable] = true;
        int q = term.Legendre.Count, d = term.LegendreDegree;
        return LogInclusion(weights[BassPart.Legendre], included) + _logDegree[q][d] + LegendreTermSpace.LogComposition(d, q);
    }

    private BassTerm? Draw(Xoshiro256PlusPlus rng, BassTermCounts counts, BassVariableWeights? weights)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Check(counts);
        int p = Space.Variables, top = Space.MaxOrder;
        int q0 = DrawFrom(rng, _logOrder, 1, top);
        Span<double> eta = p <= StackLimit ? stackalloc double[p] : new double[p];
        if (weights is null) InclusionProbabilities(counts.AppearanceCounts(BassPart.Legendre), q0, eta);
        else InclusionProbabilities(weights[BassPart.Legendre], q0, eta);
        Span<bool> chi = p <= StackLimit ? stackalloc bool[p] : new bool[p];
        int q;
        do
        {
            q = 0;
            for (int j = 0; j < p; j++)
            {
                chi[j] = rng.NextDouble() < eta[j];
                if (chi[j]) q++;
            }
        } while (q == 0 || q > top);
        int degree = DrawFrom(rng, _logDegree[q], q, Space.MaxDegree);
        int[] parts = LegendreTermSpace.DrawComposition(rng, degree, q);
        var factors = new LegendreFactor[q];
        for (int j = 0, k = 0; j < p; j++)
            if (chi[j]) factors[k] = new LegendreFactor(j, parts[k++]);
        return new BassTerm(legendre: factors);
    }

    /// <inheritdoc/>
    public double LogBirth(BassTerm term, BassTermCounts counts) => LogBirth(term, counts, 0);

    /// <inheritdoc/>
    public double LogBirthRatio(BassTerm born, BassTermCounts counts) => -Math.Log(counts.Terms + 1.0) - LogBirth(born, counts, 0);

    /// <inheritdoc/>
    public double LogDeathRatio(BassTerm dying, BassTermCounts counts) => LogBirth(dying, counts, 1) + Math.Log(counts.Terms);

    /// <summary>ln P(χ): the inclusion's exact log probability, Σ_{q₀} J(q₀)·P(χ | q₀)/P(1 ≤ q ≤ q_max | q₀), under the appearance
    /// counts <paramref name="appearances"/>, for the variables <paramref name="included"/> mark.</summary>
    public double LogInclusion(ReadOnlySpan<int> appearances, ReadOnlySpan<bool> included)
    {
        int p = Space.Variables;
        if (appearances.Length != p || included.Length != p) throw new ArgumentException("One entry per variable.", nameof(included));
        Span<double> raw = p <= StackLimit ? stackalloc double[p] : new double[p];
        for (int j = 0; j < p; j++) raw[j] = Math.Pow(appearances[j], Settings.InclusionExponent) + Settings.InclusionPseudoCount;
        return LogInclusion(raw, included);
    }

    /// <summary>ln P(χ) as <see cref="LogInclusion(ReadOnlySpan{int}, ReadOnlySpan{bool})"/> reads it, from raw weights
    /// <paramref name="raw"/> in place of (aⱼ^α + ε).</summary>
    public double LogInclusion(ReadOnlySpan<double> raw, ReadOnlySpan<bool> included)
    {
        int p = Space.Variables, top = Space.MaxOrder;
        if (raw.Length != p || included.Length != p) throw new ArgumentException("One entry per variable.", nameof(included));
        Span<double> eta = p <= StackLimit ? stackalloc double[p] : new double[p];
        Span<double> terms = stackalloc double[top];
        for (int q0 = 1; q0 <= top; q0++)
        {
            InclusionProbabilities(raw, q0, eta);
            double log = 0.0;
            for (int j = 0; j < p; j++) log += included[j] ? Math.Log(eta[j]) : SpecialFunctions.Log1p(-eta[j]);
            terms[q0 - 1] = _logOrder[q0] + log - Math.Log(PoissonBinomial.Probability(eta, 1, top));
        }
        return LogSum(terms);
    }

    /// <inheritdoc/>
    public BassTermMove ProposeChange(Xoshiro256PlusPlus rng, BassTerm term, BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(rng);
        ArgumentNullException.ThrowIfNull(term);
        Check(counts);
        Prior.Validate(term);
        double mutation = DegreeMutation;
        bool degree = mutation >= 1.0 || (mutation > 0.0 && rng.NextDouble() < mutation);
        return degree ? ProposeDegree(rng, term) : ProposeSwap(rng, term, counts.AppearanceCounts(BassPart.Legendre));
    }

    /// <summary>The degree mutation: d′ from the birth's law at the term's order, a new uniform composition, the variables kept.</summary>
    public BassTermMove ProposeDegree(Xoshiro256PlusPlus rng, BassTerm term)
    {
        int q = term.Legendre.Count, d = term.LegendreDegree;
        int next = DrawFrom(rng, _logDegree[q], q, Space.MaxDegree);
        int[] parts = LegendreTermSpace.DrawComposition(rng, next, q);
        var factors = new LegendreFactor[q];
        for (int k = 0; k < q; k++) factors[k] = new LegendreFactor(term.Legendre[k].Variable, parts[k]);
        return new BassTermMove(BassMoveKind.DegreeChange, new BassTerm(legendre: factors), LogDegreeRatio(q, d, next));
    }

    /// <summary>The degree mutation's ln q(reverse)/q(forward) from degree <paramref name="from"/> to <paramref name="to"/> at order
    /// <paramref name="order"/>: each side's degree weight and composition probability.</summary>
    public double LogDegreeRatio(int order, int from, int to)
        => _logDegree[order][from] + LegendreTermSpace.LogComposition(from, order) - _logDegree[order][to] - LegendreTermSpace.LogComposition(to, order);

    /// <summary>The variable swap: a factor uniformly, a new variable among those the term does not use ∝ w + its appearances;
    /// a null move when the term uses every variable.</summary>
    public BassTermMove ProposeSwap(Xoshiro256PlusPlus rng, BassTerm term, ReadOnlySpan<int> appearances)
    {
        int p = Space.Variables, q = term.Legendre.Count;
        if (q == p) return new BassTermMove(BassMoveKind.VariableChange, null, double.NaN);
        int pick = BassDraws.Below(rng, q);
        Span<bool> held = p <= StackLimit ? stackalloc bool[p] : new bool[p];
        foreach (LegendreFactor f in term.Legendre) held[f.Variable] = true;
        double total = 0.0;
        for (int u = 0; u < p; u++)
            if (!held[u]) total += Settings.SwapPseudoCount + appearances[u];
        double target = rng.NextDouble() * total, cumulative = 0.0;
        int chosen = -1, last = -1;
        for (int u = 0; u < p && chosen < 0; u++)
        {
            if (held[u]) continue;
            last = u;
            cumulative += Settings.SwapPseudoCount + appearances[u];
            if (target < cumulative) chosen = u;
        }
        if (chosen < 0) chosen = last;
        BassTerm swapped = Swapped(term, pick, chosen);
        return new BassTermMove(BassMoveKind.VariableChange, swapped, LogSwapRatio(term, term.Legendre[pick].Variable, chosen, appearances));
    }

    /// <summary>The swap's ln q(reverse)/q(forward) for <paramref name="term"/> with variable <paramref name="from"/> swapped for
    /// <paramref name="to"/>, under the appearance counts <paramref name="appearances"/> before the swap: the factor's pick, 1/q, is
    /// the same both ways, and each side's variable draws among the variables its term leaves out.</summary>
    public double LogSwapRatio(BassTerm term, int from, int to, ReadOnlySpan<int> appearances)
    {
        int p = Space.Variables;
        Span<bool> held = p <= StackLimit ? stackalloc bool[p] : new bool[p];
        foreach (LegendreFactor f in term.Legendre) held[f.Variable] = true;
        if (!held[from] || held[to]) throw new ArgumentException("A swap takes a held variable for one the term leaves out.", nameof(to));
        double w = Settings.SwapPseudoCount;
        double forward = 0.0, reverse = 0.0;
        for (int u = 0; u < p; u++)
        {
            if (!held[u]) forward += w + appearances[u];
            // The swapped term leaves out the others the term left out, less `to`, plus `from`, under the counts after the swap.
            if ((!held[u] && u != to) || u == from) reverse += w + appearances[u] - (u == from ? 1 : 0);
        }
        double forwardPick = (w + appearances[to]) / forward;
        double reversePick = (w + appearances[from] - 1) / reverse;
        return Math.Log(reversePick) - Math.Log(forwardPick);
    }

    // ln q_b(term) under the counts less `shift` copies of the term itself: 0 for a birth into them, 1 for the reverse of a
    // death from them.
    private double LogBirth(BassTerm term, BassTermCounts counts, int shift)
    {
        Check(counts);
        Prior.Validate(term);
        if (counts.Terms - shift < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
        int p = Space.Variables;
        int[] held = counts.AppearanceCounts(BassPart.Legendre);
        Span<int> appearances = p <= StackLimit ? stackalloc int[p] : new int[p];
        held.CopyTo(appearances);
        Span<bool> included = p <= StackLimit ? stackalloc bool[p] : new bool[p];
        foreach (LegendreFactor f in term.Legendre)
        {
            included[f.Variable] = true;
            appearances[f.Variable] -= shift;
            if (appearances[f.Variable] < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
        }
        int q = term.Legendre.Count, d = term.LegendreDegree;
        return LogInclusion(appearances, included) + _logDegree[q][d] + LegendreTermSpace.LogComposition(d, q);
    }

    private static BassTerm Swapped(BassTerm term, int pick, int variable)
    {
        var factors = new LegendreFactor[term.Legendre.Count];
        for (int k = 0; k < factors.Length; k++)
            factors[k] = k == pick ? new LegendreFactor(variable, term.Legendre[k].Degree) : term.Legendre[k];
        return new BassTerm(legendre: factors);
    }

    // ln w(i) − ln Σ w over i in [lower, upper], in an array indexed by i, entries below `lower` −∞.
    private static double[] Normalized(int lower, int upper, Func<int, double> logWeight, int length)
    {
        var log = new double[length];
        Array.Fill(log, double.NegativeInfinity);
        for (int i = lower; i <= upper; i++) log[i] = logWeight(i);
        double total = LogSum(log.AsSpan(lower, upper - lower + 1));
        for (int i = lower; i <= upper; i++) log[i] -= total;
        return log;
    }

    // An index in [lower, upper] with probability exp(log[i]): one uniform against the cumulative weights.
    private static int DrawFrom(Xoshiro256PlusPlus rng, double[] log, int lower, int upper)
    {
        if (lower == upper) return lower;
        double target = rng.NextDouble(), cumulative = 0.0;
        for (int i = lower; i < upper; i++)
        {
            cumulative += Math.Exp(log[i]);
            if (target < cumulative) return i;
        }
        return upper;
    }

    private static double LogSum(ReadOnlySpan<double> values)
    {
        double top = double.NegativeInfinity;
        foreach (double v in values) top = Math.Max(top, v);
        if (double.IsNegativeInfinity(top)) return top;
        double sum = 0.0;
        foreach (double v in values) sum += Math.Exp(v - top);
        return top + Math.Log(sum);
    }

    private void Check(BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (!ReferenceEquals(counts.Prior, Prior)) throw new ArgumentException("The counts are over another prior's term space.", nameof(counts));
    }
}
