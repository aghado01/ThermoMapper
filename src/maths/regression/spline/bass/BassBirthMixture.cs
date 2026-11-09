using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>A birth's variable weights by part, as a history birth sets them: one positive weight per variable of each
/// part the fit has.</summary>
internal sealed class BassVariableWeights
{
    private readonly double[]?[] _parts;

    internal BassVariableWeights(double[]?[] parts)
    {
        ArgumentNullException.ThrowIfNull(parts);
        _parts = parts;
    }

    /// <summary>The weights over <paramref name="part"/>'s variables.</summary>
    public ReadOnlySpan<double> this[BassPart part]
        => _parts[(int)part] ?? throw new ArgumentException($"The weights hold no {part} part.", nameof(part));
}

/// <summary>A base kernel's birth under variable weights in place of the current terms' counts (E4's history birth).</summary>
internal interface IBassWeightedBirth
{
    /// <summary>A birth's term under <paramref name="weights"/>, or null for a draw that proposes nothing.</summary>
    BassTerm? DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts, BassVariableWeights weights);

    /// <summary>ln q_b(term) under <paramref name="weights"/>, under the counts less <paramref name="shift"/> copies of the term:
    /// 0 for a birth into them, 1 for the reverse of a death from them.</summary>
    double LogBirth(BassTerm term, BassTermCounts counts, int shift, BassVariableWeights weights);
}

/// <summary>
/// A chain's inclusion history for E4's history birth (Nott and Kohn 2005, <i>Biometrika</i> 92): ĥ_u, the share of the
/// recorded states in which some term uses variable u, per part, and the birth's weights max(ĥ_u, δ).
/// <list type="bullet">
/// <item><b>D10's policy.</b> States are recorded only while <see cref="Learning"/>, which a run sets for burn-in; the
/// weights are refreshed every <see cref="HistoryBirth.Refresh"/> records, from the past only, and <see cref="Freeze"/>
/// refreshes once more and stops, so production runs a fixed kernel.</item>
/// <item><b>Before the first refresh</b> there are no weights, and the history birth is the fit's own.</item>
/// <item><b>SMC</b> sets the weights from each stage's weighted population (<see cref="SetFromPopulation"/>).</item>
/// </list>
/// One history serves a chain, the rungs of its ladder, or an SMC run's particles; it is read, not written, while they move.
/// </summary>
public sealed class BassInclusionHistory
{
    private readonly BassStructurePrior _prior;
    private readonly long[]?[] _included = new long[]?[BassStructurePrior.PartCount];        // per part, the recorded states that use each variable
    private readonly bool[]?[] _scratch = new bool[]?[BassStructurePrior.PartCount];
    private double[]?[]? _frequencies;                            // ĥ per part at the last refresh

    /// <summary>A history over <paramref name="prior"/>'s parts under <paramref name="spec"/>, empty and not learning.</summary>
    public BassInclusionHistory(BassStructurePrior prior, HistoryBirth spec)
    {
        ArgumentNullException.ThrowIfNull(prior);
        ArgumentNullException.ThrowIfNull(spec);
        _prior = prior;
        Spec = spec;
        foreach (BassPart part in Enum.GetValues<BassPart>())
        {
            if (!prior.Has(part)) continue;
            _included[(int)part] = new long[prior.Variables(part)];
            _scratch[(int)part] = new bool[prior.Variables(part)];
        }
    }

    /// <summary>The settings.</summary>
    public HistoryBirth Spec { get; }

    /// <summary>Whether states are being recorded: burn-in.</summary>
    public bool Learning { get; private set; }

    /// <summary>Whether the weights are frozen: production.</summary>
    public bool Frozen { get; private set; }

    /// <summary>The states recorded.</summary>
    public long Recorded { get; private set; }

    /// <summary>The birth's weights, max(ĥ_u, δ) per part; null before the first refresh.</summary>
    internal BassVariableWeights? Weights { get; private set; }

    /// <summary>ĥ over <paramref name="part"/>'s variables at the last refresh, or null before it: Nott and Kohn's estimate of
    /// the marginal inclusion probabilities, a readout.</summary>
    public IReadOnlyList<double>? Frequencies(BassPart part) => _frequencies?[(int)part];

    /// <summary>The reading at the last refresh (<see cref="BassInclusionReading"/>); null before it.</summary>
    public BassInclusionReading? Reading()
    {
        if (_frequencies is not { } f) return null;
        var parts = new IReadOnlyList<double>?[BassStructurePrior.PartCount];
        for (int i = 0; i < parts.Length; i++) parts[i] = (double[]?)f[i]?.Clone();
        return new BassInclusionReading(Recorded, parts);
    }

    /// <summary>Starts recording: burn-in.</summary>
    public void Learn()
    {
        if (Frozen) throw new InvalidOperationException("The history is frozen.");
        Learning = true;
    }

    /// <summary>Records the structure <paramref name="terms"/> while learning, refreshing every
    /// <see cref="HistoryBirth.Refresh"/> records; outside learning, nothing.</summary>
    public void Record(IReadOnlyList<BassTerm> terms)
    {
        ArgumentNullException.ThrowIfNull(terms);
        if (!Learning) return;
        foreach (BassPart part in Enum.GetValues<BassPart>())
        {
            if (_included[(int)part] is not { } counts) continue;
            bool[] used = _scratch[(int)part]!;
            Array.Clear(used);
            foreach (BassTerm term in terms)
                foreach (int v in BassTermCounts.Variables(term, part)) used[v] = true;
            for (int v = 0; v < counts.Length; v++)
                if (used[v]) counts[v]++;
        }
        Recorded++;
        if (Recorded % Spec.Refresh == 0) Refresh();
    }

    /// <summary>Ends learning: one refresh from every record, when there is one, then the weights hold.</summary>
    public void Freeze()
    {
        if (Frozen) throw new InvalidOperationException("The history is frozen already.");
        if (Recorded > 0) Refresh();
        Learning = false;
        Frozen = true;
    }

    /// <summary>
    /// Sets ĥ from a population: Σᵢ wᵢ·[some term of structure i uses u] / Σᵢ wᵢ, over <paramref name="structures"/> and their
    /// <paramref name="weights"/>, not all zero. An SMC stage's kernel then reads these weights for every particle.
    /// </summary>
    public void SetFromPopulation(IReadOnlyList<IReadOnlyList<BassTerm>> structures, ReadOnlySpan<double> weights)
    {
        ArgumentNullException.ThrowIfNull(structures);
        if (structures.Count != weights.Length) throw new ArgumentException("One weight per structure.", nameof(weights));
        double total = 0.0;
        foreach (double w in weights) total += w;
        if (!(total > 0.0 && double.IsFinite(total))) throw new ArgumentException("The weights must sum to a positive finite value.", nameof(weights));
        var frequencies = new double[]?[BassStructurePrior.PartCount];
        foreach (BassPart part in Enum.GetValues<BassPart>())
        {
            if (_scratch[(int)part] is not { } used) continue;
            var h = new double[used.Length];
            for (int i = 0; i < structures.Count; i++)
            {
                if (weights[i] == 0.0) continue;
                Array.Clear(used);
                foreach (BassTerm term in structures[i])
                    foreach (int v in BassTermCounts.Variables(term, part)) used[v] = true;
                for (int v = 0; v < used.Length; v++)
                    if (used[v]) h[v] += weights[i];
            }
            for (int v = 0; v < h.Length; v++) h[v] /= total;
            frequencies[(int)part] = h;
        }
        Set(frequencies);
    }

    // ĥ = counts/recorded per part, and the weights from it.
    private void Refresh()
    {
        var frequencies = new double[]?[BassStructurePrior.PartCount];
        for (int index = 0; index < BassStructurePrior.PartCount; index++)
        {
            if (_included[index] is not { } counts) continue;
            var h = new double[counts.Length];
            for (int v = 0; v < h.Length; v++) h[v] = (double)counts[v] / Recorded;
            frequencies[index] = h;
        }
        Set(frequencies);
    }

    private void Set(double[]?[] frequencies)
    {
        var weights = new double[]?[BassStructurePrior.PartCount];
        for (int index = 0; index < BassStructurePrior.PartCount; index++)
        {
            if (frequencies[index] is not { } h) continue;
            var w = new double[h.Length];
            for (int v = 0; v < w.Length; v++) w[v] = Math.Max(h[v], Spec.Floor);
            weights[index] = w;
        }
        _frequencies = frequencies;
        Weights = new BassVariableWeights(weights);
    }

    /// <summary>The restart state: the counts, the records, ĥ, and the phase.</summary>
    internal BassInclusionHistorySnapshot Snapshot()
    {
        var included = new long[]?[BassStructurePrior.PartCount];
        for (int i = 0; i < BassStructurePrior.PartCount; i++) included[i] = (long[]?)_included[i]?.Clone();
        double[]?[]? frequencies = null;
        if (_frequencies is { } f)
        {
            frequencies = new double[]?[BassStructurePrior.PartCount];
            for (int i = 0; i < BassStructurePrior.PartCount; i++) frequencies[i] = (double[]?)f[i]?.Clone();
        }
        return new BassInclusionHistorySnapshot(included, Recorded, frequencies, Learning, Frozen);
    }

    /// <summary>A history over <paramref name="prior"/> under <paramref name="spec"/> as <paramref name="snapshot"/> left it.</summary>
    internal static BassInclusionHistory Restore(BassStructurePrior prior, HistoryBirth spec, BassInclusionHistorySnapshot snapshot)
    {
        ArgumentNullException.ThrowIfNull(snapshot);
        var history = new BassInclusionHistory(prior, spec);
        for (int i = 0; i < BassStructurePrior.PartCount; i++)
        {
            if ((history._included[i] is null) != (snapshot.Included[i] is null) ||
                (history._included[i] is { } own && own.Length != snapshot.Included[i]!.Length))
                throw new ArgumentException("The snapshot is another term space's history.", nameof(snapshot));
            if (snapshot.Included[i] is { } counts) Array.Copy(counts, history._included[i]!, counts.Length);
        }
        history.Recorded = snapshot.Recorded;
        if (snapshot.Frequencies is { } f)
        {
            var copy = new double[]?[BassStructurePrior.PartCount];
            for (int i = 0; i < BassStructurePrior.PartCount; i++) copy[i] = (double[]?)f[i]?.Clone();
            history.Set(copy);
        }
        history.Learning = snapshot.Learning;
        history.Frozen = snapshot.Frozen;
        return history;
    }

    /// <summary>A copy, for a moves object that restores its own and is then given the shared one.</summary>
    internal BassInclusionHistory Copy() => Restore(_prior, Spec, Snapshot());
}

/// <summary>A <see cref="BassInclusionHistory"/>'s restart state.</summary>
internal sealed record BassInclusionHistorySnapshot(long[]?[] Included, long Recorded, double[]?[]? Frequencies, bool Learning, bool Frozen);

/// <summary>
/// E4's birth mixture over a fit's own kernel (the proposal-kernels brief): q(t | S) = (1 − ω_H)·q_base(t | S) +
/// ω_H·q_H(t | S), q_H the base birth under the history's weights (<see cref="BassInclusionHistory"/>), and q_H = q_base
/// before its first refresh.
/// <list type="bullet">
/// <item><b>Exactness.</b> ln q is the log-sum-exp of the components' exact ln q, and a death's reverse is the same sum under
/// the counts without the dying term; the base share keeps every structure the base birth reaches.</item>
/// <item><b>Draws.</b> One uniform picks the component, then the component draws.</item>
/// <item><b>Deaths and changes</b> are the base kernel's.</item>
/// </list>
/// One per chain's moves, reading the history its run shares with it.
/// </summary>
internal sealed class BassBirthMixture : IBassProposalKernel
{
    private readonly IBassProposalKernel _base;
    private readonly IBassWeightedBirth _weighted;
    private readonly double _share, _logShare, _logBase;

    /// <summary>The mixture over <paramref name="kernel"/> with the history birth <paramref name="spec"/>, reading
    /// <paramref name="history"/>.</summary>
    public BassBirthMixture(IBassProposalKernel kernel, HistoryBirth spec, BassInclusionHistory history)
    {
        ArgumentNullException.ThrowIfNull(kernel);
        ArgumentNullException.ThrowIfNull(spec);
        ArgumentNullException.ThrowIfNull(history);
        _base = kernel;
        _weighted = kernel as IBassWeightedBirth ?? throw new ArgumentException("The kernel has no weighted birth.", nameof(kernel));
        _share = spec.Share;
        _logShare = Math.Log(spec.Share);
        _logBase = SpecialFunctions.Log1p(-spec.Share);
        History = history;
    }

    /// <summary>The history whose weights the history birth reads.</summary>
    public BassInclusionHistory History { get; set; }

    /// <inheritdoc/>
    public BassStructurePrior Prior => _base.Prior;

    /// <inheritdoc/>
    public BassTerm? DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(rng);
        bool history = rng.NextDouble() < _share;
        return history && History.Weights is { } weights ? _weighted.DrawBirth(rng, counts, weights) : _base.DrawBirth(rng, counts);
    }

    /// <inheritdoc/>
    public double LogBirth(BassTerm term, BassTermCounts counts) => Mix(_base.LogBirth(term, counts), term, counts, 0);

    /// <inheritdoc/>
    public double LogBirthRatio(BassTerm born, BassTermCounts counts) => -Math.Log(counts.Terms + 1.0) - LogBirth(born, counts);

    /// <inheritdoc/>
    public double LogDeathRatio(BassTerm dying, BassTermCounts counts)
    {
        double logM = Math.Log(counts.Terms);
        return Mix(_base.LogDeathRatio(dying, counts) - logM, dying, counts, 1) + logM;
    }

    /// <inheritdoc/>
    public BassTermMove ProposeChange(Xoshiro256PlusPlus rng, BassTerm term, BassTermCounts counts) => _base.ProposeChange(rng, term, counts);

    // ln[(1 − ω)·q_base + ω·q_H] from ln q_base, or ln q_base itself before the history's first refresh.
    private double Mix(double logBase, BassTerm term, BassTermCounts counts, int shift)
    {
        if (History.Weights is not { } weights) return logBase;
        double a = _logBase + logBase, b = _logShare + _weighted.LogBirth(term, counts, shift, weights);
        double top = Math.Max(a, b);
        return double.IsNegativeInfinity(top) ? top : top + Math.Log(Math.Exp(a - top) + Math.Exp(b - top));
    }
}
