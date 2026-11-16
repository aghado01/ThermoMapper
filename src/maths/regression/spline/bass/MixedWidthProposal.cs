using System;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>Accumulates proposals over the disjoint mixed-width kinds. Births choose available kinds uniformly and adapt within
/// the chosen kind; deaths choose among all terms. Width flips share the change slot with the ordinary directional change.</summary>
public sealed class MixedWidthProposal : IBassProposalKernel
{
    private readonly DirectionalTermProposal _directional;

    /// <summary>The mixed proposal under the configured pseudo-counts and directional concentration.</summary>
    public MixedWidthProposal(BassConfig config, BassStructurePrior prior)
    {
        if (!prior.MixedWidths) throw new ArgumentException("The prior needs mixed widths.", nameof(prior));
        Prior = prior;
        _directional = DirectionalTermProposal.From(config, prior);
    }

    /// <inheritdoc/>
    public BassStructurePrior Prior { get; }

    private BassStructurePrior.Part? Part(BassPart part) => Prior[part];

    private int DrawOrder(Xoshiro256PlusPlus rng, BassStructurePrior.Part part, BassTermCounts counts, int minimum)
    {
        double w = _directional.OrderWeight;
        double u = rng.NextDouble() * (w * (part.Cap - minimum + 1) + counts.Terms);
        for (int k = minimum; k < part.Cap; k++)
        {
            u -= w + counts.Order(part.Kind, k);
            if (u < 0.0) return k;
        }
        return part.Cap;
    }

    /// <inheritdoc/>
    public BassTerm? DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts)
    {
        Check(counts);
        int kind = BassDraws.Below(rng, Prior.AvailableKinds);
        BassTermCounts own = counts.ForKind(kind);
        ProjectionDirection? direction = null;
        if (kind != 2)
        {
            var part = Part(BassPart.Directional)!;
            int a = DrawOrder(rng, part, own, 1);
            int[] variables = BassTermProposal.DrawVariablesOf(rng, _directional.VariableWeight, part, own.AppearanceCounts(part.Kind), a);
            direction = Prior.Directional!.DrawDirection(rng, variables);
        }
        var gate = Part(BassPart.Categorical);
        int k = gate is null ? 0 : DrawOrder(rng, gate, own, kind == 2 ? 1 : 0);
        var cells = new CategoricalFactor[k];
        if (k > 0)
        {
            int[] variables = BassTermProposal.DrawVariablesOf(rng, _directional.VariableWeight, gate!, own.AppearanceCounts(BassPart.Categorical), k);
            for (int f = 0; f < k; f++) cells[f] = gate!.DrawSubset(rng, variables[f], Prior.Subsets);
        }
        return Prior.MixedTermWithKnot(rng, kind, direction, cells);
    }

    /// <inheritdoc/>
    public double LogBirth(BassTerm term, BassTermCounts counts)
    {
        Check(counts);
        Prior.Validate(term);
        int kind = BassStructurePrior.KindOf(term);
        BassTermCounts own = counts.ForKind(kind);
        double log = -Math.Log(Prior.AvailableKinds);
        foreach (BassPart p in new[] { BassPart.Directional, BassPart.Categorical })
        {
            var part = Part(p);
            if (part is null || (p == BassPart.Directional && kind == 2)) continue;
            int minimum = p == BassPart.Directional || kind == 2 ? 1 : 0;
            int k = BassStructurePrior.Order(term, p);
            double w = _directional.OrderWeight;
            log += Math.Log(w + own.Order(p, k)) - Math.Log(w * (part.Cap - minimum + 1) + own.Terms);
            if (k == 0) continue;
            log += BassTermProposal.LogVariablesOf(_directional.VariableWeight, term, p, part, own.AppearanceCounts(p), 0);
            log += part.LogFactors(term);
        }
        return log;
    }

    /// <inheritdoc/>
    public double LogBirthRatio(BassTerm born, BassTermCounts counts) => -Math.Log(counts.Terms + 1.0) - LogBirth(born, counts);

    /// <inheritdoc/>
    public double LogDeathRatio(BassTerm dying, BassTermCounts counts)
    {
        var reduced = counts.Clone();
        reduced.Remove(dying);
        return LogBirth(dying, reduced) + Math.Log(counts.Terms);
    }

    /// <inheritdoc/>
    public BassTermMove ProposeChange(Xoshiro256PlusPlus rng, BassTerm term, BassTermCounts counts)
    {
        Check(counts);
        Prior.Validate(term);
        if (term.Directional is { } h && rng.NextDouble() < 0.5)
            return new(BassMoveKind.WidthChange, new BassTerm(categorical: term.Categorical,
                directional: new DirectionalHinge(h.Direction, h.Knot, h.Columns == 1 ? 4 : 1)), 0.0);
        BassTermMove move = _directional.ProposeChange(rng, term, counts);
        if (move.Term?.Directional is { } moved && term.Directional is { } old)
            return move with { Term = new BassTerm(categorical: move.Term.Categorical,
                directional: new DirectionalHinge(moved.Direction, moved.Knot, old.Columns)) };
        return move;
    }

    private void Check(BassTermCounts counts)
    {
        if (!ReferenceEquals(counts.Prior, Prior)) throw new ArgumentException("The counts belong to another prior.", nameof(counts));
    }
}
