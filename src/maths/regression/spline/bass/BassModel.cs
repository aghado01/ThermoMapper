using System;
using Maths.Regression.Spline.Bars;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One fit's immutable parts, assembled once from the config and the data and shared by every chain: the resolution,
/// the input transform and the prepared inputs, the basis, the admission rules, the prior on structures, the
/// proposals and the collapsed target. Each chain adds its own mutable design (<see cref="BassStructureMoves"/>).
/// </summary>
public sealed class BassModel
{
    private BassModel(BassConfig config, BassData data, BassResolved resolved, BassInputTransform transform, BassInputs inputs,
                      BassBasis basis, BassAdmission admission, BassStructurePrior prior, BassTermProposal proposal, BassTarget target)
    {
        Config = config;
        Data = data;
        Resolved = resolved;
        Transform = transform;
        Inputs = inputs;
        Basis = basis;
        Admission = admission;
        Prior = prior;
        Proposal = proposal;
        Target = target;
    }

    /// <summary><paramref name="config"/> on <paramref name="data"/>, its target tempered in <paramref name="form"/>; null takes
    /// the config's (<see cref="TemperingSpec.Form"/> or <see cref="SmcSpec.Form"/>), the integrated form when it declares
    /// none.</summary>
    public static BassModel Create(BassConfig config, BassData data, TemperedForm? form = null)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(data);
        BassResolved resolved = BassResolved.Resolve(config, data);
        BassInputTransform transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        var basis = new BassBasis(inputs, config.Structure.Degree);
        BassAdmission admission = BassAdmission.From(config, resolved);
        BassStructurePrior prior = BassStructurePrior.From(config, resolved, inputs);
        BassTermProposal proposal = BassTermProposal.From(config, prior);
        var target = new BassTarget(config, resolved, prior, form ?? config.Sampling.Tempering?.Form ?? config.Sampling.Smc?.Form ?? TemperedForm.Integrated);
        return new(config, data, resolved, transform, inputs, basis, admission, prior, proposal, target);
    }

    /// <summary>The config.</summary>
    public BassConfig Config { get; }

    /// <summary>The data.</summary>
    public BassData Data { get; }

    /// <summary>What the config resolved to on the data.</summary>
    public BassResolved Resolved { get; }

    /// <summary>The input transform, fitted on the data.</summary>
    public BassInputTransform Transform { get; }

    /// <summary>The prepared training inputs.</summary>
    public BassInputs Inputs { get; }

    /// <summary>The basis.</summary>
    public BassBasis Basis { get; }

    /// <summary>D06's admission rules.</summary>
    public BassAdmission Admission { get; }

    /// <summary>The prior on structures (D01, D02).</summary>
    public BassStructurePrior Prior { get; }

    /// <summary>D01's proposals.</summary>
    public BassTermProposal Proposal { get; }

    /// <summary>The collapsed target (D03, D05).</summary>
    public BassTarget Target { get; }

    /// <summary>M_max, the basis cap.</summary>
    public int MaxBasis => Resolved.MaxBasis;

    /// <summary>A new intercept-only design over the data, one per chain.</summary>
    public BassDesign NewDesign() => new(Basis, Data, Admission);
}
