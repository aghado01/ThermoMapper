using System;
using Maths.Regression.GPrior;
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
                      BassBasis basis, BassAdmission admission, BassStructurePrior prior, BassTermProposal? proposal, IBassProposalKernel kernel,
                      BassTarget target, IBassObservationModel observation, WorkingResponse working)
    {
        Config = config;
        Observation = observation;
        Working = working;
        Data = data;
        Resolved = resolved;
        Transform = transform;
        Inputs = inputs;
        Basis = basis;
        Admission = admission;
        Prior = prior;
        _proposal = proposal;
        Kernel = kernel;
        Target = target;
    }

    /// <summary><paramref name="config"/> on <paramref name="data"/>, its target tempered in <paramref name="form"/>; null takes
    /// the config's (<see cref="TemperingSpec.Form"/> or <see cref="SmcSpec.Form"/>), the integrated form when it declares
    /// none. The noise is the config's: <see cref="VarianceMeanMixture"/> for a declared <see cref="BassConfig.Noise"/>, Gaussian
    /// otherwise.</summary>
    public static BassModel Create(BassConfig config, BassData data, TemperedForm? form = null)
    {
        ArgumentNullException.ThrowIfNull(config);
        return Create(config, data, config.Noise is { } noise ? new VarianceMeanMixture(noise, config.Sigma2) : GaussianObservation.Instance, form);
    }

    /// <summary><see cref="Create(BassConfig, BassData, TemperedForm?)"/> with <paramref name="observation"/> as the noise
    /// the Gaussian core fits (F-X1), for a config that declares none; Gaussian BASS is the overload above.</summary>
    public static BassModel Create(BassConfig config, BassData data, IBassObservationModel observation, TemperedForm? form = null)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(observation);
        if (config.Noise is { } declared && !(observation is VarianceMeanMixture mixture && mixture.Noise == declared))
            throw new ArgumentException("The config declares its noise, so the observation model is the config's.", nameof(observation));
        WorkingResponse working = observation.Working(data.Response);
        BassResolved resolved = BassResolved.Resolve(config, data, working);
        BassInputTransform transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        var basis = new BassBasis(inputs, config.Structure.Degree);
        BassAdmission admission = BassAdmission.From(config, resolved);
        BassStructurePrior prior = BassStructurePrior.From(config, resolved, inputs);
        BassTermProposal? proposal = config.Structure.Legendre is null ? BassTermProposal.From(config, prior) : null;
        IBassProposalKernel kernel = proposal ?? (IBassProposalKernel)new LegendreTermProposal(prior, config.Structure.Legendre!);
        var target = new BassTarget(config, resolved, prior, form ?? config.Sampling.Tempering?.Form ?? config.Sampling.Smc?.Form ?? TemperedForm.Integrated, working);
        return new(config, data, resolved, transform, inputs, basis, admission, prior, proposal, kernel, target, observation, working);
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

    private readonly BassTermProposal? _proposal;

    /// <summary>D01's proposals over the hinges; refused for a polynomial-chaos fit, which draws from KHAOS's (<see cref="Kernel"/>).</summary>
    public BassTermProposal Proposal
        => _proposal ?? throw new InvalidOperationException("A Legendre basis draws from KHAOS's proposal (Kernel), not D01's.");

    /// <summary>The proposal kernel the structure moves draw from (F-X4): D01's over hinges, KHAOS's over a Legendre term space.</summary>
    public IBassProposalKernel Kernel { get; }

    /// <summary>The collapsed target (D03, D05).</summary>
    public BassTarget Target { get; }

    /// <summary>The observation model: what the Gaussian core fits of the response (F-X1).</summary>
    public IBassObservationModel Observation { get; }

    /// <summary>What the observation model hands the Gaussian core: the working response and weights every design starts on.</summary>
    public WorkingResponse Working { get; }

    /// <summary>
    /// Whether the priors read the working response as it stands: no <see cref="BassPriorReference"/>, with level centring or a σ² prior
    /// that reads s². A working response that moves is not data, so such a prior would read parameters, with conditionals that share no
    /// joint, and a chain of this model refuses to move (<see cref="BassStructureMoves.ProposeFollow"/>). The routes that read nothing
    /// are a flat or zero intercept with a σ² prior that reads no s² (Jeffreys, a declared inverse gamma, a held σ²).
    /// </summary>
    public bool PriorsReadResponse => Working.Reference is null && (Config.Intercept == GPriorCentre.Mean || Config.Sigma2.ReadsResponse);

    /// <summary>M_max, the basis cap.</summary>
    public int MaxBasis => Resolved.MaxBasis;

    /// <summary>A new intercept-only design over the data, one per chain, on <paramref name="working"/>; null is the model's.</summary>
    public BassDesign NewDesign(WorkingResponse? working = null)
        => new(Basis, Data, Admission, Config.Coefficients.Form is RidgeForm, working ?? Working, Config.Intercept, Config.Coefficients.Form is RawRidgeForm);
}
