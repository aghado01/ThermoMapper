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
                      BassTarget target, IBassObservationModel observation, WorkingResponse working, WorkingResponse baseReading, RandomInterceptBlock? block)
    {
        RandomIntercept = block;
        Config = config;
        Observation = observation;
        Working = working;
        BaseReading = baseReading;
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
        if (data.Response is MultiOutputResponse && (observation is not GaussianObservation || config.RandomIntercept is not null))
            throw new ArgumentException("Several outputs take Gaussian noise without random intercepts in E34's slice 1; the noise by row and the random "
                                        + "intercepts (decisions 6 and 7) arrive in slice 4.", nameof(config));
        if (config.Sampling.Births?.Residual is not null && data.Response is not ScalarResponse)
            throw new ArgumentException("E4's residual birth (Sampling.Births.Residual) reads a scalar response's residual; a direct functional "
                                        + "response's, or several outputs', is a matrix, which no score here reads.", nameof(config));
        if (config.Sampling.Births?.Residual is not null && !System.Linq.Enumerable.Any(data.Predictors, p => p is ContinuousPredictor))
            throw new ArgumentException("E4's residual birth (Sampling.Births.Residual) adds continuous factors, and the data hold no continuous "
                                        + "predictor.", nameof(config));
        WorkingResponse working = observation.Working(data.Response);
        WorkingResponse baseReading = working;
        RandomInterceptBlock? block = null;
        if (config.RandomIntercept is { } intercepts)
        {
            config.Validate();
            CategoricalPredictor grouping = data.Grouping
                ?? throw new ArgumentException("The config declares random intercepts (RandomIntercept), and the data hold no grouping.", nameof(data));
            block = RandomInterceptBlock.Of(grouping, intercepts);
            // A drawn s moves the whitened reading, so the priors read the response's own moments, frozen, as a fit without the block
            // reads them: a prior that read the reading at s would read a parameter (PriorsReadResponse).
            if (intercepts.ScalePrior is not null && working.Reference is null && working.IsResponseItself && working.Latents is null)
                baseReading = new WorkingResponse(working.Response, output: working.Output, meaning: working.Meaning, priorGram: working.PriorGram,
                    reference: config.Intercept == GPriorCentre.Zero
                    ? BassPriorReference.About(working.Values, 0.0)
                    : BassPriorReference.Moments(working.Values));
            working = block.Whitened(baseReading);
        }
        BassResolved resolved = BassResolved.Resolve(config, data, working);
        BassInputTransform transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        var basis = new BassBasis(inputs, config.Structure.Degree);
        BassAdmission admission = BassAdmission.From(config, resolved);
        BassStructurePrior prior = BassStructurePrior.From(config, resolved, inputs);
        BassTermProposal? proposal = config.Structure.Legendre is null ? BassTermProposal.From(config, prior) : null;
        IBassProposalKernel kernel = proposal ?? (IBassProposalKernel)new LegendreTermProposal(prior, config.Structure.Legendre!);
        var target = new BassTarget(config, resolved, prior, form ?? config.Sampling.Tempering?.Form ?? config.Sampling.Smc?.Form ?? TemperedForm.Integrated, working);
        return new(config, data, resolved, transform, inputs, basis, admission, prior, proposal, kernel, target, observation, working, baseReading, block);
    }

    /// <summary>The random-intercept block (E21) the config declares, integrated out of every marginal by whitening, which the working
    /// response carries (<see cref="WorkingResponse.Whitening"/>); null for a config without one. Under a drawn s it is the block at the
    /// chain's start, and a state's reading carries the block at its own s.</summary>
    public RandomInterceptBlock? RandomIntercept { get; }

    /// <summary>Whether the random-intercept block's scale is drawn (<see cref="BassRandomIntercept.ScalePrior"/>).</summary>
    public bool ScaleDrawn => RandomIntercept?.Spec.ScalePrior is not null;

    // Whether a state's reading moves with it: a noise's latents, or a drawn scale's whitening.
    internal bool ReadingMoves => Observation is VarianceMeanMixture || ScaleDrawn;

    // The reading before the block's whitening, with the priors' frozen reference: what each scale's whitening takes.
    internal WorkingResponse BaseReading { get; }

    // The block at a draw's scale: the model's at a held s, or at the s the draw carries.
    internal RandomInterceptBlock? BlockAt(BassDraw draw) => RandomIntercept is { } block && draw.GroupScale is { } s ? block.At(s) : RandomIntercept;

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
    public bool PriorsReadResponse => Working.Reference is null && (Config.Intercept == GPriorCentre.Mean ||
        (Data.Response is MultiOutputResponse ? Config.OutputCovariance?.ReadsResponse ?? true : Config.Sigma2.ReadsResponse));

    /// <summary>M_max, the basis cap.</summary>
    public int MaxBasis => Resolved.MaxBasis;

    /// <summary>A new intercept-only design over the data, one per chain, on <paramref name="working"/>; null is the model's.</summary>
    public BassDesign NewDesign(WorkingResponse? working = null)
        => new(Basis, Data, Admission, Config.Coefficients.Form is RidgeForm, working ?? Working, Config.Intercept, BassPriorPrecision.Of(Config.Coefficients.Form));
}
