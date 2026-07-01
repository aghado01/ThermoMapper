using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Profiling.Signals;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Execution.Sinks;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

using Graphs.Models.Potts;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Configuration for <see cref="FixedGridSweepStrategy"/>: a user-supplied
/// temperature grid plus per-point sampler parameters.
/// </summary>
/// <remarks>
/// <para>The user owns the grid — no log-spacing assumption, no
/// auto-bracketing. Order is preserved end-to-end; the assembled
/// <see cref="SweepProfile"/> simply records what was measured at each
/// temperature, in whatever order was given. Pair with a dense linspace
/// (Tier-0 χ(T) survey), a hand-picked diagnostic grid, or the output of
/// any pre-existing analysis.</para>
///
/// <para><b>Why this exists.</b> The adaptive scheduler's coarse→dense
/// pipeline is opinionated about where to spend cycles — useful when you
/// want one good T, hostile when you want to <i>see</i> the curves
/// before deciding what "good" means. A well-sampled fixed-grid sweep
/// is the natural input for cross-temperature analyses (pseudo-transition
/// detection on the way to T_c) and surrogate-model fitting.</para>
/// </remarks>
public sealed class FixedGridSweepConfig
{
    /// <summary>
    /// User-supplied temperature grid. Order preserved; no sorting,
    /// dedup, or log-spacing applied. Must be non-empty.
    /// </summary>
    public required IReadOnlyList<double> Temperatures { get; init; }

    /// <summary>
    /// Independent replicas per temperature. Each replica is a fresh
    /// PottsModel with a deterministic per-task seed (round=0,
    /// replica-indexed). Their per-T measurements feed into
    /// <see cref="SweepProfile"/> as additional data points at the
    /// same T — the profile aggregator averages them.
    /// </summary>
    public int Replicas { get; init; } = 1;

    /// <summary>Burn-in and measurement cycles per sweep probe.</summary>
    public RunBudget SweepBudget { get; init; } = new(200, 1000);

    /// <summary>
    /// Sampler-level Potts configuration. See
    /// <see cref="PottsModelConfig"/>.
    /// </summary>
    public PottsModelConfig Sampler { get; init; } = new();

    /// <summary>Burn-in and measurement cycles for the chosen-T equilibrium pass.</summary>
    public RunBudget EquilibriumBudget { get; init; } = new(1000, 5000);

    /// <summary>
    /// What sufficient-statistics the sweep collects at every T.
    /// Defaults to <see cref="AccumulationSpec.None"/> (cheap: scalar
    /// moments only). Set to <see cref="AccumulationSpec.Currencies"/> for
    /// a rich sweep that avoids a second equilibrium pass at chosen-T and
    /// makes bond-observables-over-T available.
    /// </summary>
    public AccumulationSpec Accumulation { get; init; } = AccumulationSpec.None;

    /// <summary>
    /// Which susceptibility channel drives chosen-T peak-finding. All three
    /// (FK cluster / FK reduced / magnetization) are always assembled into the
    /// <see cref="SweepProfile"/> as channels for free comparison; this selects the
    /// primary. Defaults to FK cluster (SW-native, lower variance).
    /// </summary>
    public SusceptibilityKind SusceptibilityKind { get; init; } = SusceptibilityKind.FkCluster;

    /// <summary>
    /// Root seed for reproducible per-task RNG derivation. Null draws
    /// from OS entropy (non-reproducible).
    /// </summary>
    public int? BaseSeed { get; init; }

    /// <summary>Worker-budget policy for the executor's flat-task pool (defaults to the auto policy).</summary>
    public WorkerBudgetPolicy Parallelism { get; init; } = new();

    /// <summary>
    /// Optional directory for per-task checkpoint persistence (and
    /// resume on re-run). When non-null, each (T, replica) probe and
    /// the chosen-T equilibrium pass are written as <c>.spcx</c> files
    /// via <see cref="Execution.Sinks.SpcxDiskFrameSink"/>; re-running the
    /// strategy against the same directory skips tasks whose checkpoints
    /// already exist. When null, the sweep runs in-memory only (no
    /// persistence, no resume).
    /// </summary>
    public string? CheckpointDirectory { get; init; }
}

/// <summary>
/// Brute-force sweep strategy: runs every point on a user-supplied
/// temperature grid (with optional replicas), assembles a
/// <see cref="SweepProfile"/>, picks chosen-T via the supplied analyzer,
/// and mints the <see cref="Affinities"/> (and optionally
/// <see cref="Alignments"/>) currencies at that T.
/// </summary>
/// <remarks>
/// <para><b>Use case.</b> The diagnostic / research workhorse, and the
/// current bedrock sweep: "just give me χ(T) on this exact grid." Use when
/// you want to <i>see</i> the susceptibility (and label-entropy,
/// specific-heat, bond-frequency) curves across the whole T-domain, or when
/// downstream analysis (pseudo-transition detection, surrogate fitting)
/// needs a well-sampled trajectory.</para>
///
/// <para><b>Parallelization.</b> Flattens (T × replica) into a single
/// task list dispatched against the resolved <see cref="FixedGridSweepConfig.Parallelism"/>
/// budget through the shared executor pool.</para>
///
/// <para><b>Two-stage efficiency.</b> When
/// <see cref="FixedGridSweepConfig.Accumulation"/> is
/// <see cref="AccumulationSpec.None"/> (the default), a second
/// equilibrium pass at chosen-T mints the currencies — one extra run at the
/// end. When <c>Accumulation.Affinities</c> is already true (rich sweep),
/// the chosen-T run from the sweep itself is reused and no second pass is
/// needed.</para>
/// </remarks>
public sealed class FixedGridSweepStrategy : ISweepStrategy
{
    private readonly FixedGridSweepConfig _config;

    public FixedGridSweepStrategy(FixedGridSweepConfig config)
    {
        ArgumentNullException.ThrowIfNull(config);
        if (config.Temperatures.Count == 0)
            throw new ArgumentException(
                "FixedGridSweepConfig.Temperatures must be non-empty.", nameof(config));
        _config = config;
    }

    /// <inheritdoc />
    public SweepResult Run(CsrGraph graph, ISignalAnalyzer? analyzer = null)
    {
        var cfg = _config;
        var signalAnalyzer = analyzer ?? new ChiPeakSignalAnalyzer();
        var sw = Stopwatch.StartNew();

        if (graph.NodeCount < 2)
            return SweepKernel.BuildTrivialResult(graph, cfg.Sampler.Q, sw.Elapsed);

        var executor = new SpcExecutor();
        var execOptions = new SpcExecutionOptions
        {
            FrameSink   = ResolveSink(cfg.CheckpointDirectory),
            Parallelism = cfg.Parallelism,
        };

        // ── Stage 1: the (T × replica) sweep ─────────────────────────────
        int tCount = cfg.Temperatures.Count;
        var sweepTasks = new List<SpcTaskSpec>(tCount * cfg.Replicas);
        for (int rep = 0; rep < cfg.Replicas; rep++)
        {
            for (int tIdx = 0; tIdx < tCount; tIdx++)
            {
                double T = cfg.Temperatures[tIdx];
                sweepTasks.Add(new SpcTaskSpec
                {
                    Temperature    = T,
                    ReplicaIndex   = rep,
                    Budget         = cfg.SweepBudget,
                    Q              = cfg.Sampler.Q,
                    Accumulation   = cfg.Accumulation,
                    BaseSeed       = SpcSeedHelper.Derive(cfg.BaseSeed, T, replica: rep, round: 0),
                    CheckpointPath = CheckpointFor(cfg.CheckpointDirectory, T, rep, schedule: 0),
                });
            }
        }

        var stageTimings = new Dictionary<string, TimeSpan>(StringComparer.Ordinal);
        var stageSw = Stopwatch.StartNew();
        SpcBatchResult sweepBatch = executor.RunBatch(graph, sweepTasks, default, execOptions);
        var sweepRuns = HydrateAlignedRuns(sweepBatch, sweepTasks, "sweep");
        stageTimings["sweep"] = stageSw.Elapsed;

        var profile = SweepProfile.From(sweepRuns, cfg.SusceptibilityKind);
        SpPlateauResult plateau = SpcProfileAnalysis.SpPlateau(profile);
        double chosenT  = plateau.TClus;
        double stability = SpcProfileAnalysis.ComputeStability(profile);

        // ── Stage 2: mint chosen-T currencies ────────────────────────────
        // When the sweep already collected Affinities at every T, reuse the
        // chosen-T run (replica-0). Otherwise run one equilibrium pass with
        // Currencies to mint both per-edge precursors.
        stageSw.Restart();
        Accumulator eqAcc;

        if (cfg.Accumulation.Affinities)
        {
            // Rich sweep — find replica-0 run at chosenT.
            SpcRunResult? richRun = sweepRuns
                .FirstOrDefault(r => Math.Abs(r.Accumulator.Temperature - chosenT) < 1e-9
                                  && r.Accumulator.ReplicaIndex == 0);

            if (richRun is null)
            {
                // Shouldn't happen on a well-formed probe batch, but fall through
                // gracefully to a re-run rather than throwing.
                eqAcc = RunEquilibriumPass(graph, cfg, chosenT, executor, execOptions);
            }
            else
            {
                eqAcc = richRun.Accumulator;
            }
        }
        else
        {
            eqAcc = RunEquilibriumPass(graph, cfg, chosenT, executor, execOptions);
        }

        stageTimings["equilibrium"] = stageSw.Elapsed;

        Affinities   chosenAffinities   = SwCurrencies.ToAffinities(eqAcc);
        Alignments?  chosenAlignments   = eqAcc.SpinAgreementCount is null  ? null : SwCurrencies.ToAlignments(eqAcc);
        CoMembership? chosenCoMembership = eqAcc.CoMembershipCount is null  ? null : SwCurrencies.ToCoMembership(eqAcc);

        int totalCycles = sweepTasks.Count * (cfg.SweepBudget.BurnIn + cfg.SweepBudget.Cycles)
                        + cfg.EquilibriumBudget.BurnIn + cfg.EquilibriumBudget.Cycles;

        ProfileCriteria criteria = signalAnalyzer.Analyze(profile);

        return new SweepResult
        {
            Summary = new SweepSummary
            {
                SubgraphNodes     = graph.NodeCount,
                SubgraphEdges     = graph.Targets.Length / 2,
                ChosenTemperature = chosenT,
                StabilityScore    = stability,
                TotalCyclesUsed   = totalCycles,
                EarlyStopped      = false,
                Elapsed           = sw.Elapsed,
                StageTimings      = stageTimings,
                Profile           = profile,
            },
            SweepRuns        = sweepRuns,
            ProfileCriteria  = criteria,
            Graph            = graph,
            ChosenAffinities    = chosenAffinities,
            ChosenAlignments    = chosenAlignments,
            ChosenCoMembership  = chosenCoMembership,
        };
    }

    private Accumulator RunEquilibriumPass(
        CsrGraph             graph,
        FixedGridSweepConfig cfg,
        double               chosenT,
        SpcExecutor          executor,
        SpcExecutionOptions  execOptions)
    {
        var eqTask = new SpcTaskSpec
        {
            Temperature    = chosenT,
            ReplicaIndex   = 0,
            Budget         = cfg.EquilibriumBudget,
            Q              = cfg.Sampler.Q,
            Accumulation   = AccumulationSpec.Currencies,
            BaseSeed       = SpcSeedHelper.Derive(cfg.BaseSeed, chosenT, replica: 0, round: -1),
            CheckpointPath = CheckpointFor(cfg.CheckpointDirectory, chosenT, replica: 0, schedule: 1),
        };
        SpcBatchResult eqBatch = executor.RunBatch(graph, new[] { eqTask }, default, execOptions);
        SpcRunResult eqRun = eqBatch.AlignedRuns[0]
            ?? throw new InvalidOperationException(
                "Chosen-T equilibrium task produced no result. This indicates a sink " +
                "configuration issue — disk sink + a missing checkpoint path, or a hydration failure.");
        return eqRun.Accumulator;
    }

    /// <summary>
    /// Pulls task-aligned results out of <see cref="SpcBatchResult.AlignedRuns"/>
    /// and throws if any slot is null — covers both fresh-run and resume paths.
    /// </summary>
    private static List<SpcRunResult> HydrateAlignedRuns(
        SpcBatchResult batch,
        IReadOnlyList<SpcTaskSpec> tasks,
        string stageLabel)
    {
        var list = new List<SpcRunResult>(tasks.Count);
        for (int i = 0; i < tasks.Count; i++)
        {
            SpcRunResult? r = batch.AlignedRuns[i];
            if (r is null)
            {
                var task = tasks[i];
                throw new InvalidOperationException(
                    $"FixedGridSweepStrategy {stageLabel} task {i} " +
                    $"(T={task.Temperature.ToString("G6", CultureInfo.InvariantCulture)}, " +
                    $"replica={task.ReplicaIndex}) produced no result and no cache hit.");
            }
            list.Add(r);
        }
        return list;
    }

    private static IFrameSink ResolveSink(string? checkpointDirectory)
        => checkpointDirectory is null
            ? NullFrameSink.Instance
            : SpcxDiskFrameSink.Instance;

    // Flat, strategy-agnostic checkpoint bag: one file per work item, identity
    // (T, replica, schedule) carried by the name.
    private static string? CheckpointFor(string? root, double T, int replica, int schedule)
    {
        if (root is null) return null;
        string fileName = $"T_{T.ToString("F5", CultureInfo.InvariantCulture)}_rep_{replica}_s{schedule}.spcx";
        return Path.Combine(root, fileName);
    }
}
