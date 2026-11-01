using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// One grid point read by M1's method (<see cref="BassMixingCostTests"/>), generalized to the cost brief's grid: the
/// chain burned by <see cref="BassCost.Burned"/>, then each part of a sweep at its state, round-robin with the sweep's
/// batches, read in µs per call or in bytes per call (<see cref="BassCost.Meter"/>).
/// <list type="bullet">
///   <item><b>The parts</b> are M1's: a birth's draw, ratio, column, cross products, bordered Cholesky and design, the
///   whole birth move; a death's and a change's proposal; the target; the three commits, the factor's rebuild amortized
///   over its cadence; τ's block with its evaluations; one composition; and the sweep, a step and a composition. A direct
///   functional point's column and cross products are its case and output factors and <see cref="KhatriRao.CentredBorder"/>
///   with the two output-weighted cross products, as the design takes them. A categorical point adds a subset's
///   membership test over the cases; every point adds an empty part, the timer's own overhead.</item>
///   <item><b>Its frequencies</b> per sweep come from the sweep's batches, each on a copy of the chain restored at the
///   state with its own streams, as M1 takes them.</item>
///   <item><b>Its analytic counts</b> per call, in multiply-adds (or element passes, for a column), follow the code's
///   loops: the column n(ō + 4) (the fill, ō kernels, the support count, the mean, the centring); its cross products
///   n(M + 3) (M border dots, zᵀz, zᵀỹ, zᵀr), direct 3M(n_x + m) + 2(N + n_x) (every current factor re-centred, two
///   passes each, its entry, and the two weighted products); the bordered proposal 1.5M² (Border's forward and back
///   solves, the quadratic form's forward solve); a commit NM (Refit's explicit residual); the rebuild M³/3 (the strict
///   factor and diag(A⁻¹)); a death 2M²; a change a column, n(M + 2) and 3.5M²; a composition M²/2 (the
///   innovation's back solve).</item>
/// </list>
/// It reads; it asserts only that the timed proposals committed nothing (<see cref="BassMixing.Untouched"/>), the chain
/// stayed where it was, and every reading is finite where it is defined.
/// </summary>
internal static class BassCostSweep
{
    /// <summary>A point's reading: each part's nine batch values, the frequencies and the analytic counts.</summary>
    internal sealed class Reading(BassCost.Point point)
    {
        public BassCost.Point Point { get; } = point;
        public int BurnSweeps, M, Cap, Rebuild;
        public double HeldOrder, PoolOrder, EmptyShare, BurnSeconds, ReadSeconds;
        public readonly Dictionary<string, double[]> Parts = [];
        public readonly Dictionary<string, double> Units = [];
        public double BirthsPerSweep, NonEmptyPerSweep, DeathsPerSweep, ChangesPerSweep, AdmittedPerSweep, RefusedPerSweep,
                      BirthAccepts, DeathAccepts, ChangeAccepts, TauEvaluations, SweepEndM, ThreadBytesPerSweep, MembershipCases = double.NaN;
        public long RebuildsFired, CommitsTimed;

        public double Of(string name) => Parts.TryGetValue(name, out double[]? v) ? BassCost.Median(v) : double.NaN;
        public double IqrOf(string name) => Parts.TryGetValue(name, out double[]? v) ? BassCost.Iqr(v) : double.NaN;

        // At the cap the design refuses a birth before its column and the chain proposes none, so the birth's residuals
        // are not read there.
        public bool AtCap => M == Cap;
        public double Candidate => AtCap ? 0.0 : (Of("birth move") - Of("draw")) / (1.0 - EmptyShare) - Of("ratio") - Of("design");
        public double DesignOther => AtCap ? 0.0 : Of("design") - Of("column") - Of("cross") - Of("cholesky");
        public double Commit(string kind) => Of($"commit {kind}") - (CommitsTimed > 0 ? (double)RebuildsFired / CommitsTimed : 0.0) * Of("rebuild")
                                             + (Rebuild > 0 ? Of("rebuild") / Rebuild : 0.0);
        public double RebuildsPerSweep => Rebuild > 0 ? (BirthAccepts + DeathAccepts + ChangeAccepts) / Rebuild : 0.0;
        public double ColumnsPerSweep => NonEmptyPerSweep + ChangesPerSweep;

        /// <summary>The sweep's parts per sweep, each at its frequency, then the rest: µs (or bytes) per sweep.</summary>
        public (string Name, double Value)[] PerSweep()
        {
            (string, double)[] parts =
            [
                ("birth: draw", BirthsPerSweep * Of("draw")),
                ("birth: ratio", NonEmptyPerSweep * Of("ratio")),
                ("birth: column", NonEmptyPerSweep * Of("column")),
                ("birth: cross products", NonEmptyPerSweep * Of("cross")),
                ("birth: Cholesky", NonEmptyPerSweep * Of("cholesky")),
                ("birth: design other", NonEmptyPerSweep * DesignOther),
                ("birth: candidate", NonEmptyPerSweep * Candidate),
                ("death proposal", Weigh(DeathsPerSweep, Of("death"))),
                ("change proposal", Weigh(ChangesPerSweep, Of("change"))),
                ("target", AdmittedPerSweep * Of("target")),
                ("commit", Weigh(BirthAccepts, Commit("birth")) + Weigh(DeathAccepts, Commit("death")) + Weigh(ChangeAccepts, Commit("change"))),
                ("τ block", Of("tau")),
                ("composition", Of("compose")),
            ];
            return [.. parts, ("rest", Of("sweep") - parts.Sum(x => x.Item2))];
        }

        // A part the chain never takes here costs nothing; one the state cannot read (a death or change at M = 0, which the
        // sweeps leave) is counted as nothing too, and named in Unread.
        private static double Weigh(double frequency, double value) => frequency == 0.0 || double.IsNaN(value) ? 0.0 : frequency * value;

        /// <summary>The parts the sweeps take that the state could not read: a death and a change at M = 0.</summary>
        public string[] Unread => [.. new[] { ("death", DeathsPerSweep), ("change", ChangesPerSweep) }
                                         .Where(q => q.Item2 > 0 && double.IsNaN(Of(q.Item1))).Select(q => q.Item1)];
    }

    /// <summary>Reads <paramref name="x"/> under <paramref name="meter"/>.</summary>
    internal static Reading Read(BassCost.Point x, BassCost.Meter meter)
    {
        var clock = Stopwatch.StartNew();
        (BassModel model, BassChain chain, int burn) = BassCost.Burned(x);
        BassState state = chain.Current;
        int m = state.Count;
        var reading = new Reading(x) { BurnSweeps = burn, M = m, Cap = model.MaxBasis, BurnSeconds = clock.Elapsed.TotalSeconds };
        clock.Restart();
        reading.HeldOrder = m == 0 ? 0.0 : state.Terms.Average(BassCost.OrderOf);
        bool functional = model.Data.Response is FunctionalResponse;
        int n = model.Data.Response.Cases, locations = model.Data.Response.Locations, rows = n * locations;

        var untouched = new BassMixing.Untouched(model, chain);
        BassDesignSnapshot snapshot = untouched.Before;
        BassStructureMoves structure = chain.Structure;
        BassDesign design = structure.DesignAt(state);
        BassTermCounts counts = BassMixing.Counts(model, state);
        IReadOnlyList<IRjMove<BassState>> moves = structure.Moves;
        IRjMove<BassState> birth = BassMixing.Move(moves, "birth"), death = BassMixing.Move(moves, "death"), change = BassMixing.Move(moves, "change");
        var rng = new Xoshiro256PlusPlus(seed: 113);

        // The pool: births as the chain would draw them here, each one's columns, border, cross products and b.
        var terms = new List<BassTerm>();
        int draws = 0;
        while (terms.Count < BassCost.Pool)
        {
            draws++;
            if (model.Proposal.DrawBirth(rng, counts) is { } t) terms.Add(t);
        }
        reading.EmptyShare = 1.0 - (double)BassCost.Pool / draws;
        reading.PoolOrder = terms.Average(BassCost.OrderOf);
        double[][] current = snapshot.Columns!.Cases;
        double[][]? outputs = snapshot.Columns.Outputs;
        double[] residual = design.Residual.ToArray();
        double[] centred = [.. model.Data.Response.Values.ToArray().Select(v => v - design.ResponseMean)];
        // As the design: a repeated term is refused before its column, an unsupported one after it.
        (double[] A, double[]? H)?[] columns = [.. terms.Select(t => state.Terms.Contains(t) ? null : Column(model, t))];
        var borders = new double[BassCost.Pool][];
        var diagonals = new double[BassCost.Pool];
        var bs = new double[BassCost.Pool][];
        for (int k = 0; k < BassCost.Pool; k++)
        {
            if (columns[k] is not var (a, h)) continue;
            if (h is null)
            {
                borders[k] = [.. current.Select(c => PairwiseSum.Dot(a, c))];
                diagonals[k] = PairwiseSum.Dot(a, a);
            }
            else if (m > 0) borders[k] = KhatriRao.CentredBorder(current, outputs!, a, h, out diagonals[k]);
            else columns[k] = null;   // a direct design with no terms: its self product only, left unread
            if (columns[k] is not null) bs[k] = [.. snapshot.Columns.CrossResponses, Cross(a, h, centred, locations)];
        }

        int i = 0;
        double sink = 0.0;
        // The design's factor, or at its capacity the factor its growth would make, a strict factorization of the same Gram.
        var factor = (BorderedCholesky)design.Factor;
        if (m >= factor.Capacity)
        {
            var gram = new double[m, m];
            for (int r = 0; r < m; r++)
                for (int c = 0; c < m; c++) gram[r, c] = factor.Entry(r, c);
            factor = new BorderedCholesky(m + 1, factor.RebuildInterval);
            factor.Decompose(gram);
        }
        reading.Rebuild = snapshot.Factor.RebuildInterval;
        // The target's evaluation of admitted candidates, each a fresh evaluation as the chain makes one.
        var candidates = new List<BassState>();
        for (int k = 0; k < 20 * BassCost.Pool && candidates.Count < BassCost.Pool; k++)
            foreach (IRjMove<BassState> move in moves)
                if (move.Propose(state, rng) is { } q && model.Target.Refusal(q.Candidate.Fit) is null) candidates.Add(q.Candidate);

        // The commits, forced on proposals from the state, the design restored before each and after every batch. A commit
        // whose accept fell on the factor's rebuild cadence (a rebuild leaves its drift recorded, a strict accept clears
        // it) has the rebuild taken out afterwards; every commit carries the rebuild's amortized share instead.
        BassTerm[] supported = [.. terms.Where((t, k) => columns[k] is not null)];
        long[] fired = new long[2];
        Func<int, double> Commit(Func<int, DesignProposal> propose)
        {
            int k = 0;
            Func<int, double> each = BassCost.EachBatch(meter, read =>
            {
                design.Restore(snapshot);
                if (!propose(k++).Admitted) return -1;
                long t0 = read();
                design.Accept();
                long spent = read() - t0;
                var held = (BorderedCholesky)design.Factor;
                if (held.UpdatesSinceRebuild == 0 && held.LastDrift is not null) fired[0]++;
                fired[1]++;
                return spent;
            });
            return round =>
            {
                double mean = each(round);
                design.Restore(snapshot);
                return mean;
            };
        }

        // τ's block on a copy of the chain restored at the state, its target's evaluations counted.
        BassChainSnapshot fixedChain = chain.Snapshot();
        var counting = new Counting(model.Target);
        var tau = new BassTauBlock(counting);
        BassChain moved = BassChain.Restore(model, fixedChain);
        long updates = 0;
        var table = new BassTermTable();

        // The sweep, a step and a composition: each batch on a copy of the chain restored at the state with its own
        // streams, so every batch starts there; the moves' frequencies over the read batches.
        int[] seeds = SeedTree.Derive(127, 2 * (BassCost.Repeats + 1));
        int sweepBatch = BassCost.Scaled(BassCost.SweepBatch);
        long sweeps = 0, births = 0, empties = 0, deathsTaken = 0, changes = 0, refused = 0, bornAccepted = 0, diedAccepted = 0,
             changedAccepted = 0, endM = 0, threadBytes = 0;
        double Sweeps(int round)
        {
            (ulong a0, ulong a1, ulong a2, ulong a3) = new Xoshiro256PlusPlus(seeds[2 * (round + 1)]).SaveState();
            BassChainSnapshot reseeded = fixedChain with
            {
                Chain = fixedChain.Chain with { Rng0 = a0, Rng1 = a1, Rng2 = a2, Rng3 = a3 },
                BlockRng = new Xoshiro256PlusPlus(seeds[2 * (round + 1) + 1]).SaveState(),
            };
            BassChain run = BassChain.Restore(model, reseeded);
            IReadOnlyList<MoveStat> before = run.Moves.MoveStats();
            long emptyBefore = run.Structure.EmptyBirths, refusedBefore = run.Structure.Refusals.Total;
            (long p0, long th0) = BassCost.Allocated();
            long t0 = Stopwatch.GetTimestamp();
            for (int k = 0; k < sweepBatch; k++)
            {
                run.Step();
                run.Compose(rng, table);
            }
            long ticks = Stopwatch.GetTimestamp() - t0;
            (long p1, long th1) = BassCost.Allocated();
            double value = meter == BassCost.Meter.Time ? ticks * 1e6 / Stopwatch.Frequency / sweepBatch : (double)(p1 - p0) / sweepBatch;
            if (round < 0) return value;   // the warm-up round
            IReadOnlyList<MoveStat> after = run.Moves.MoveStats();
            sweeps += sweepBatch;
            births += Delta(after, before, "birth", false);
            deathsTaken += Delta(after, before, "death", false);
            changes += Delta(after, before, "change", false);
            bornAccepted += Delta(after, before, "birth", true);
            diedAccepted += Delta(after, before, "death", true);
            changedAccepted += Delta(after, before, "change", true);
            empties += run.Structure.EmptyBirths - emptyBefore;
            refused += run.Structure.Refusals.Total - refusedBefore;
            endM += run.Current.Count;
            threadBytes += th1 - th0;
            return value;
        }

        var timers = new List<(string Name, Func<int, double> Batch)>
        {
            ("empty", BassCost.Batch(meter, static () => { })),
            ("draw", BassCost.Batch(meter, () => model.Proposal.DrawBirth(rng, counts))),
            ("ratio", BassCost.Batch(meter, () => model.Proposal.LogBirthRatio(terms[i++ % BassCost.Pool], counts))),
            ("column", BassCost.Batch(meter, () =>
            {
                BassTerm t = terms[i++ % BassCost.Pool];
                if (!state.Terms.Contains(t)) Column(model, t);
            })),
            ("cross", BassCost.Batch(meter, () =>
            {
                if (columns[i++ % BassCost.Pool] is var (a, h)) sink += CrossProducts(a, h, current, outputs, centred, residual, locations);
            })),
            ("cholesky", BassCost.Batch(meter, () =>
            {
                int k = i++ % BassCost.Pool;
                if (columns[k] is not null) sink += BorderedProposal(factor, borders[k], diagonals[k], bs[k]);
            })),
            ("design", BassCost.Batch(meter, () => design.ProposeBirth(terms[i++ % BassCost.Pool]))),
            ("birth move", BassCost.Batch(meter, () => birth.Propose(state, rng))),
            ("target", BassCost.Batch(meter, () => sink += Temperature.Posterior.LogDensity(model.Target.Evaluate(candidates[i++ % candidates.Count])))),
            ("compose", BassCost.Batch(meter, () => structure.Compose(state, rng, table))),
            ("tau", BassCost.Batch(meter, () => { tau.Update(moved.Moves, rng); updates++; })),
            ("sweep", Sweeps),
        };
        if (m > 0)
        {
            timers.Add(("death", BassCost.Batch(meter, () => death.Propose(state, rng))));
            timers.Add(("change", BassCost.Batch(meter, () => change.Propose(state, rng))));
            timers.Add(("rebuild", BassCost.EachBatch(meter, read =>
            {
                BorderedCholesky copy = snapshot.Factor.Copy();
                long t0 = read();
                copy.Rebuild();
                return read() - t0;
            })));
            timers.Add(("commit death", Commit(k => design.ProposeDeath(k % m))));
            if (supported.Length > 0) timers.Add(("commit change", Commit(k => design.ProposeChange(k % m, supported[k % supported.Length]))));
        }
        // At the cap a birth commits from one term below it: a death committed first, unread.
        if (supported.Length > 0)
            timers.Add(("commit birth", Commit(k =>
            {
                if (m == model.MaxBasis)
                {
                    design.ProposeDeath(k % m);
                    design.Accept();
                }
                return design.ProposeBirth(supported[k % supported.Length]);
            })));
        // A categorical point: one subset of half the first categorical's levels, its membership test over the cases.
        if (x.Surface == BassCost.Surface.Categorical)
        {
            var subset = new CategoricalFactor(0, Enumerable.Range(0, x.L).Where(l => l % 2 == 0));
            int[] codes = [.. Enumerable.Range(0, n).Select(c => c % x.L)];   // the generator's balanced codes
            reading.MembershipCases = n;
            timers.Add(("membership", BassCost.Batch(meter, () => sink += Membership(subset, codes), batch: 20)));
        }

        foreach ((string name, double[] values) in BassCost.Interleave(timers)) reading.Parts[name] = values;
        untouched.Check();
        reading.RebuildsFired = fired[0];
        reading.CommitsTimed = fired[1];
        reading.TauEvaluations = updates > 0 ? (double)counting.Evaluations / updates : double.NaN;
        reading.SweepEndM = (double)endM / BassCost.Repeats;
        reading.NonEmptyPerSweep = (double)births / sweeps;
        reading.BirthsPerSweep = (double)(births + empties) / sweeps;
        reading.DeathsPerSweep = (double)deathsTaken / sweeps;
        reading.ChangesPerSweep = (double)changes / sweeps;
        reading.AdmittedPerSweep = (double)(births + deathsTaken + changes - refused) / sweeps;
        reading.RefusedPerSweep = (double)refused / sweeps;
        reading.BirthAccepts = (double)bornAccepted / sweeps;
        reading.DeathAccepts = (double)diedAccepted / sweeps;
        reading.ChangeAccepts = (double)changedAccepted / sweeps;
        reading.ThreadBytesPerSweep = (double)threadBytes / sweeps;
        reading.ReadSeconds = clock.Elapsed.TotalSeconds;

        // The analytic counts per call, at the state's M and the pool's mean order.
        double o = reading.PoolOrder, mm = m;
        double column = functional
            ? n * (terms.Average(t => t.Continuous.Count + t.Categorical.Count) + 2) + locations * (terms.Average(t => t.Functional.Count) + 2)
            : n * (o + 4);
        reading.Units["column"] = column;
        reading.Units["cross"] = functional ? 3 * mm * (n + locations) + 2.0 * (rows + n) : n * (mm + 3);
        reading.Units["cholesky"] = 1.5 * mm * mm;
        reading.Units["commit"] = (double)rows * mm;
        reading.Units["rebuild"] = mm * mm * mm / 3;
        reading.Units["death"] = 2 * mm * mm;
        reading.Units["change"] = column + (functional ? reading.Units["cross"] : n * (mm + 2)) + 3.5 * mm * mm;
        reading.Units["compose"] = mm * mm / 2;
        reading.Units["membership"] = n;

        Assert.Same(state, chain.Current);
        Assert.True(double.IsFinite(sink));
        foreach ((string name, double value) in reading.PerSweep()) Assert.True(double.IsFinite(value), $"{x.Key}: {name} is {value}");
        return reading;
    }

    // A move's attempts or acceptances between two readings of the chain's statistics.
    private static long Delta(IReadOnlyList<MoveStat> after, IReadOnlyList<MoveStat> before, string key, bool accepted)
    {
        MoveStat a = after.Single(s => s.Key == key), b = before.Single(s => s.Key == key);
        return accepted ? a.Accepted - b.Accepted : a.Attempts - b.Attempts;
    }

    // The new term's columns as the design evaluates them (BassDesign's private Evaluate, step for step): scalar, the case
    // column, its support, then centred; direct, the case-side and output-side factors as evaluated, each supported. Null
    // when a part falls short of the support rule. The test-side replicas are compiled fully optimized, as the engine's
    // hot paths are by the time they are read.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal static (double[] A, double[]? H)? Column(BassModel model, BassTerm term)
    {
        BassBasis basis = model.Basis;
        basis.Validate(term);
        var a = new double[basis.Inputs.Cases.Cases];
        basis.CaseColumn(term, basis.Inputs.Cases, a);
        BassAdmission rules = model.Admission;
        if (!Supported(rules, a, rules.SupportMinimum)) return null;
        if (basis.Inputs.Outputs is not { } at)
        {
            double mean = PairwiseSum.Of(a) / a.Length;
            for (int k = 0; k < a.Length; k++) a[k] -= mean;
            return (a, null);
        }
        var h = new double[at.Locations];
        basis.OutputColumn(term, at, h);
        if (!Supported(rules, h, rules.SupportMinimumFunctional ?? 0.0)) return null;
        return (a, h);
    }

    private static bool Supported(BassAdmission rules, double[] part, double minimum)
    {
        if (rules.Support == SupportKind.ColumnMass) return PairwiseSum.Of(part) >= minimum;
        int nonzero = 0;
        foreach (double v in part)
            if (v != 0.0) nonzero++;
        return nonzero >= minimum;
    }

    // zᵀv for the new column and a zero-mean v over all N rows: scalar z·v; direct aᵀVh under independent noise.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static double Cross(double[] a, double[]? h, double[] v, int locations)
    {
        if (h is null) return PairwiseSum.Dot(a, v);
        var vh = new double[a.Length];
        for (int i = 0; i < a.Length; i++) vh[i] = PairwiseSum.Dot(v.AsSpan(i * locations, locations), h);
        return PairwiseSum.Dot(a, vh);
    }

    // The design's cross products for a birth's columns: its border against every current column and its self product
    // (scalar dots; direct, the centred border), then its products with ỹ and the residual.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static double CrossProducts(double[] a, double[]? h, double[][] current, double[][]? outputs, double[] centred, double[] residual, int locations)
    {
        double sum = 0.0;
        if (h is null)
        {
            foreach (double[] c in current) sum += PairwiseSum.Dot(a, c);
            sum += PairwiseSum.Dot(a, a);
        }
        else
        {
            double[] border = KhatriRao.CentredBorder(current, outputs!, a, h, out double diagonal);
            sum += diagonal + (border.Length > 0 ? border[0] : 0.0);
        }
        return sum + Cross(a, h, centred, locations) + Cross(a, h, residual, locations);
    }

    // The factor's bordered proposal and, from a positive definite one, its quadratic form, as the design takes them.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static double BorderedProposal(BorderedCholesky factor, double[] border, double diagonal, double[] b)
        => factor.ProposeBirth(border, diagonal).PositiveDefinite ? factor.ProposedInverseQuadraticForm(b) : 0.0;

    // A subset's membership test over every case's level, as a categorical kernel takes it per case.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static double Membership(CategoricalFactor subset, int[] codes)
    {
        int inside = 0;
        foreach (int c in codes)
            if (subset.Contains(c)) inside++;
        return inside;
    }

    // τ's block's target, its evaluations counted.
    private sealed class Counting(BassTarget target) : ISamplingTarget<BassState>
    {
        public long Evaluations { get; private set; }

        public ITargetEvaluation Evaluate(BassState state)
        {
            Evaluations++;
            return target.Evaluate(state);
        }
    }
}
