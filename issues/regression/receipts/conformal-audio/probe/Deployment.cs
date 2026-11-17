using Maths.Distributions;
using Maths.Regression;
using Maths.Regression.Conformal;
using Maths.Rng;

namespace ConformalAudio;

internal static class Transport
{
    internal static void Run()
    {
        var config = App.Settings.GetProperty("Transport"); int n = App.I(config, "Groups"); double alpha = App.D(config, "Alpha");
        double[] factorial = new double[n + 1]; for (int i = 1; i <= n; i++) factorial[i] = factorial[i - 1] + Math.Log(i);
        double probabilitySum = 0, weightedMiss = 0, unweightedMiss = 0, weightedSize = 0, unweightedSize = 0, essMean = 0;
        double literalMiss = 0, literalDifferenceMass = 0; int differences = 0, states = 0;
        using var file = new StreamWriter(Path.Combine(App.Output, "transport-counts.jsonl"));
        CandidateLogDensity candidates = new(new double[,] { { -1, 0, -1, 0 } }, new[] { new[] { -1.0, 1.0 }, new[] { -1.0, 1.0 } }, []) { Forms = ConformalForms.Split };
        for (int positiveHigh = 0; positiveHigh <= n; positiveHigh++)
            for (int negativeHigh = 0; negativeHigh <= n - positiveHigh; negativeHigh++)
            {
                int low = n - positiveHigh - negativeHigh, high = positiveHigh + negativeHigh;
                double probability = Math.Exp(factorial[n] - factorial[low] - factorial[positiveHigh] - factorial[negativeHigh] + low * Math.Log(0.95) + high * Math.Log(0.025));
                double[] scores = Enumerable.Repeat(0.0, low + positiveHigh).Concat(Enumerable.Repeat(-1.0, negativeHigh)).ToArray();
                double[] integers = Enumerable.Repeat(10.0, low).Concat(Enumerable.Repeat(190.0, high)).ToArray();
                double[] literal = Enumerable.Repeat(10.0 / 19, low).Concat(Enumerable.Repeat(10.0, high)).ToArray();
                var cal = App.Seam(scores);
                var count = SplitConformal.Density(cal, candidates, alpha);
                var exact = SplitConformal.Density(cal, candidates, alpha, new ShiftWeights(integers, [10, 10, 190, 190]));
                var raw = SplitConformal.Density(cal, candidates, alpha, new ShiftWeights(literal, [10.0 / 19, 10.0 / 19, 10, 10]));
                bool[] wa = exact.Sets.Select(s => s.AcceptsPlain(0)).ToArray(), ua = count.Sets.Select(s => s.AcceptsPlain(0)).ToArray(), la = raw.Sets.Select(s => s.AcceptsPlain(0)).ToArray();
                for (int h = 0; h < 2; h++)
                {
                    int candidate = h == 0 ? 10 : 190;
                    bool expected = 10L * (190L * negativeHigh + candidate) > 10L * low + 190L * high + candidate;
                    App.Require(wa[h] == expected, "Integer-weight enumeration.");
                    App.Require(ua[h] == (negativeHigh + 1 > Math.Floor(alpha * (n + 1))), "Unweighted enumeration.");
                    App.Require(exact.Sets[h].AcceptsPlain(1) && count.Sets[h].AcceptsPlain(1) && raw.Sets[h].AcceptsPlain(1), "Positive sign always accepted.");
                }
                if (!wa.SequenceEqual(la)) { differences++; literalDifferenceMass += probability; }
                probabilitySum += probability; states++;
                weightedMiss += probability * (wa[1] ? 0 : 0.25); unweightedMiss += probability * (ua[1] ? 0 : 0.25); literalMiss += probability * (la[1] ? 0 : 0.25);
                weightedSize += probability * (1 + 0.5 * ((wa[0] ? 1 : 0) + (wa[1] ? 1 : 0)));
                unweightedSize += probability * (1 + 0.5 * ((ua[0] ? 1 : 0) + (ua[1] ? 1 : 0)));
                double ess = Math.Pow(10.0 * low + 190.0 * high, 2) / (100.0 * low + 36100.0 * high); essMean += probability * ess;
                App.Line(file, new { Low = low, PositiveHigh = positiveHigh, NegativeHigh = negativeHigh, Probability = probability,
                    WeightedNegativeAccepted = wa, UnweightedNegativeAccepted = ua, LiteralNegativeAccepted = la, Ess = ess });
            }
        App.Require(Math.Abs(probabilitySum - 1) < 1e-12 && weightedMiss <= alpha + 1e-12, "Transport probability/coverage.");
        var failed = SplitConformal.Density(App.Seam(new double[n]), candidates, alpha, new ShiftWeights(Enumerable.Repeat(1.0, n).ToArray(), [1, 1, 1, 1]));
        App.Require(!failed.Sets[0].AcceptsPlain(0) && failed.Sets[0].AcceptsPlain(1), "Conditional-law failure.");
        App.Write("transport-summary.json", new { Passed = true, States = states, ProbabilitySum = probabilitySum,
            WeightedMarginalMiss = weightedMiss, WeightedHighStratumMiss = 2 * weightedMiss, WeightedLowStratumMiss = 0,
            UnweightedMarginalMiss = unweightedMiss, UnweightedHighStratumMiss = 2 * unweightedMiss, UnweightedLowStratumMiss = 0,
            WeightedCardinality = weightedSize, WeightedFullAlphabetFrequency = weightedSize - 1,
            UnweightedCardinality = unweightedSize, UnweightedFullAlphabetFrequency = unweightedSize - 1,
            ExpectedCalibrationEss = essMean, LiteralMarginalMiss = literalMiss, LiteralBoundaryDifferenceStates = differences,
            LiteralBoundaryDifferenceProbability = literalDifferenceMass, FailedConditionalLawMiss = 1, FailedConditionalLawCardinality = 1 });
    }
}

internal static class Risk
{
    private sealed record Conversation(double[] Prefix, double[] Budget, double[] Completion);
    private static Conversation[] Generate(int seed, int n)
    {
        return SeedTree.Derive(seed, n).Select(s => {
            int[] child = SeedTree.Derive(s, 2); var rng = new Xoshiro256PlusPlus(child[0]);
            int length = rng.NextDouble() < 0.1 ? 0 : 1 + rng.NextInt(8); double remaining = 0.6 * rng.NextDouble();
            var turn = new Xoshiro256PlusPlus(child[1]); var prefix = new double[length]; var budget = new double[length]; var completion = new double[length];
            for (int j = 0; j < length; j++) { prefix[j] = 0.5 + turn.NextDouble(); budget[j] = 0.45 + 0.1 * turn.NextDouble(); completion[j] = prefix[j] + 0.1 + remaining + 0.05 * turn.NextDouble(); }
            return new Conversation(prefix, budget, completion);
        }).ToArray();
    }
    private static (double Loss, double StartRate, double TotalDelay, int Starts) Evaluate(Conversation c, double? delay)
    {
        if (c.Prefix.Length == 0) return (0, 0, 0, 0);
        int early = 0, starts = 0; double totalDelay = 0;
        for (int t = 0; t < c.Prefix.Length; t++)
        {
            if (delay is null || 0.1 + delay > c.Budget[t]) continue;
            double start = c.Prefix[t] + 0.1 + delay.Value;
            starts++; totalDelay += 0.1 + delay.Value; if (start < c.Completion[t]) early++;
        }
        return (early / (double)c.Prefix.Length, starts / (double)c.Prefix.Length, totalDelay, starts);
    }
    internal static void Run(int seed)
    {
        var config = App.Settings.GetProperty("Risk"); double alpha = App.D(config, "Alpha");
        double?[] delays = App.Doubles(config, "Delays").Select(v => (double?)v).Append(null).ToArray();
        var equality = RiskControl.Mean(new double[,] { { 1, 0.5, 0 }, { 1, 0, 0 } }, 1, 0.5);
        var fallback = RiskControl.Mean(new double[1, 3], 1, 0.1);
        App.Require(equality.Index == 1 && equality.Met && equality.Reading == 0.5 && fallback.Index == 2 && !fallback.Met, "Risk equality/fallback controls.");
        var witness = new Conversation([0, 1.9], [1, 0.1], [1, 2]);
        var before = Evaluate(witness, 0); var after = Evaluate(witness, 0.01);
        App.Require(before.Loss == 0.5 && after.Loss == 0.5 && before.Loss / before.StartRate == 0.5 && after.Loss / after.StartRate == 1,
            "Eligible-turn versus among-start denominator.");
        App.Require(Evaluate(witness, null).Loss == 0 && Evaluate(new Conversation([], [], []), 0).Loss == 0, "Abstain/empty convention.");
        int[] seeds = SeedTree.Derive(seed, App.I(config, "Splits")); var risks = new List<double>(); var starts = new List<double>(); var selected = new List<int>();
        using var output = new StreamWriter(Path.Combine(App.Output, "risk-splits.jsonl"));
        using var raw = new StreamWriter(Path.Combine(App.Output, "risk-conversations.jsonl"));
        double sumDelays = 0; int countStarts = 0, emptyCount = 0;
        for (int rep = 0; rep < seeds.Length; rep++)
        {
            int[] child = SeedTree.Derive(seeds[rep], 2); Conversation[] calibration = Generate(child[0], App.I(config, "CalibrationGroups")), test = Generate(child[1], App.I(config, "TestGroups"));
            var matrix = new double[calibration.Length, delays.Length];
            for (int i = 0; i < calibration.Length; i++) for (int j = 0; j < delays.Length; j++)
            {
                matrix[i, j] = Evaluate(calibration[i], delays[j]).Loss;
                App.Require(matrix[i, j] >= 0 && matrix[i, j] <= 1 && (j == 0 || matrix[i, j] <= matrix[i, j - 1]), "Bounded antitone recorded loss.");
            }
            RiskThreshold result = RiskControl.Mean(matrix, 1, alpha); int expected = delays.Length - 1; bool met = false; double last = 0;
            for (int j = 0; j < delays.Length; j++) { double sum = 1; for (int i = 0; i < calibration.Length; i++) sum += matrix[i, j]; last = sum / (calibration.Length + 1); if (last <= alpha) { expected = j; met = true; break; } }
            App.Require(result.Index == expected && result.Met == met && result.Reading == last, "Independent risk threshold.");
            double risk = 0, startRate = 0; int localStarts = 0; double localDelay = 0;
            for (int i = 0; i < test.Length; i++)
            {
                var e = Evaluate(test[i], delays[result.Index]); risk += e.Loss / test.Length; startRate += e.StartRate / test.Length;
                localStarts += e.Starts; localDelay += e.TotalDelay; emptyCount += test[i].Prefix.Length == 0 ? 1 : 0;
                App.Line(raw, new { Split = rep + 1, Conversation = i, test[i].Prefix, test[i].Budget, test[i].Completion,
                    Policy = result.Index, e.Loss, e.StartRate, e.Starts });
            }
            risks.Add(risk); starts.Add(startRate); selected.Add(result.Index); sumDelays += localDelay; countStarts += localStarts;
            App.Line(output, new { Split = rep + 1, Seed = seeds[rep], Policy = result.Index, result.Reading, result.Met,
                TestRisk = risk, TestStartRate = startRate, Starts = localStarts, DelaySum = localDelay });
        }
        App.Write("risk-summary.json", new { Passed = true, ExpectedEligibleTurnLoss = App.MeanInterval(risks.ToArray()),
            EqualConversationStartFraction = App.MeanInterval(starts.ToArray()), PooledDelayAmongStarts = countStarts > 0 ? sumDelays / countStarts : (double?)null,
            TotalStarts = countStarts, EmptyTestConversations = emptyCount, PolicyCounts = selected.GroupBy(x => x).Select(g => new { Policy = g.Key, Count = g.Count() }).ToArray(),
            EqualityControl = equality, FallbackControl = fallback,
            DenominatorWitness = new { EligibleLossBefore = before.Loss, EligibleLossAfter = after.Loss, AmongStartEarlyBefore = before.Loss / before.StartRate, AmongStartEarlyAfter = after.Loss / after.StartRate,
                AllAbstainEligibleLoss = 0, EmptyEligibleLoss = 0 } });
    }
}

internal static class Alarms
{
    private sealed class HalfBet : IBettingFunction
    {
        public double LogBet(double p) => p <= 0.5 ? Math.Log(2) : double.NegativeInfinity;
        public void Observe(double p) { }
        public bool IntegratesToOne => true;
    }
    internal static void Run(int seed)
    {
        var config = App.Settings.GetProperty("Alarm"); int horizon = App.I(config, "Horizon"), changeAt = App.I(config, "ChangeAt");
        int[] masters = SeedTree.Derive(seed, 4); var summaries = new List<object>();
        using var file = new StreamWriter(Path.Combine(App.Output, "alarm-streams.jsonl"));
        foreach (string law in new[] { "continuous-null", "tied-null", "all-ties-null", "noise-change" })
        {
            int index = Array.IndexOf(new[] { "continuous-null", "tied-null", "all-ties-null", "noise-change" }, law);
            int n = App.I(config, law == "noise-change" ? "ChangeStreams" : "Streams"), alarms = 0, pre = 0, censored = 0;
            var delays = new List<int>(); int[] seeds = SeedTree.Derive(masters[index], n);
            for (int i = 0; i < n; i++)
            {
                int[] child = SeedTree.Derive(seeds[i], 2); var rng = new Xoshiro256PlusPlus(child[0]);
                var alarm = new MartingaleAlarm(AlarmProcedure.Ville, App.D(config, "Threshold"));
                var monitor = new ConformalTestMartingale(new PowerBet(App.D(config, "PowerEpsilon")), child[1], alarm);
                int? first = null;
                for (int t = 1; t <= horizon; t++)
                {
                    double score = law == "all-ties-null" ? 0 : law == "tied-null" ? rng.NextInt(2) :
                        -Math.Abs(Normal.SampleStandard(rng)) * (law == "noise-change" && t >= changeAt ? App.D(config, "ChangedNoiseSd") : 1);
                    monitor.Add(score);
                    if (alarm.AlarmTimes.Count > 0) { first = t; break; }
                }
                if (first is int time) { alarms++; if (time < changeAt) pre++; else delays.Add(time - changeAt); } else censored++;
                App.Line(file, new { Law = law, Stream = i, Seed = seeds[i], FirstAlarm = first, LogCapital = monitor.LogCapital });
            }
            summaries.Add(new { Law = law, EverAlarm = App.Interval(alarms, n), PreChangeAlarms = law == "noise-change" ? pre : (int?)null,
                DetectedAfterChange = law == "noise-change" ? delays.Count : (int?)null, Censored = censored,
                MedianDetectedDelay = delays.Count > 0 ? 0.5 * (delays.Order().ElementAt((delays.Count - 1) / 2) + delays.Order().ElementAt(delays.Count / 2)) : (double?)null });
        }
        var negative = new List<object>();
        foreach (double p in new[] { 0.25, 0.75 }) { var capital = new BettingCapital(new HalfBet()); capital.Add(p); capital.Add(p); double final = Math.Exp(capital.LogCapital);
            App.Require(final == (p == 0.25 ? 4 : 0), "Repeated-uniform capital witness."); negative.Add(new { P = p, Capital = final, Probability = 0.5 }); }
        App.Write("alarm-summary.json", new { Passed = true, Horizon = horizon, Threshold = App.D(config, "Threshold"),
            VilleBound = 1 / App.D(config, "Threshold"), Streams = summaries, DependentUniformControl = negative, NegativeControlCrossingProbability = 0.5 });
    }
}
