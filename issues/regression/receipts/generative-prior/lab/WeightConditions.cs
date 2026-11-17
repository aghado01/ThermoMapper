using System.Text.Json;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;

namespace GenerativePrior;

internal sealed record WeightSettings(int MasterSeed, int[] Dimensions, int[] Ratios, int Repetitions,
    int TrainingCases, double ActiveCases, double LocalProbability, double[] Radii, int Orientations,
    double[] Angles, double DecisionEpsilon, double ReferenceTolerance, double WitnessStandardErrors);

internal static class WeightConditions
{
    private sealed record Pair(double[] X, double[] Y, double Angle, string Name);
    private sealed record Reading(int K, int N, int Repetition, int Seed, string Law, double Radius,
        double GramGap, double MaxGap, string MaximizingPair, double WitnessDiagonal2,
        double WitnessGap, double MaxEigenvalueCheck, double MaxKnotCheck);

    private static double Dot(double[] a, double[] b) => PairwiseSum.Dot(a, b);

    private static double[] Unit(Xoshiro256PlusPlus rng, int k)
    {
        double[] x = Enumerable.Range(0, k).Select(_ => Normal.SampleStandard(rng)).ToArray();
        double norm = Math.Sqrt(Dot(x, x));
        return x.Select(v => v / norm).ToArray();
    }

    private static double[,] Target(Pair pair)
    {
        int k = pair.X.Length;
        var q = new double[k, k];
        double c = Math.Cos(pair.Angle), s = Math.Sin(pair.Angle);
        double[] v = s > 1e-12 ? pair.Y.Zip(pair.X, (y, x) => (y - c * x) / s).ToArray() : new double[k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++)
            {
                double swap = c * pair.X[i] * pair.X[j] +
                    s * (pair.X[i] * v[j] + v[i] * pair.X[j]) - c * v[i] * v[j];
                q[i, j] = (i == j ? (Math.PI - pair.Angle) / (2 * Math.PI) : 0) +
                    (s > 1e-12 ? s * swap / (2 * Math.PI) : 0);
            }
        return q;
    }

    private static (double Norm, double Error) Spectral(double[,] matrix, bool check)
    {
        EigenResult result = DenseEigen.DecomposeSymmetric(matrix);
        if (!result.Convergence.Converged) throw new ArithmeticException("Dense eigensolver did not converge.");
        double norm = result.Eigenvalues.Max(Math.Abs), error = 0;
        if (check)
        {
            double[] reference = Eigen.DecomposeSymmetric(matrix).Eigenvalues.Order().ToArray();
            double[] actual = result.Eigenvalues.Order().ToArray();
            error = actual.Zip(reference, (a, b) => Math.Abs(a - b)).Max();
            int k = matrix.GetLength(0);
            for (int v = 0; v < k; v++)
                for (int i = 0; i < k; i++)
                {
                    double applied = 0;
                    for (int j = 0; j < k; j++) applied += matrix[i, j] * result.Eigenvectors[v][j];
                    error = Math.Max(error, Math.Abs(applied - result.Eigenvalues[v] * result.Eigenvectors[v][i]));
                }
        }
        return (norm, error);
    }

    private static double[,] Gram(double[][] rows, int n, double scale, Pair? pair, double[] knots, double radius)
    {
        int k = rows[0].Length;
        var gram = new double[k, k];
        for (int r = 0; r < n; r++)
        {
            double[] w = rows[r];
            if (pair is not null && !(radius * Dot(w, pair.X) > knots[r] && radius * Dot(w, pair.Y) > knots[r])) continue;
            for (int i = 0; i < k; i++)
                for (int j = 0; j <= i; j++) gram[i, j] += scale * w[i] * w[j];
        }
        for (int i = 0; i < k; i++)
            for (int j = 0; j < i; j++) gram[j, i] = gram[i, j];
        return gram;
    }

    private static List<Pair> Pairs(WeightSettings s, int k, Xoshiro256PlusPlus rng)
    {
        var pairs = new List<Pair>();
        var e1 = new double[k]; e1[0] = 1;
        var e2 = new double[k]; e2[1] = 1;
        pairs.Add(new(e1, e1, 0, "e1-e1"));
        pairs.Add(new(e1, e2, Math.PI / 2, "e1-e2"));
        pairs.Add(new(e1, e1.Select(v => -v).ToArray(), Math.PI, "e1-minus-e1"));
        foreach (double angle in s.Angles)
            for (int p = 0; p < s.Orientations; p++)
            {
                double[] x = Unit(rng, k), v = Unit(rng, k);
                double projection = Dot(x, v);
                for (int j = 0; j < k; j++) v[j] -= projection * x[j];
                double norm = Math.Sqrt(Dot(v, v));
                for (int j = 0; j < k; j++) v[j] /= norm;
                double[] y = x.Zip(v, (a, b) => Math.Cos(angle) * a + Math.Sin(angle) * b).ToArray();
                if (angle == 0) y = (double[])x.Clone();
                if (angle == Math.PI) y = x.Select(a => -a).ToArray();
                pairs.Add(new(x, y, angle, $"angle-{angle:R}-orientation-{p}"));
            }
        return pairs;
    }

    public static int Run(WeightSettings s, string output)
    {
        int[] masters = SeedTree.Derive(s.MasterSeed, s.Dimensions.Length);
        var readings = new List<Reading>();
        using var controls = new StreamWriter(Path.Combine(output, "control-matrices.jsonl"));
        for (int d = 0; d < s.Dimensions.Length; d++)
        {
            int k = s.Dimensions[d], maximum = k * s.Ratios.Max();
            int[] seeds = SeedTree.Derive(masters[d], s.Repetitions);
            for (int rep = 0; rep < s.Repetitions; rep++)
            {
                int[] branches = SeedTree.Derive(seeds[rep], 6);
                var designRng = new Xoshiro256PlusPlus(branches[4]);
                double[][] columns = Enumerable.Range(0, k).Select(_ =>
                    Enumerable.Range(0, s.TrainingCases).Select(_ => Normal.SampleStandard(designRng)).ToArray()).ToArray();
                var data = new BassData(columns.Select((v, j) => (BassPredictor)new ContinuousPredictor($"x{j}", v)).ToArray(),
                    new ScalarResponse(columns[0]));
                BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
                var space = new DirectionalTermSpace(inputs, Math.Min(3, k), 1 - s.ActiveCases / s.TrainingCases, s.LocalProbability);
                var gauss = new double[maximum][];
                var sphere = new double[maximum][];
                var sparse = new double[maximum][];
                var knots = new double[maximum];
                var gaussianRng = new Xoshiro256PlusPlus(branches[0]);
                var sphereRng = new Xoshiro256PlusPlus(branches[1]);
                var directionRng = new Xoshiro256PlusPlus(branches[2]);
                var knotRng = new Xoshiro256PlusPlus(branches[3]);
                double knotError = 0;
                for (int r = 0; r < maximum; r++)
                {
                    gauss[r] = Enumerable.Range(0, k).Select(_ => Normal.SampleStandard(gaussianRng)).ToArray();
                    sphere[r] = Unit(sphereRng, k);
                    int active = 1 + directionRng.NextInt(Math.Min(3, k));
                    int[] variables = Enumerable.Range(0, k).ToArray();
                    for (int j = 0; j < active; j++)
                    {
                        int index = j + directionRng.NextInt(k - j);
                        (variables[j], variables[index]) = (variables[index], variables[j]);
                    }
                    ProjectionDirection direction = space.DrawDirection(directionRng, variables[..active]);
                    sparse[r] = Enumerable.Range(0, k).Select(direction.WeightOn).ToArray();
                    KnotRange range = space.KnotRangeOf(direction) ?? throw new ArithmeticException("Missing knot range.");
                    knots[r] = DirectionalTermSpace.DrawKnot(knotRng, range);
                    double[] projection = new double[s.TrainingCases];
                    for (int i = 0; i < projection.Length; i++)
                        for (int j = 0; j < k; j++)
                        {
                            var metric = inputs.Standardizing(j)!.Value;
                            projection[i] += sparse[r][j] * (metric.Offset + metric.Slope * inputs.Cases.Continuous(j)[i]);
                        }
                    Array.Sort(projection);
                    double index7 = (projection.Length - 1) * space.KnotLevel;
                    int lower = (int)Math.Floor(index7);
                    double upper = (1 - (index7 - lower)) * projection[lower] + (index7 - lower) * projection[lower + 1];
                    double low = upper - (upper - projection[0]) / s.LocalProbability;
                    knotError = Math.Max(knotError, Math.Max(Math.Abs(upper - range.Upper), Math.Abs(low - range.Lower)));
                }
                List<Pair> pairs = Pairs(s, k, new Xoshiro256PlusPlus(branches[5]));
                var laws = new[] { "gaussian", "sphere", "sparse", "sparse-threshold" };
                foreach (int ratio in s.Ratios)
                    foreach (string law in laws)
                        foreach (double radius in law == "sparse-threshold" ? s.Radii : new[] { 1.0 })
                        {
                            int n = k * ratio;
                            double[][] rows = law == "gaussian" ? gauss : law == "sphere" ? sphere : sparse;
                            double scale = law == "gaussian" ? 1.0 / n : (double)k / n;
                            double[] biases = law == "sparse-threshold" ? knots : new double[maximum];
                            double[,] full = Gram(rows, n, scale, null, biases, radius);
                            for (int j = 0; j < k; j++) full[j, j] -= 1;
                            var fullEigen = Spectral(full, true);
                            double maximumGap = 0, diagonal = 0, witnessGap = 0, eigenError = fullEigen.Error;
                            string maximumPair = "";
                            foreach (Pair pair in pairs)
                            {
                                double[,] gated = Gram(rows, n, scale, pair, biases, radius);
                                if (pair.Name == "e1-e1") diagonal = gated[1, 1];
                                double[,] target = Target(pair);
                                for (int i = 0; i < k; i++)
                                    for (int j = 0; j < k; j++) gated[i, j] -= target[i, j];
                                var eigen = Spectral(gated, pair.Name.StartsWith("e1"));
                                eigenError = Math.Max(eigenError, eigen.Error);
                                if (pair.Name == "e1-e1") witnessGap = eigen.Norm;
                                if (eigen.Norm > maximumGap) { maximumGap = eigen.Norm; maximumPair = pair.Name; }
                                if (k == 8 && ratio == s.Ratios.Max() && pair.Name.StartsWith("e1"))
                                {
                                    double[][] matrix = Enumerable.Range(0, k).Select(i =>
                                        Enumerable.Range(0, k).Select(j => gated[i, j]).ToArray()).ToArray();
                                    controls.WriteLine(JsonSerializer.Serialize(new { k, n, Repetition = rep + 1, Law = law,
                                        radius, Pair = pair.Name, pair.X, pair.Y, Residual = matrix, eigen.Norm }));
                                }
                            }
                            if (eigenError > s.ReferenceTolerance * Math.Max(1, maximumGap) || knotError > 1e-10)
                                throw new ArithmeticException("Independent spectral or knot check failed.");
                            readings.Add(new(k, n, rep + 1, seeds[rep], law, radius, fullEigen.Norm, maximumGap,
                                maximumPair, diagonal, witnessGap, eigenError, knotError));
                        }
                Console.WriteLine($"G3 k={k}, repetition {rep + 1}/{s.Repetitions} completed.");
            }
        }
        File.WriteAllLines(Path.Combine(output, "readings.jsonl"), readings.Select(row => JsonSerializer.Serialize(row)));
        var groups = readings.GroupBy(r => new { r.K, r.N, r.Law, r.Radius }).Select(g => new {
            g.Key.K, g.Key.N, Ratio = g.Key.N / g.Key.K, g.Key.Law, g.Key.Radius,
            MedianGap = Median(g.Select(r => r.MaxGap)), MinimumGap = g.Min(r => r.MaxGap), MaximumGap = g.Max(r => r.MaxGap),
            MedianGramGap = Median(g.Select(r => r.GramGap)), MeanWitnessDiagonal2 = g.Average(r => r.WitnessDiagonal2),
            SampledEpsilonReached = Median(g.Select(r => r.MaxGap)) <= s.DecisionEpsilon
        }).ToArray();
        var witnesses = readings.Where(r => r.N == r.K * s.Ratios.Max() && r.Law != "sparse-threshold")
            .GroupBy(r => new { r.K, r.Law }).Select(g => {
                double mean = g.Average(r => r.WitnessDiagonal2);
                double se = Math.Sqrt(g.Sum(r => Math.Pow(r.WitnessDiagonal2 - mean, 2)) / (g.Count() - 1) / g.Count());
                double expected = g.Key.Law == "sparse" ? (Math.Min(3, g.Key.K) - 1.0) / (4 * (g.Key.K - 1)) : 0.5;
                return new { g.Key.K, g.Key.Law, Mean = mean, StandardError = se, Expected = expected,
                    Passed = Math.Abs(mean - expected) <= s.WitnessStandardErrors * se + 1e-10 };
            }).ToArray();
        LabProgram.Write(Path.Combine(output, "summary.json"), new {
            LabProgram.Revision, LabProgram.SettingsSha256, Groups = groups, Witnesses = witnesses,
            NumericalControlsPassed = true, WitnessGatesPassed = witnesses.All(w => w.Passed)
        });
        return witnesses.All(w => w.Passed) ? 0 : 1;
    }

    private static double Median(IEnumerable<double> values)
    {
        double[] sorted = values.Order().ToArray();
        return 0.5 * (sorted[(sorted.Length - 1) / 2] + sorted[sorted.Length / 2]);
    }
}
