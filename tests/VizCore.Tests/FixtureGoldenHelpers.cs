using System;
using System.Collections.Generic;
using System.Linq;
using Graphs;
using Graphs.Primitives;
using Graphs.Proximity;
using Graphs.Spectral;
using Maths.Geometry;
using Maths.LinAlg;
using TDA.Ph;
using Viz;

namespace VizCore.Tests;

internal sealed record DraftSpectralFixture(
    double[][] Points,
    CsrGraph Graph,
    IReadOnlyList<EigenPair> Eigenpairs,
    int[] ComponentLabels,
    int ComponentCount,
    string ExtractionPath,
    int FiedlerEigenIndex,
    int SecondEigenEigenIndex,
    NodeSignalLayer Fiedler,
    NodeSignalLayer SecondEigen,
    LineFieldLayer FiedlerGradient,
    double[] ArcParameter,
    double[] LoopParameter);

internal sealed record MobiusModeSelection(
    IReadOnlyList<int> CandidateEigenIndices,
    IReadOnlyList<int> CandidateSignFlipCounts,
    int SelectedEigenIndex,
    int SelectedSignFlipCount,
    double[] SelectedValues);

internal static class FixtureGoldenHelpers
{
    public const int SpectralSeed = 12345;
    public const int ConnectedFixtureK = 40;

    public static CsrGraph BuildDraftGraph(
        double[][] points,
        int k,
        bool ensureConnected = true)
    {
        int n = points.Length;
        Func<int, int, double> pairDistance = (i, j) => EuclideanDistance(points[i], points[j]);

        NeighborSelection selection = ProximityGraph.SelectKnn(n, k, pairDistance);
        if (ensureConnected)
            selection = ConnectivityRepair.EnsureConnected(selection, n, pairDistance);

        return BuildCombinatorialGraph(selection, n);
    }

    private static CsrGraph BuildCombinatorialGraph(NeighborSelection selection, int nodeCount)
    {
        var edges = new List<Edge>();
        var seen = new HashSet<long>();

        for (int source = 0; source < nodeCount; source++)
        {
            foreach (Neighbor neighbor in selection.AllNeighbors[source])
            {
                int lo = Math.Min(source, neighbor.Index);
                int hi = Math.Max(source, neighbor.Index);
                long key = ((long)lo * nodeCount) + hi;
                if (!seen.Add(key))
                    continue;

                edges.Add(new Edge(lo, hi, 1.0));
            }
        }

        return CsrGraph.FromEdges(edges.ToArray(), nodeCount);
    }

    public static DraftSpectralFixture RunDraftSpectral(
        double[][] points,
        int eigenCount,
        int k,
        int seed,
        bool ensureConnected = true,
        LaplacianType laplacianType = LaplacianType.Combinatorial)
    {
        CsrGraph graph = BuildDraftGraph(points, k, ensureConnected);

        return RunDraftSpectral(points, graph, eigenCount, seed, laplacianType);
    }

    public static DraftSpectralFixture RunDraftSpectral(
        double[][] points,
        CsrGraph graph,
        int eigenCount,
        int seed,
        LaplacianType laplacianType = LaplacianType.Combinatorial)
    {
        int[] componentLabels = BuildComponentLabels(graph, out int componentCount);

        IReadOnlyList<EigenPair> eigenpairs;
        string extractionPath;

        if (componentCount == 1)
        {
            eigenpairs = Spectral.ComputeBottomK(
                graph,
                seed: seed,
                k: Math.Max(eigenCount + componentCount + 1, 8),
                lapType: laplacianType,
                solverKind: SolverKind.Auto);
            extractionPath = "spectral-bottomk";
        }
        else
        {
            eigenpairs = ComputeDenseEigenpairs(graph, laplacianType);
            extractionPath = "helper-dense-global";
        }

        int[] nonTrivialIndices = SelectModesAfterNullSpace(eigenpairs.Count, componentCount, eigenCount);

        var nonTrivial = nonTrivialIndices.Select(index => eigenpairs[index]).ToArray();
        int fiedlerEigenIndex = nonTrivialIndices.Length > 0 ? nonTrivialIndices[0] : -1;
        int secondEigenIndex = nonTrivialIndices.Length > 1 ? nonTrivialIndices[1] : -1;

        double[] fiedler = nonTrivial.Length > 0
            ? (double[])nonTrivial[0].Vector.Clone()
            : new double[points.Length];

        double[] second = nonTrivial.Length > 1
            ? (double[])nonTrivial[1].Vector.Clone()
            : new double[points.Length];

        double[] gradient = fiedlerEigenIndex >= 0
            ? SpectralBridge.LineFieldFromEigenvector(graph, points, fiedler)
            : new double[points.Length * points[0].Length];

        return new DraftSpectralFixture(
            Points: points,
            Graph: graph,
            Eigenpairs: eigenpairs,
            ComponentLabels: componentLabels,
            ComponentCount: componentCount,
            ExtractionPath: extractionPath,
            FiedlerEigenIndex: fiedlerEigenIndex,
            SecondEigenEigenIndex: secondEigenIndex,
            Fiedler: new NodeSignalLayer("fiedler", fiedler, ScalarSource.Eigenfunction),
            SecondEigen: new NodeSignalLayer("eigen-2", second, ScalarSource.Eigenfunction),
            FiedlerGradient: new LineFieldLayer("fiedler-gradient", gradient, points.Length, points[0].Length, LineFieldSource.SpectralGradient),
            ArcParameter: BuildAngularParameter(points),
            LoopParameter: BuildAngularParameter(points));
    }

    private static IReadOnlyList<EigenPair> ComputeDenseEigenpairs(CsrGraph graph, LaplacianType laplacianType)
    {
        double[,] laplacian = BuildDenseLaplacian(graph, laplacianType);
        EigenResult eig = Eigen.DecomposeSymmetric(laplacian);
        var pairs = new List<EigenPair>(eig.Eigenvalues.Length);

        for (int i = 0; i < eig.Eigenvalues.Length; i++)
        {
            var vector = new double[graph.NodeCount];
            Array.Copy(eig.Eigenvectors[i], vector, graph.NodeCount);
            pairs.Add(new EigenPair(eig.Eigenvalues[i], vector));
        }

        pairs.Sort(static (left, right) => left.Lambda.CompareTo(right.Lambda));
        return pairs;
    }

    private static double[,] BuildDenseLaplacian(CsrGraph graph, LaplacianType lapType)
    {
        int n = graph.NodeCount;
        var laplacian = new double[n, n];
        var degree = new double[n];

        for (int i = 0; i < n; i++)
        {
            int start = graph.RowPointers[i];
            int end = graph.RowPointers[i + 1];
            double sum = 0.0;
            for (int edge = start; edge < end; edge++)
                sum += graph.Weights[edge] > 0.0 ? graph.Weights[edge] : 1.0;

            degree[i] = sum;
        }

        if (lapType == LaplacianType.NormalizedSymmetric)
        {
            var invSqrtDegree = new double[n];
            for (int i = 0; i < n; i++)
                invSqrtDegree[i] = degree[i] > 1e-12 ? 1.0 / Math.Sqrt(degree[i]) : 0.0;

            for (int i = 0; i < n; i++)
            {
                laplacian[i, i] = 1.0;
                if (invSqrtDegree[i] == 0.0)
                    continue;

                int start = graph.RowPointers[i];
                int end = graph.RowPointers[i + 1];
                for (int edge = start; edge < end; edge++)
                {
                    int target = graph.Targets[edge];
                    double weight = graph.Weights[edge];
                    laplacian[i, target] -= invSqrtDegree[i] * weight * invSqrtDegree[target];
                }
            }

            return laplacian;
        }

        for (int i = 0; i < n; i++)
        {
            laplacian[i, i] = degree[i];
            int start = graph.RowPointers[i];
            int end = graph.RowPointers[i + 1];
            for (int edge = start; edge < end; edge++)
            {
                int target = graph.Targets[edge];
                laplacian[i, target] -= graph.Weights[edge];
            }
        }

        return laplacian;
    }

    public static NodeSignalLayer[] BuildPerComponentFiedlers(
        double[][] points,
        CsrGraph graph,
        int eigenCount,
        int seed,
        LaplacianType laplacianType = LaplacianType.Combinatorial)
    {
        int[] labels = BuildComponentLabels(graph, out int componentCount);
        var lifted = new NodeSignalLayer[componentCount];

        for (int component = 0; component < componentCount; component++)
        {
            bool[] mask = labels.Select(label => label == component).ToArray();
            CsrGraph subgraph = graph.InducedSubgraph(mask, out int[] newToOld, out _);
            double[][] subPoints = newToOld.Select(index => points[index]).ToArray();
            var subFixture = RunDraftSpectral(subPoints, subgraph, eigenCount, seed, laplacianType);

            var full = new double[points.Length];
            double[] local = subFixture.Fiedler.Values.ToArray();
            for (int i = 0; i < newToOld.Length; i++)
                full[newToOld[i]] = local[i];

            lifted[component] = new NodeSignalLayer($"component-{component}-fiedler", full, ScalarSource.Eigenfunction);
        }

        return lifted;
    }

    private static int[] SelectModesAfterNullSpace(int eigenpairCount, int componentCount, int eigenCount)
    {
        if (eigenpairCount <= componentCount) return Array.Empty<int>();

        int take = Math.Min(eigenCount, eigenpairCount - componentCount);
        return Enumerable.Range(componentCount, take).ToArray();
    }

    public static int CountZeroModes(IReadOnlyList<EigenPair> eigenpairs, double tolerance = 1e-7) =>
        eigenpairs.Count(pair => Math.Abs(pair.Lambda) < tolerance);

    public static string FormatFirstEigenvalues(IReadOnlyList<EigenPair> eigenpairs, int count = 10)
    {
        int take = Math.Min(count, eigenpairs.Count);
        return string.Join(", ", Enumerable.Range(0, take).Select(i => $"{i}:{eigenpairs[i].Lambda:G17}"));
    }

    public static string FormatLastEigenvalues(IReadOnlyList<EigenPair> eigenpairs, int count = 10)
    {
        int take = Math.Min(count, eigenpairs.Count);
        int start = Math.Max(0, eigenpairs.Count - take);
        return string.Join(", ", Enumerable.Range(start, eigenpairs.Count - start).Select(i => $"{i}:{eigenpairs[i].Lambda:G17}"));
    }

    public static string FormatComponentSummary(int[] componentLabels)
    {
        return string.Join(", ",
            componentLabels
                .GroupBy(label => label)
                .OrderBy(group => group.Key)
                .Select(group => $"component {group.Key} -> [{group.First()}] count={group.Count()}"));
    }

    public static string FormatLargestComponentSizes(CsrGraph graph, int take = 3)
    {
        int[] labels = BuildComponentLabels(graph, out _);
        int[] sizes = labels
            .GroupBy(label => label)
            .Select(group => group.Count())
            .OrderByDescending(size => size)
            .Take(take)
            .ToArray();

        return sizes.Length == 0 ? string.Empty : string.Join(",", sizes);
    }

    public static int GetUndirectedEdgeCount(CsrGraph graph) => graph.Targets.Length / 2;

    public static double GetMinimumEdgeWeight(CsrGraph graph) =>
        graph.Weights.Length == 0 ? double.NaN : graph.Weights.Min();

    public static string FormatOrderedSigns(double[] values, int[] ordered, double epsilon = 1e-6)
    {
        return string.Join(",", ordered.Select(index => SignToken(values[index], epsilon)));
    }

    public static double MaxAbsOnComponent(double[] values, int[] componentLabels, int component)
    {
        double max = 0.0;
        for (int i = 0; i < values.Length; i++)
        {
            if (componentLabels[i] != component) continue;
            max = Math.Max(max, Math.Abs(values[i]));
        }

        return max;
    }

    public static double MaxAbsOffComponent(double[] values, int[] componentLabels, int component)
    {
        double max = 0.0;
        for (int i = 0; i < values.Length; i++)
        {
            if (componentLabels[i] == component) continue;
            max = Math.Max(max, Math.Abs(values[i]));
        }

        return max;
    }

    public static double SpearmanRank(double[] left, double[] right)
    {
        if (left.Length != right.Length)
            throw new ArgumentException("Input vectors must have the same length.");

        double[] leftRank = BuildRanks(left);
        double[] rightRank = BuildRanks(right);
        return Pearson(leftRank, rightRank);
    }

    public static int[] BinHistogram(double[] values, int bins)
    {
        if (values.Length == 0) return Array.Empty<int>();
        if (bins <= 0) throw new ArgumentOutOfRangeException(nameof(bins));

        double min = values.Min();
        double max = values.Max();
        var histogram = new int[bins];

        if (max - min <= 1e-12)
        {
            histogram[0] = values.Length;
            return histogram;
        }

        double scale = bins / (max - min);
        foreach (double value in values)
        {
            int bucket = (int)((value - min) * scale);
            if (bucket == bins) bucket--;
            histogram[bucket]++;
        }

        return histogram;
    }

    public static MobiusModeSelection SelectMobiusModeWithMinimumSignFlips(
        DraftSpectralFixture draft,
        int candidateCount = 5,
        double epsilon = 1e-6)
    {
        int[] ordered = Enumerable.Range(0, draft.Points.Length).OrderBy(i => draft.LoopParameter[i]).ToArray();
        int[] candidateIndices = SelectModesAfterNullSpace(draft.Eigenpairs.Count, draft.ComponentCount, candidateCount);

        if (candidateIndices.Length == 0)
        {
            return new MobiusModeSelection(
                CandidateEigenIndices: Array.Empty<int>(),
                CandidateSignFlipCounts: Array.Empty<int>(),
                SelectedEigenIndex: -1,
                SelectedSignFlipCount: int.MaxValue,
                SelectedValues: Array.Empty<double>());
        }

        var candidateFlipCounts = new int[candidateIndices.Length];
        int bestOffset = 0;

        for (int i = 0; i < candidateIndices.Length; i++)
        {
            double[] values = draft.Eigenpairs[candidateIndices[i]].Vector;
            candidateFlipCounts[i] = CountOrderedSignFlips(values, ordered, epsilon);

            if (candidateFlipCounts[i] < candidateFlipCounts[bestOffset])
                bestOffset = i;
        }

        int selectedEigenIndex = candidateIndices[bestOffset];
        return new MobiusModeSelection(
            CandidateEigenIndices: candidateIndices,
            CandidateSignFlipCounts: candidateFlipCounts,
            SelectedEigenIndex: selectedEigenIndex,
            SelectedSignFlipCount: candidateFlipCounts[bestOffset],
            SelectedValues: (double[])draft.Eigenpairs[selectedEigenIndex].Vector.Clone());
    }

    public static string FormatCandidateModeFlipCounts(MobiusModeSelection selection)
    {
        return string.Join(", ",
            selection.CandidateEigenIndices
                .Zip(selection.CandidateSignFlipCounts, static (index, count) => $"{index}:{count}"));
    }

    private static int CountOrderedSignFlips(double[] values, int[] ordered, double epsilon)
    {
        var signs = new List<int>(ordered.Length);

        for (int j = 0; j < ordered.Length; j++)
        {
            double value = values[ordered[j]];
            if (Math.Abs(value) < epsilon)
                continue;

            int sign = Math.Sign(value);
            if (signs.Count == 0 || signs[^1] != sign)
                signs.Add(sign);
        }

        if (signs.Count <= 1)
            return 0;

        int cyclicCrossings = 0;
        for (int j = 0; j < signs.Count; j++)
        {
            if (signs[j] != signs[(j + 1) % signs.Count])
                cyclicCrossings++;
        }

        return Math.Max(0, cyclicCrossings - 1);
    }

    public static double[] LocalChordFromNeighbors(double[][] points, CsrGraph graph, int index)
    {
        int start = graph.RowPointers[index];
        int end = graph.RowPointers[index + 1];
        int dim = points[index].Length;
        var chord = new double[dim];

        for (int edge = start; edge < end; edge++)
        {
            int target = graph.Targets[edge];
            for (int d = 0; d < dim; d++)
                chord[d] += points[target][d] - points[index][d];
        }

        NormalizeInPlace(chord);
        return chord;
    }

    public static double UnorientedAngleBetween(double[] first, double[] second)
    {
        double denom = Norm(first) * Norm(second);
        if (denom <= 1e-12) return 0.0;
        double cos = Math.Abs(Dot(first, second) / denom);
        cos = Math.Clamp(cos, -1.0, 1.0);
        return Math.Acos(cos);
    }

    public static double UnorientedAngleSigned(double[] first, double[] second)
    {
        double denom = Norm(first) * Norm(second);
        if (denom <= 1e-12) return 0.0;

        double dot = Dot(first, second) / denom;
        double altDot = Dot(first, Negate(second)) / denom;
        dot = Math.Clamp(dot, -1.0, 1.0);
        altDot = Math.Clamp(altDot, -1.0, 1.0);

        double best = Math.Abs(dot) >= Math.Abs(altDot) ? dot : altDot;
        return Math.Acos(Math.Clamp(best, -1.0, 1.0));
    }

    public static double AverageNormalComponent(LineFieldLayer field, double[][] points, CsrGraph graph)
    {
        double total = 0.0;
        for (int i = 0; i < points.Length; i++)
        {
            double[] direction = GetDirection(field, i);
            double[] chord = LocalChordFromNeighbors(points, graph, i);
            double parallel = Dot(direction, chord);
            total += Math.Abs(Math.Sqrt(Math.Max(0.0, Dot(direction, direction) - parallel * parallel)));
        }

        return points.Length == 0 ? 0.0 : total / points.Length;
    }

    public static double AverageRadialAlignment(LineFieldLayer field, double[][] points)
    {
        if (points.Length == 0)
            return 0.0;

        int dim = points[0].Length;
        var centroid = new double[dim];
        for (int i = 0; i < points.Length; i++)
            for (int d = 0; d < dim; d++)
                centroid[d] += points[i][d];

        double invCount = 1.0 / points.Length;
        for (int d = 0; d < dim; d++)
            centroid[d] *= invCount;

        double total = 0.0;
        int count = 0;
        for (int i = 0; i < points.Length; i++)
        {
            double[] direction = GetDirection(field, i);
            double[] radial = Subtract(points[i], centroid);
            double dirNorm = Norm(direction);
            double radialNorm = Norm(radial);
            if (dirNorm <= 1e-12 || radialNorm <= 1e-12)
                continue;

            total += Math.Abs(Dot(direction, radial) / (dirNorm * radialNorm));
            count++;
        }

        return count == 0 ? 0.0 : total / count;
    }

    public static double[] GetDirection(LineFieldLayer layer, int index)
    {
        double[] flat = layer.Directions.ToArray();
        var direction = new double[layer.D];
        Array.Copy(flat, index * layer.D, direction, 0, layer.D);
        return direction;
    }

    public static double[] GetDirection(double[] flat, int d, int index)
    {
        var direction = new double[d];
        Array.Copy(flat, index * d, direction, 0, d);
        return direction;
    }

    // Oracle for the geometric tangent on a planar arc: at arc parameter θ the tangent is (-sin θ, cos θ, 0, ...).
    // Uses generator ground truth (arc parameter) rather than the neighbor centroid, which is inward-radial on curved manifolds.
    // ambientDim must match the line-field dimension; extra dims beyond the XY plane are padded with zero.
    public static double[] ArcTangentFromParameter(double arcParam, int ambientDim)
    {
        var tangent = new double[ambientDim];
        tangent[0] = -Math.Sin(arcParam);
        if (ambientDim > 1) tangent[1] = Math.Cos(arcParam);
        return tangent;
    }

    private static double[] BuildAngularParameter(double[][] points)
    {
        var parameter = new double[points.Length];
        for (int i = 0; i < points.Length; i++)
        {
            double angle = Math.Atan2(points[i][1], points[i][0]);
            parameter[i] = angle < 0.0 ? angle + 2.0 * Math.PI : angle;
        }
        return parameter;
    }

    private static int[][] BuildAdjacency(CsrGraph graph)
    {
        var adjacency = new int[graph.NodeCount][];
        for (int i = 0; i < graph.NodeCount; i++)
        {
            int start = graph.RowPointers[i];
            int length = graph.RowPointers[i + 1] - start;
            var neighbors = new int[length];
            Array.Copy(graph.Targets, start, neighbors, 0, length);
            adjacency[i] = neighbors;
        }

        return adjacency;
    }

    private static int[] BuildComponentLabels(CsrGraph graph, out int componentCount)
    {
        int n = graph.NodeCount;
        var labels = Enumerable.Repeat(-1, n).ToArray();
        int next = 0;

        for (int root = 0; root < n; root++)
        {
            if (labels[root] >= 0) continue;
            var queue = new Queue<int>();
            queue.Enqueue(root);
            labels[root] = next;

            while (queue.Count > 0)
            {
                int node = queue.Dequeue();
                for (int edge = graph.RowPointers[node]; edge < graph.RowPointers[node + 1]; edge++)
                {
                    int target = graph.Targets[edge];
                    if (labels[target] >= 0) continue;
                    labels[target] = next;
                    queue.Enqueue(target);
                }
            }

            next++;
        }

        componentCount = next;
        return labels;
    }

    private static double[] BuildRanks(double[] values)
    {
        var ordered = values
            .Select((value, index) => (value, index))
            .OrderBy(static pair => pair.value)
            .ToArray();

        var ranks = new double[values.Length];
        int start = 0;
        while (start < ordered.Length)
        {
            int end = start + 1;
            while (end < ordered.Length && ordered[end].value == ordered[start].value)
                end++;

            double rank = (start + end - 1) * 0.5 + 1.0;
            for (int i = start; i < end; i++)
                ranks[ordered[i].index] = rank;

            start = end;
        }

        return ranks;
    }

    private static double Pearson(double[] left, double[] right)
    {
        double meanLeft = left.Average();
        double meanRight = right.Average();
        double num = 0.0;
        double denLeft = 0.0;
        double denRight = 0.0;

        for (int i = 0; i < left.Length; i++)
        {
            double dx = left[i] - meanLeft;
            double dy = right[i] - meanRight;
            num += dx * dy;
            denLeft += dx * dx;
            denRight += dy * dy;
        }

        double den = Math.Sqrt(denLeft * denRight);
        return den <= 1e-12 ? 0.0 : num / den;
    }

    private static double EuclideanDistance(double[] left, double[] right)
    {
        double sum = 0.0;
        for (int i = 0; i < left.Length; i++)
        {
            double diff = left[i] - right[i];
            sum += diff * diff;
        }

        return Math.Sqrt(sum);
    }

    private static double Dot(double[] left, double[] right)
    {
        double sum = 0.0;
        for (int i = 0; i < left.Length; i++)
            sum += left[i] * right[i];
        return sum;
    }

    private static double Norm(double[] vector) => Math.Sqrt(Dot(vector, vector));

    private static double[] Subtract(double[] left, double[] right)
    {
        var result = new double[left.Length];
        for (int i = 0; i < left.Length; i++)
            result[i] = left[i] - right[i];
        return result;
    }

    private static double[] Negate(double[] vector)
    {
        var result = new double[vector.Length];
        for (int i = 0; i < vector.Length; i++)
            result[i] = -vector[i];
        return result;
    }

    private static string SignToken(double value, double epsilon)
    {
        if (Math.Abs(value) < epsilon) return "0";
        return value > 0.0 ? "+" : "-";
    }

    private static void NormalizeInPlace(double[] vector)
    {
        double norm = Norm(vector);
        if (norm <= 1e-12) return;
        for (int i = 0; i < vector.Length; i++)
            vector[i] /= norm;
    }
}
