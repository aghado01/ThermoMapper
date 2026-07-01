using System;
using Maths.LinAlg;

namespace TDA.Mapper.Diagnostics;

/// <summary>
/// Reports local intrinsic dimension per mapper node using PCA participation ratio,
/// $(\sum \lambda_i)^2 / \sum \lambda_i^2$, over the node's member feature vectors.
/// </summary>
public readonly record struct IntrinsicDimensionReport(
    int NodeCount,
    double[] PerNodeIntrinsicDim,
    double MedianDim,
    double MeanDim,
    int EffectivelyOneDimensionalNodes,
    int EffectivelyHigherDimNodes);

public static class IntrinsicDimension
{
    public static IntrinsicDimensionReport From(
        MapperResult result,
        double[][] originalPoints,
        double lowDimThreshold = 1.2,
        double highDimThreshold = 2.5)
    {
        ArgumentNullException.ThrowIfNull(originalPoints);

        int nodeCount = result.Nodes.Count;
        var perNode = new double[nodeCount];
        var finiteValues = new double[nodeCount];
        int finiteCount = 0;
        int effectivelyOneDimensional = 0;
        int effectivelyHigherDimensional = 0;
        double finiteSum = 0.0;

        for (int i = 0; i < nodeCount; i++)
        {
            var node = result.Nodes[i];
            if (node.MemberIndices.Length < 2)
            {
                perNode[i] = double.NaN;
                continue;
            }

            var subset = new double[node.MemberIndices.Length][];
            for (int j = 0; j < node.MemberIndices.Length; j++)
            {
                int memberIndex = node.MemberIndices[j];
                if ((uint)memberIndex >= (uint)originalPoints.Length)
                    throw new ArgumentOutOfRangeException(nameof(originalPoints),
                        $"Member index {memberIndex} is outside the original point array.");

                subset[j] = originalPoints[memberIndex];
            }

            double dimension = ComputeParticipationRatio(subset);
            perNode[i] = dimension;
            finiteValues[finiteCount++] = dimension;
            finiteSum += dimension;

            if (dimension < lowDimThreshold)
                effectivelyOneDimensional++;
            if (dimension > highDimThreshold)
                effectivelyHigherDimensional++;
        }

        double median = double.NaN;
        double mean = double.NaN;

        if (finiteCount > 0)
        {
            Array.Sort(finiteValues, 0, finiteCount);
            median = ComputeMedian(finiteValues, finiteCount);
            mean = finiteSum / finiteCount;
        }

        return new IntrinsicDimensionReport(
            NodeCount: nodeCount,
            PerNodeIntrinsicDim: perNode,
            MedianDim: median,
            MeanDim: mean,
            EffectivelyOneDimensionalNodes: effectivelyOneDimensional,
            EffectivelyHigherDimNodes: effectivelyHigherDimensional);
    }

    private static double ComputeParticipationRatio(double[][] data)
    {
        int dimension = data[0].Length;
        var pca = Pca.Compute(data, numComponents: dimension, center: true, whiten: false);

        double lambdaSum = 0.0;
        double lambdaSquaredSum = 0.0;
        for (int i = 0; i < pca.Eigenvalues.Length; i++)
        {
            double eigenvalue = pca.Eigenvalues[i];
            if (eigenvalue <= 0.0)
                continue;

            lambdaSum += eigenvalue;
            lambdaSquaredSum += eigenvalue * eigenvalue;
        }

        if (lambdaSquaredSum <= 0.0)
            return double.NaN;

        return (lambdaSum * lambdaSum) / lambdaSquaredSum;
    }

    private static double ComputeMedian(double[] sortedValues, int count)
    {
        int mid = count / 2;
        return (count & 1) == 0
            ? (sortedValues[mid - 1] + sortedValues[mid]) / 2.0
            : sortedValues[mid];
    }
}
