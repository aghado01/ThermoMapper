using System;
using System.Collections.Generic;
using Clustering.Primitives;

namespace Clustering.Evaluation
{
    internal static class EvaluationHelpers
    {
        /// <summary>
        /// The assigned-only view of a point set: the rows whose label is not
        /// <see cref="Assignment.Unassigned"/> — MATLAB's <c>~Missing</c> subset.
        /// Internal indices score over this view so an unassigned point distorts
        /// neither centroids nor silhouettes. Returns the inputs unchanged when
        /// nothing is unassigned (fast path); otherwise compacted parallel copies
        /// whose rows are shared by reference (evaluators only read them).
        /// </summary>
        public static (double[][] Data, int[] Labels) AssignedSubset(double[][] data, int[] labels)
        {
            int assigned = CountAssigned(labels);
            if (assigned == labels.Length)
                return (data, labels);

            var dataOut = new double[assigned][];
            var labelsOut = new int[assigned];
            int k = 0;
            for (int i = 0; i < labels.Length; i++)
            {
                if (labels[i] == Assignment.Unassigned) continue;
                dataOut[k] = data[i];
                labelsOut[k] = labels[i];
                k++;
            }
            return (dataOut, labelsOut);
        }

        /// <summary>
        /// The assigned-only view of a predicted/reference label pair: the points
        /// whose <i>predicted</i> label is not <see cref="Assignment.Unassigned"/>.
        /// External indices score agreement over this subset (reference labels are
        /// ground truth, assumed fully labelled). Validates the pair, then compacts;
        /// returns the inputs unchanged when nothing is unassigned.
        /// </summary>
        public static (int[] Predicted, int[] Reference) AssignedByPredicted(int[] predicted, int[] reference)
        {
            ArgumentNullException.ThrowIfNull(predicted);
            ArgumentNullException.ThrowIfNull(reference);
            if (predicted.Length != reference.Length)
                throw new ArgumentException(
                    $"predicted length ({predicted.Length}) does not match reference length ({reference.Length}).");

            int assigned = CountAssigned(predicted);
            if (assigned == predicted.Length)
                return (predicted, reference);

            var predOut = new int[assigned];
            var refOut = new int[assigned];
            int k = 0;
            for (int i = 0; i < predicted.Length; i++)
            {
                if (predicted[i] == Assignment.Unassigned) continue;
                predOut[k] = predicted[i];
                refOut[k] = reference[i];
                k++;
            }
            return (predOut, refOut);
        }

        private static int CountAssigned(int[] labels)
        {
            int assigned = 0;
            for (int i = 0; i < labels.Length; i++)
                if (labels[i] != Assignment.Unassigned) assigned++;
            return assigned;
        }

        public static void ValidateInputs(double[][] data, int[] labels)
        {
            if (data is null)
                throw new ArgumentNullException(nameof(data));
            if (labels is null)
                throw new ArgumentNullException(nameof(labels));
            if (data.Length != labels.Length)
                throw new ArgumentException("Data and label arrays must have the same length.", nameof(labels));
            if (data.Length == 0)
                throw new ArgumentException("Data must contain at least one point.", nameof(data));

            int dim = data[0]?.Length ?? 0;
            if (dim == 0)
                throw new ArgumentException("Data rows must contain at least one dimension.", nameof(data));

            for (int i = 1; i < data.Length; i++)
            {
                if (data[i] is null)
                    throw new ArgumentException($"Data row {i} is null.", nameof(data));
                if (data[i].Length != dim)
                    throw new ArgumentException("All data rows must have the same dimensionality.", nameof(data));
            }
        }

        public static int[] MapLabelsToDense(int[] labels, out int[] uniqueLabels)
        {
            var labelMap = new Dictionary<int, int>(labels.Length);
            var mapped = new int[labels.Length];
            int next = 0;

            for (int i = 0; i < labels.Length; i++)
            {
                int label = labels[i];
                if (!labelMap.TryGetValue(label, out int mappedLabel))
                {
                    mappedLabel = next++;
                    labelMap[label] = mappedLabel;
                }
                mapped[i] = mappedLabel;
            }

            uniqueLabels = new int[next];
            foreach (var kvp in labelMap)
            {
                uniqueLabels[kvp.Value] = kvp.Key;
            }

            return mapped;
        }

        public static int[] CountClusters(int[] mappedLabels, int clusterCount)
        {
            var counts = new int[clusterCount];
            for (int i = 0; i < mappedLabels.Length; i++)
            {
                counts[mappedLabels[i]]++;
            }
            return counts;
        }

        public static double[][] ComputeCentroids(double[][] data, int[] mappedLabels, int clusterCount)
        {
            int n = data.Length;
            int dim = data[0].Length;
            var centroids = new double[clusterCount][];
            var counts = new int[clusterCount];

            for (int c = 0; c < clusterCount; c++)
            {
                centroids[c] = new double[dim];
            }

            for (int i = 0; i < n; i++)
            {
                int cluster = mappedLabels[i];
                counts[cluster]++;
                for (int d = 0; d < dim; d++)
                {
                    centroids[cluster][d] += data[i][d];
                }
            }

            for (int c = 0; c < clusterCount; c++)
            {
                if (counts[c] == 0)
                    continue;
                double inv = 1.0 / counts[c];
                for (int d = 0; d < dim; d++)
                {
                    centroids[c][d] *= inv;
                }
            }

            return centroids;
        }
    }
}
