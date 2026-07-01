using System;
using Maths.Rng;
using Maths.LinAlg;

namespace Maths.Geometry.DimReduction
{
    public delegate double SubspaceCostFunction(double[][] projectionMatrix);

    public static class Spred
    {
        /// <summary>
        /// Shape-Preserving Dimensionality Reduction (SPRED).
        /// Finds a linear projection matrix that minimizes a given topological cost function 
        /// (typically Wasserstein distance between persistent homology diagrams) using
        /// Simulated Annealing on the Stiefel Manifold.
        /// Based on Kisung You's 2106.02096.
        /// </summary>
        /// <param name="data">High dimensional data</param>
        /// <param name="targetDim">Target lower dimension (e.g., 2 or 3 for visualization)</param>
        /// <param name="costFunction">A function that takes a candidate projection matrix and returns its topological cost</param>
        /// <param name="maxIters">Max iterations for simulated annealing</param>
        /// <param name="seed">Optional RNG seed for reproducible annealing; null draws OS entropy.</param>
        /// <returns>The optimized targetDim x originalDim projection matrix</returns>
        public static double[][] Compute(double[][] data, int targetDim, SubspaceCostFunction costFunction, int maxIters = 1000, int? seed = null)
        {
            int n = data.Length;
            if (n == 0) throw new ArgumentException("Empty data");
            int d = data[0].Length;

            // 1. Initialize with PCA (to start with high variance explained)
            var pca = Pca.Compute(data, targetDim, center: true, whiten: false);
            double[][] currentProj = pca.Components; // targetDim x d
            double currentCost = costFunction(currentProj);
            
            double[][] bestProj = CopyMatrix(currentProj);
            double bestCost = currentCost;

            var rng = new Xoshiro256PlusPlus(seed);
            double initialTemp = 1.0;

            for (int iter = 0; iter < maxIters; iter++)
            {
                double temp = initialTemp * Math.Pow(0.99, iter); // Cooling schedule

                // 2. Propose a new projection matrix on the Stiefel Manifold
                double[][] proposal = GenerateStiefelProposal(currentProj, rng, stepSize: temp * 0.1);

                // 3. Evaluate cost
                double proposalCost = costFunction(proposal);

                // 4. Acceptance probability
                if (proposalCost < currentCost || rng.NextDouble() < Math.Exp((currentCost - proposalCost) / temp))
                {
                    currentProj = proposal;
                    currentCost = proposalCost;

                    if (currentCost < bestCost)
                    {
                        bestProj = CopyMatrix(currentProj);
                        bestCost = currentCost;
                    }
                }
            }

            return bestProj;
        }

        private static double[][] GenerateStiefelProposal(double[][] current, Xoshiro256PlusPlus rng, double stepSize)
        {
            int k = current.Length;
            int d = current[0].Length;

            // Perturb each frame vector, packing into a column-major block [row0 | row1 | ...]
            // so each length-d Stiefel vector lands as one column for the Orthonormalize kernel.
            double[] block = new double[k * d];
            for (int i = 0; i < k; i++)
                for (int j = 0; j < d; j++)
                    block[i * d + j] = current[i][j] + (rng.NextDouble() * 2 - 1) * stepSize;

            // Re-impose the Stiefel constraint with the audited double-pass modified Gram-Schmidt
            // kernel, replacing a single-pass loop that shed orthogonality over long annealing runs.
            MatrixOps.Orthonormalize(block, d, k);

            double[][] proposal = new double[k][];
            for (int i = 0; i < k; i++)
            {
                proposal[i] = new double[d];
                Array.Copy(block, i * d, proposal[i], 0, d);
            }
            return proposal;
        }

        private static double[][] CopyMatrix(double[][] matrix)
        {
            double[][] copy = new double[matrix.Length][];
            for (int i = 0; i < matrix.Length; i++)
            {
                copy[i] = new double[matrix[i].Length];
                Array.Copy(matrix[i], copy[i], matrix[i].Length);
            }
            return copy;
        }
    }
}
