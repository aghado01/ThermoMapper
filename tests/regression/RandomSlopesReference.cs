using System;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;

namespace Maths.Regression.Tests;

// Full level/case matrices: no production contrast basis, thin factorization or whitening.
internal static class RandomSlopesReference
{
    internal static (BassConfig Config, BassData Data) Fixture(RandomSlopeScale scaling, bool intercept = true, bool drawn = true)
    {
        int[] counts = [3,4,5,6,7];
        int[] codes = [.. counts.SelectMany((count,g) => Enumerable.Repeat(g+1,count))];
        double[] x = [.. Enumerable.Range(0,codes.Length).Select(i => ((11*i)%31+0.25)/32)];
        double[] y = [.. x.Select((v,i) => 0.8+1.2*v+0.4*Math.Sin(1.7*i)+0.2*codes[i])];
        var data = new BassData([new ContinuousPredictor("x",x,new InputDomain(0,1))], new ScalarResponse(y),
            [new CategoricalPredictor("g",codes,["unused","a","b","c","d","e"])]);
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxBasis = 1, MaxInteraction = 1, Legendre = new LegendreBasis { MaxDegree = 2 } },
            Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(0.7)),
            Count = new BassCountPrior(1.5,1), Intercept = GPriorCentre.Zero, Sigma2 = Sigma2Prior.InverseGamma(3,2),
            RandomIntercepts = intercept ? [new("g",0.9,ScalePrior: drawn ? new(0.8) : null)] : null,
            RandomSlopes = [new("g","x",scaling,0.6,drawn ? new(1.2) : null)],
            Sampling = new BassSampling { Chains = 1, BurnIn = 500, Samples = 500, MasterSeed = 210701 },
        };
        return (config,data);
    }

    internal static double[] Input(BassData data, BassRandomSlope slope)
    {
        ContinuousPredictor input = data.Predictors.OfType<ContinuousPredictor>().Single(p => p.Name == slope.Input);
        double[] x = input.Values.ToArray();
        double min = input.Domain?.Lower ?? x.Min(), max = input.Domain?.Upper ?? x.Max();
        double[] unit = [.. x.Select(v => (v-min)/(max-min))];
        if (slope.InputScale == RandomSlopeScale.Unit) return unit;
        double mean = unit.Average(), deviation = Math.Sqrt(unit.Select(v => (v-mean)*(v-mean)).Average());
        return [.. unit.Select(v => (v-mean)/deviation)];
    }

    internal static Block[] Blocks(BassModel model, double[]? scales = null)
    {
        scales ??= model.RandomEffects!.Scales.ToArray();
        return [.. model.RandomEffects!.Specifications.Select((spec,k) => spec is BassRandomSlope slope
            ? Slope(model.Data.Groupings.Single(g => g.Name == slope.Grouping).Codes.ToArray(), Input(model.Data,slope), scales[k])
            : Intercept(model.Data.Groupings.Single(g => g.Name == spec.Grouping).Codes.ToArray(), scales[k]))];
    }

    internal sealed record Block(double[,] Design, double[,] Projector, double Scale);

    internal static Block Slope(int[] codes, double[] z, double scale)
    {
        int[] levels = [.. codes.Distinct().Order()];
        var design = new double[codes.Length, levels.Length];
        var first = new double[levels.Length];
        var second = new double[levels.Length];
        for (int i = 0; i < codes.Length; i++)
        {
            int g = Array.IndexOf(levels, codes[i]);
            design[i,g] = z[i];
            first[g] += z[i];
            second[g] += z[i] * z[i];
        }
        return new(design, Projector(first, second), scale);
    }

    internal static Block Intercept(int[] codes, double scale)
    {
        int[] levels = [.. codes.Distinct().Order()];
        var design = new double[codes.Length, levels.Length];
        var counts = new double[levels.Length];
        for (int i = 0; i < codes.Length; i++)
        {
            int g = Array.IndexOf(levels, codes[i]);
            design[i,g] = 1;
            counts[g]++;
        }
        return new(design, Projector(counts, new double[levels.Length]), scale);
    }

    // Subtract the row-space projector, using the analytic inverse of the two-row Gram.
    // A dependent second row uses the one-row projector; the zero design has no constraints.
    private static double[,] Projector(double[] first, double[] second)
    {
        int l = first.Length;
        double aa = Dot(first, first), bb = Dot(second, second), ab = Dot(first, second);
        double det = aa * bb - ab * ab;
        var p = new double[l,l];
        for (int g = 0; g < l; g++)
            for (int h = 0; h < l; h++)
            {
                double removed = det > 1e-12 * aa * bb
                    ? (bb*first[g]*first[h] - ab*(first[g]*second[h] + second[g]*first[h]) + aa*second[g]*second[h]) / det
                    : aa > 0 ? first[g]*first[h]/aa : bb > 0 ? second[g]*second[h]/bb : 0;
                p[g,h] = (g == h ? 1 : 0) - removed;
            }
        return p;
    }

    internal static double[,] Covariance(params Block[] blocks)
    {
        int n = blocks[0].Design.GetLength(0);
        var v = new double[n,n];
        for (int i = 0; i < n; i++) v[i,i] = 1;
        foreach (Block block in blocks)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    for (int g = 0; g < block.Projector.GetLength(0); g++)
                        for (int h = 0; h < block.Projector.GetLength(1); h++)
                            v[i,j] += block.Scale*block.Scale*block.Design[i,g]*block.Projector[g,h]*block.Design[j,h];
        return v;
    }

    internal static (double[][] Mean, double[,] Covariance) Conditional(Block[] blocks, double[] residual, double sigma2)
    {
        double[,] v = Covariance(blocks);
        var factor = new DenseCholesky(residual.Length);
        factor.Decompose(v);
        int[] offsets = new int[blocks.Length+1];
        for (int k = 0; k < blocks.Length; k++) offsets[k+1] = offsets[k] + blocks[k].Projector.GetLength(0);
        int l = offsets[^1], n = residual.Length;
        var cross = new double[l][];
        var prior = new double[l,l];
        for (int k = 0; k < blocks.Length; k++)
        {
            Block block = blocks[k];
            for (int g = 0; g < block.Projector.GetLength(0); g++)
            {
                int a = offsets[k]+g;
                cross[a] = new double[n];
                for (int h = 0; h < block.Projector.GetLength(1); h++)
                {
                    double covariance = block.Scale*block.Scale*block.Projector[g,h];
                    prior[a,offsets[k]+h] = covariance;
                    for (int i = 0; i < n; i++) cross[a][i] += covariance*block.Design[i,h];
                }
            }
        }
        double[] solved = factor.Solve(residual);
        double[][] means = [.. Enumerable.Range(0, blocks.Length).Select(k =>
            Enumerable.Range(offsets[k], offsets[k+1]-offsets[k]).Select(g => Dot(cross[g], solved)).ToArray())];
        var conditional = new double[l,l];
        for (int a = 0; a < l; a++)
        {
            double[] solution = factor.Solve(cross[a]);
            for (int b = 0; b < l; b++) conditional[a,b] = sigma2*(prior[a,b]-Dot(cross[b], solution));
        }
        return (means, conditional);
    }

    internal static double LogMarginal(BassModel model, BassTerm[] terms, Block[] blocks, double tau)
        => LogMarginal(model,terms,Covariance(blocks),tau);

    internal static double LogMarginal(BassModel model, BassTerm[] terms, double[,] covariance, double tau)
    {
        double[,] v = (double[,])covariance.Clone();
        int n = v.GetLength(0);
        var columns = new double[terms.Length+1][];
        columns[0] = [.. Enumerable.Repeat(1.0,n)];
        for (int j = 0; j < terms.Length; j++)
        {
            columns[j+1] = new double[n];
            model.Basis.CaseColumn(terms[j], model.Inputs.Cases, columns[j+1]);
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                foreach (double[] column in columns) v[i,j] += column[i]*column[j]/tau;
        var factor = new DenseCholesky(n);
        factor.Decompose(v);
        double[] y = model.Data.Response.Values.ToArray();
        double energy = Dot(y, factor.Solve(y));
        const double shape = 3, rate = 2;
        return SpecialFunctions.LogGamma(shape+0.5*n)-SpecialFunctions.LogGamma(shape) + shape*Math.Log(rate)
            -(shape+0.5*n)*Math.Log(rate+0.5*energy)-0.5*n*Math.Log(2*Math.PI)-0.5*factor.LogDet;
    }

    internal static double Dot(double[] a, double[] b) => a.Zip(b, (x,y) => x*y).Sum();
}
