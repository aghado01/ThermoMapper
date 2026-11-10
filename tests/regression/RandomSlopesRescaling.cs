using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

// Symmetric, homogeneous relative comparisons. A vector's paired RMS supplies its own scale near zero;
// positive scalar parameters use their own magnitude. There is no absolute floor or physical-unit constant.
internal static class RandomSlopesRescaling
{
    internal static double Relative(ReadOnlySpan<double> left, ReadOnlySpan<double> right)
    {
        Assert.Equal(left.Length,right.Length);
        double magnitude = 0;
        for (int i = 0; i < left.Length; i++)
        {
            Assert.True(double.IsFinite(left[i]) && double.IsFinite(right[i]));
            magnitude = Math.Max(magnitude,Math.Max(Math.Abs(left[i]),Math.Abs(right[i])));
        }
        if (magnitude == 0) return 0; // both vectors exactly zero
        double squares = 0;
        for (int i = 0; i < left.Length; i++) squares += (left[i]/magnitude)*(left[i]/magnitude)+(right[i]/magnitude)*(right[i]/magnitude);
        double rms = magnitude*Math.Sqrt(squares/(2*left.Length)), worst = 0;
        for (int i = 0; i < left.Length; i++)
        {
            double scale = Math.Max(rms,Math.Max(Math.Abs(left[i]),Math.Abs(right[i])));
            worst = Math.Max(worst,Math.Abs(left[i]-right[i])/scale);
        }
        return worst;
    }

    internal static double Positive(double left,double right)
    {
        Assert.True(left > 0 && right > 0 && double.IsFinite(left) && double.IsFinite(right));
        return Math.Abs(left-right)/Math.Max(left,right);
    }

    internal static void Bits(ReadOnlySpan<double> left,ReadOnlySpan<double> right)
    {
        Assert.Equal(left.Length,right.Length);
        for (int i = 0; i < left.Length; i++) Assert.Equal(BitConverter.DoubleToInt64Bits(left[i]),BitConverter.DoubleToInt64Bits(right[i]));
    }

    internal sealed record Sweep(string Structure,MoveStat[] Moves,MoveStat[] Blocks,BassKindCount[] Kinds,long EmptyBirths,
        (ulong,ulong,ulong,ulong) MoveRng,(ulong,ulong,ulong,ulong) BlockRng,double[] Scales,double Tau);

    internal sealed class Trace
    {
        internal List<Sweep> Sweeps { get; } = [];
        internal void Observe(BassChain chain)
        {
            var moves = chain.Moves.Snapshot(); var blocks = chain.Scan.Snapshot();
            Sweeps.Add(new(BassPriorSamplerTests.Key(chain.Current.Terms),chain.Moves.MoveStats().ToArray(),chain.Scan.BlockStats().ToArray(),
                chain.Structure.Kinds.Entries.ToArray(),chain.Structure.EmptyBirths,
                (moves.Rng0,moves.Rng1,moves.Rng2,moves.Rng3),(blocks.Rng0,blocks.Rng1,blocks.Rng2,blocks.Rng3),
                chain.Current.GroupScales.ToArray(),chain.Current.Tau));
        }
    }

    internal static (int BurnRejected,int SamplingRejected,double Worst) Compare(Trace left,Trace right,int burn,int samples,bool drawn,bool exact)
    {
        Assert.Equal(burn+samples,left.Sweeps.Count); Assert.Equal(left.Sweeps.Count,right.Sweeps.Count);
        int burnRejected = 0, samplingRejected = 0; double worst = 0;
        for (int t = 0; t < left.Sweeps.Count; t++)
        {
            Sweep a = left.Sweeps[t], b = right.Sweeps[t];
            Assert.Equal(a.Structure,b.Structure); Assert.Equal(a.Moves,b.Moves); Assert.Equal(a.Blocks,b.Blocks);
            Assert.Equal(a.Kinds,b.Kinds); Assert.Equal(a.EmptyBirths,b.EmptyBirths);
            Assert.Equal(a.MoveRng,b.MoveRng); Assert.Equal(a.BlockRng,b.BlockRng);
            // This fixture has held τ, Gaussian noise and one structure step plus one MH step per drawn scale.
            // Each block is called exactly once per sweep, so a 0/1 accept-counter increment is the individual decision.
            Assert.Equal(drawn ? 3 : 1,a.Blocks.Length);
            for (int k = 0; k < a.Blocks.Length; k++)
            {
                Assert.Equal(t+1,a.Blocks[k].Attempts);
                long previous = t == 0 ? 0 : left.Sweeps[t-1].Blocks[k].Accepted;
                long decision = a.Blocks[k].Accepted-previous;
                Assert.InRange(decision,0,1);
                if (decision == 0) { if (t < burn) burnRejected++; else samplingRejected++; }
            }
            for (int k = 0; k < a.Scales.Length; k++) worst = Math.Max(worst,Positive(a.Scales[k],b.Scales[k]));
            worst = Math.Max(worst,Positive(a.Tau,b.Tau));
            if (exact) { Bits(a.Scales,b.Scales); Bits([a.Tau],[b.Tau]); }
        }
        Assert.True(burnRejected > 0 && samplingRejected > 0);
        if (drawn)
            for (int k = 1; k < 3; k++)
            {
                Assert.True(left.Sweeps[burn-1].Blocks[k].Accepted < burn); // a rejected scale proposal during burn-in
                Assert.True(left.Sweeps[^1].Blocks[k].Accepted-left.Sweeps[burn-1].Blocks[k].Accepted < samples);
            }
        return (burnRejected,samplingRejected,worst);
    }
}
