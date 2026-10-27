using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The incremental factor against fresh strict factorizations. Births grown from empty reproduce the strict factor
/// bit for bit; randomized birth, death and change sequences track a fresh <see cref="DenseCholesky"/> of the same
/// Gram (factor, log-determinant, diag(A⁻¹)); proposals predict what acceptance does and, rejected, leave the accepted
/// factor's bits alone; leave-one-out pivots match each column's squared sine against the others, computed through
/// the thin SVD; and a death's downdate of diag(A⁻¹) survives the cancellation a near-duplicate column causes.
/// </summary>
public sealed class BorderedCholeskyTests
{
    [Fact]
    public void Births_ReproduceTheStrictFactor_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 701);
        List<double[]> columns = RandomColumns(60, 25, rng);
        double[,] g = Gram(columns);

        var f = new BorderedCholesky(25, rebuildInterval: 0);
        for (int j = 0; j < 25; j++) Birth(f, g, j);
        var strict = new DenseCholesky(25);
        strict.Decompose(g);

        for (int i = 0; i < 25; i++)
            for (int j = 0; j <= i; j++) Assert.Equal(strict.Factor(i, j), f.Factor(i, j));
        Assert.Equal(strict.EquilibratedCondition(), f.EquilibratedCondition());
        Assert.True(Math.Abs(strict.LogDet - f.LogDet) <= 1e-12 * Math.Abs(strict.LogDet));
    }

    [Fact]
    public void RandomUpdates_TrackFreshFactorizations()
    {
        const int rows = 80, capacity = 20, steps = 400;
        var rng = new Xoshiro256PlusPlus(seed: 702);
        var active = new List<double[]>();
        var f = new BorderedCholesky(capacity, rebuildInterval: 0);
        var probe = new double[capacity + 1];

        for (int step = 0; step < steps; step++)
        {
            int n = active.Count;
            int kind = n < 3 ? 0 : n == capacity ? 1 : rng.NextInt(3);   // 0 birth, 1 death, 2 change
            double logDetBefore = f.LogDet;
            var next = new List<double[]>(active);
            FactorProposal p;
            if (kind == 0)
            {
                double[] c = RandomColumns(rows, 1, rng)[0];
                next.Add(c);
                p = f.ProposeBirth(Border(active, c), Dot(c, c));
            }
            else
            {
                int m = rng.NextInt(n);
                next.RemoveAt(m);
                if (kind == 1) p = f.ProposeDeath(m);
                else
                {
                    double[] c = RandomColumns(rows, 1, rng)[0];
                    var remaining = new List<double[]>(next);
                    next.Add(c);
                    p = f.ProposeChange(m, Border(remaining, c), Dot(c, c));
                }
            }
            Assert.True(p.PositiveDefinite);

            for (int i = 0; i < next.Count; i++) probe[i] = rng.NextDouble() - 0.5;
            ReadOnlySpan<double> b = probe.AsSpan(0, next.Count);
            double proposedForm = f.ProposedInverseQuadraticForm(b);
            f.Accept();
            active = next;

            double[,] g = Gram(active);
            var fresh = new DenseCholesky(capacity);
            fresh.Decompose(g);
            Assert.True(fresh.Status.PositiveDefinite);
            int k = active.Count;
            Assert.Equal(k, f.Dimension);

            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < k; j++) Assert.Equal(g[i, j], f.Entry(i, j));
                double scale = Math.Sqrt(g[i, i]);
                for (int j = 0; j <= i; j++)
                    Assert.True(Math.Abs(f.Factor(i, j) - fresh.Factor(i, j)) <= 1e-12 * scale, $"step {step}: L[{i},{j}]");
                var e = new double[k];
                e[i] = 1.0;
                double inverse = fresh.Solve(e)[i];
                Assert.True(Math.Abs(f.InverseDiagonal(i) - inverse) <= 1e-11 * inverse, $"step {step}: (A⁻¹)[{i},{i}]");
            }
            Assert.True(Math.Abs(f.LogDet - fresh.LogDet) <= 1e-11 * Math.Max(1.0, Math.Abs(fresh.LogDet)), $"step {step}: ln|A|");
            Assert.True(Math.Abs(p.LogDetChange - (fresh.LogDet - logDetBefore)) <= 1e-10, $"step {step}: ln|A| change");
            Assert.True(Math.Abs(p.MinLeaveOneOutPivot - f.MinLeaveOneOutPivot()) <= 1e-12 * f.MinLeaveOneOutPivot(), $"step {step}: min leave-one-out");
            double form = fresh.InverseQuadraticForm(b);
            Assert.True(Math.Abs(proposedForm - form) <= 1e-11 * form, $"step {step}: bᵀA⁻¹b {proposedForm:R} vs {form:R}");
        }

        FactorDrift drift = f.Rebuild();
        Assert.True(drift.PositiveDefinite);
        Assert.True(drift.Factor <= 1e-12, $"factor drift {drift.Factor:E3}");
        Assert.True(drift.InverseDiagonal <= 1e-11, $"diag(A⁻¹) drift {drift.InverseDiagonal:E3}");
        var strict = new DenseCholesky(capacity);
        strict.Decompose(Gram(active));
        for (int i = 0; i < active.Count; i++)
            for (int j = 0; j <= i; j++) Assert.Equal(strict.Factor(i, j), f.Factor(i, j));
    }

    [Fact]
    public void RejectedProposals_LeaveTheAcceptedFactorUntouched()
    {
        var rng = new Xoshiro256PlusPlus(seed: 703);
        List<double[]> columns = RandomColumns(30, 8, rng);
        var f = new BorderedCholesky(10, rebuildInterval: 0);
        double[,] g = Gram(columns);
        for (int j = 0; j < 6; j++) Birth(f, g, j);
        Assert.Equal(FactorChange.Death, f.ProposeDeath(2).Change);
        f.Accept();                                         // a factor that has been through a rotation
        List<double[]> active = [columns[0], columns[1], columns[3], columns[4], columns[5]];

        double[] before = Snapshot(f);
        f.ProposeBirth(Border(active, columns[6]), Dot(columns[6], columns[6]));
        f.ProposeDeath(0);
        f.ProposeDeath(4);
        f.ProposeChange(1, Border([columns[0], columns[3], columns[4], columns[5]], columns[7]), Dot(columns[7], columns[7]));
        FactorProposal singular = f.ProposeBirth(Border(active, columns[3]), Dot(columns[3], columns[3]));
        Assert.Equal(before, Snapshot(f));
        Assert.Equal(singular, f.Pending);
        Assert.True(singular.MinLeaveOneOutPivot < 1e-12);

        f.Rebuild();
        Assert.Null(f.Pending);
        Assert.Throws<InvalidOperationException>(f.Accept);
    }

    [Fact]
    public void NonPositiveDefiniteBirth_IsReported_AndRefused()
    {
        var f = new BorderedCholesky(3);
        Assert.True(f.ProposeBirth(ReadOnlySpan<double>.Empty, 1.0).PositiveDefinite);
        f.Accept();
        FactorProposal singular = f.ProposeBirth([1.0], 1.0);                 // pivot exactly 0
        Assert.False(singular.PositiveDefinite);
        Assert.Equal(0.0, singular.RelativePivot);
        Assert.Throws<InvalidOperationException>(f.Accept);
        FactorProposal indefinite = f.ProposeBirth([1.0], 0.5);               // pivot −0.5
        Assert.False(indefinite.PositiveDefinite);
        Assert.True(indefinite.MinLeaveOneOutPivot < 0.0);
        Assert.Throws<InvalidOperationException>(f.Accept);
        Assert.Equal(1, f.Dimension);
    }

    [Fact]
    public void FailedDecompose_StopsAtThePivot_AndClearRestores()
    {
        var f = new BorderedCholesky(3);
        f.Decompose(new double[,] { { 1.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } });
        Assert.Equal(1, f.FailedAt);
        Assert.False(f.Status.PositiveDefinite);
        Assert.Throws<InvalidOperationException>(() => f.ProposeDeath(0));
        f.Clear();
        Assert.Equal(0, f.Dimension);
        Assert.Equal(0.0, f.LogDet);
        Assert.True(f.ProposeBirth(ReadOnlySpan<double>.Empty, 2.0).PositiveDefinite);
    }

    [Fact]
    public void DeathOfANearDuplicate_KeepsItsPartnersInverseDiagonalAccurate()
    {
        // b₂ = b₁ + 10⁻⁶·noise: (A⁻¹)₁₁ is about 10¹², and drops to about 1 when b₂ dies. The downdate alone would keep
        // about four digits; the recomputation it triggers keeps them all.
        var rng = new Xoshiro256PlusPlus(seed: 704);
        List<double[]> columns = RandomColumns(40, 3, rng);
        double[] noise = RandomColumns(40, 1, rng)[0];
        for (int i = 0; i < 40; i++) columns[1][i] = columns[0][i] + 1e-6 * noise[i];
        double[,] g = Gram(columns);

        var f = new BorderedCholesky(3, rebuildInterval: 0);
        for (int j = 0; j < 3; j++) Birth(f, g, j);
        Assert.True(f.InverseDiagonal(0) > 1e10);
        f.ProposeDeath(1);
        f.Accept();

        var fresh = new DenseCholesky(2);
        fresh.Decompose(Gram([columns[0], columns[2]]));
        for (int i = 0; i < 2; i++)
        {
            var e = new double[2];
            e[i] = 1.0;
            double inverse = fresh.Solve(e)[i];
            Assert.True(Math.Abs(f.InverseDiagonal(i) - inverse) <= 1e-13 * inverse, $"(A⁻¹)[{i},{i}] {f.InverseDiagonal(i):R} vs {inverse:R}");
        }
    }

    [Fact]
    public void LeaveOneOutPivot_IsTheSquaredSineAgainstTheOthers()
    {
        // Each column's residual against an orthonormal basis of the others, from the thin SVD.
        const int rows = 50, m = 8;
        var rng = new Xoshiro256PlusPlus(seed: 705);
        List<double[]> columns = RandomColumns(rows, m, rng);
        double[] noise = RandomColumns(rows, 1, rng)[0];
        for (int i = 0; i < rows; i++) columns[5][i] = columns[2][i] + 1e-3 * noise[i];   // s ≈ 10⁻⁶
        double[,] g = Gram(columns);
        var f = new BorderedCholesky(m, rebuildInterval: 0);
        for (int j = 0; j < m; j++) Birth(f, g, j);

        for (int j = 0; j < m; j++)
        {
            var others = new double[rows, m - 1];
            for (int c = 0, t = 0; c < m; c++)
            {
                if (c == j) continue;
                for (int i = 0; i < rows; i++) others[i, t] = columns[c][i];
                t++;
            }
            ThinSvdResult svd = ThinSvd.Decompose(others);
            Assert.Equal(m - 1, svd.NumericalRank);
            double[] residual = (double[])columns[j].Clone();
            foreach (double[] u in svd.LeftVectors)
            {
                double projection = Dot(u, columns[j]);
                for (int i = 0; i < rows; i++) residual[i] -= projection * u[i];
            }
            double sine2 = Dot(residual, residual) / Dot(columns[j], columns[j]);
            double s = f.LeaveOneOutPivot(j);
            Assert.True(Math.Abs(s - sine2) <= 1e-8 * sine2, $"column {j}: {s:R} vs {sine2:R}");
        }
    }

    [Fact]
    public void LeaveOneOutPivots_DoNotDependOnOrder_WhereTheNewestPivotDoes()
    {
        // D06's counterexample: b₁ = e₁, b₂ = e₁ + 10⁻³e₂, b₃ = e₂ + 10⁻³e₃. Added last, b₃'s pivot is 10⁻⁶ and b₁'s
        // 10⁻¹²; the leave-one-out pivots are (10⁻¹², 10⁻¹², 10⁻⁶) whichever is added last.
        double[][] b = [[1.0, 0.0, 0.0], [1.0, 1e-3, 0.0], [0.0, 1.0, 1e-3]];
        var newest = new double[2];
        var sorted = new double[2][];
        int[][] orders = [[0, 1, 2], [1, 2, 0]];
        for (int o = 0; o < 2; o++)
        {
            var columns = new List<double[]>();
            foreach (int c in orders[o]) columns.Add(b[c]);
            var f = new BorderedCholesky(3);
            double[,] g = Gram(columns);
            for (int j = 0; j < 3; j++) Birth(f, g, j);
            newest[o] = f.RelativePivot(2);
            sorted[o] = [f.LeaveOneOutPivot(0), f.LeaveOneOutPivot(1), f.LeaveOneOutPivot(2)];
            Array.Sort(sorted[o]);
        }
        Assert.True(Math.Abs(newest[0] - 1e-6) <= 1e-9, $"order 123: {newest[0]:R}");
        Assert.True(Math.Abs(newest[1] - 1e-12) <= 1e-14, $"order 231: {newest[1]:R}");
        for (int t = 0; t < 3; t++)
            Assert.True(Math.Abs(sorted[0][t] - sorted[1][t]) <= 1e-3 * sorted[0][t], $"leave-one-out {t}: {sorted[0][t]:R} vs {sorted[1][t]:R}");
        Assert.True(Math.Abs(sorted[0][0] - 1e-12) <= 1e-14);
        Assert.True(Math.Abs(sorted[0][2] - 1e-6) <= 1e-9);
    }

    [Fact]
    public void Change_PutsTheNewColumnLast()
    {
        var rng = new Xoshiro256PlusPlus(seed: 706);
        List<double[]> columns = RandomColumns(25, 5, rng);
        var f = new BorderedCholesky(5);
        double[,] g = Gram(columns);
        for (int j = 0; j < 4; j++) Birth(f, g, j);

        List<double[]> remaining = [columns[0], columns[2], columns[3]];
        FactorProposal p = f.ProposeChange(1, Border(remaining, columns[4]), Dot(columns[4], columns[4]));
        Assert.Equal(FactorChange.Change, p.Change);
        Assert.Equal(1, p.Index);
        f.Accept();

        double[,] expected = Gram([columns[0], columns[2], columns[3], columns[4]]);
        var fresh = new DenseCholesky(4);
        fresh.Decompose(expected);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                Assert.Equal(expected[i, j], f.Entry(i, j));
                Assert.True(Math.Abs(f.Factor(i, j) - fresh.Factor(i, j)) <= 1e-13 * Math.Sqrt(expected[i, i]));
            }
    }

    [Fact]
    public void PeriodicRebuild_ResetsTheCount_AndRecordsTheDrift()
    {
        var rng = new Xoshiro256PlusPlus(seed: 707);
        List<double[]> columns = RandomColumns(30, 12, rng);
        double[,] g = Gram(columns);
        var f = new BorderedCholesky(12, rebuildInterval: 5);
        for (int j = 0; j < 4; j++) Birth(f, g, j);
        Assert.Equal(4, f.UpdatesSinceRebuild);
        Assert.Null(f.LastDrift);

        f.ProposeDeath(1);
        f.Accept();                                   // the fifth update triggers a rebuild
        Assert.Equal(0, f.UpdatesSinceRebuild);
        Assert.True(f.LastDrift is { PositiveDefinite: true } drift && drift.Factor <= 1e-13 && drift.InverseDiagonal <= 1e-12);

        f.ProposeDeath(0);
        f.Accept();
        Assert.Equal(1, f.UpdatesSinceRebuild);
    }

    [Fact]
    public void Solves_AndInnovations_HoldAfterRotations()
    {
        var rng = new Xoshiro256PlusPlus(seed: 708);
        List<double[]> columns = RandomColumns(40, 9, rng);
        double[,] g9 = Gram(columns);
        var f = new BorderedCholesky(9, rebuildInterval: 0);
        for (int j = 0; j < 9; j++) Birth(f, g9, j);
        f.ProposeDeath(3);
        f.Accept();
        f.ProposeDeath(0);
        f.Accept();
        List<double[]> active = [columns[1], columns[2], columns[4], columns[5], columns[6], columns[7], columns[8]];
        double[,] a = Gram(active);
        int n = active.Count;

        double[] b = Vector(n, rng);
        (double[] x, double logAbsDet) = LuSolve(a, b);
        AssertClose(x, f.Solve(b), 1e-11);
        Assert.True(Math.Abs(f.LogDet - logAbsDet) <= 1e-11 * Math.Max(1.0, Math.Abs(logAbsDet)));

        var cov = new double[n, n];
        var e = new double[n];
        var v = new double[n];
        for (int k = 0; k < n; k++)
        {
            Array.Clear(e);
            e[k] = 1.0;
            f.SampleInnovation(e, v);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cov[i, j] += v[i] * v[j];
        }
        AssertIdentity(Multiply(a, cov), 1e-11);
        Assert.True(f.EquilibratedCondition() >= 1.0);
    }

    private static void Birth(BorderedCholesky f, double[,] g, int j)
    {
        var border = new double[j];
        for (int k = 0; k < j; k++) border[k] = g[j, k];
        Assert.True(f.ProposeBirth(border, g[j, j]).PositiveDefinite);
        f.Accept();
    }

    private static List<double[]> RandomColumns(int rows, int count, Xoshiro256PlusPlus rng)
    {
        var columns = new List<double[]>(count);
        for (int c = 0; c < count; c++) columns.Add(Vector(rows, rng));
        return columns;
    }

    private static double[,] Gram(List<double[]> columns)
    {
        int k = columns.Count;
        var g = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j <= i; j++) g[i, j] = g[j, i] = Dot(columns[i], columns[j]);
        return g;
    }

    private static double[] Border(List<double[]> columns, double[] c)
    {
        var u = new double[columns.Count];
        for (int k = 0; k < columns.Count; k++) u[k] = Dot(columns[k], c);
        return u;
    }

    private static double Dot(double[] a, double[] b) => PairwiseSum.Dot(a, b);

    private static double[] Snapshot(BorderedCholesky f)
    {
        int n = f.Dimension;
        var s = new List<double> { f.LogDet, f.Status.MinPivot, f.Status.MaxPivot, f.Status.MinRelativePivot };
        for (int i = 0; i < n; i++)
        {
            s.Add(f.InverseDiagonal(i));
            for (int j = 0; j < n; j++)
            {
                s.Add(f.Factor(i, j));
                s.Add(f.Entry(i, j));
            }
        }
        return [.. s];
    }
}
