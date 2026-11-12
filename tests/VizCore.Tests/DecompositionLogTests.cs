using System.Threading.Tasks;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The decomposition log records each kernel's decompositions inside its scope, and nothing outside one: a scope belongs to the
/// async flow that opened it, nested scopes each record what runs inside them, and work on tasks the flow starts is recorded too.
/// </summary>
public sealed class DecompositionLogTests
{
    [Fact]
    public async Task RecordsEachKernelInsideItsScope()
    {
        var symmetric = new double[,] { { 2.0, 1.0 }, { 1.0, 3.0 } };
        var tall = new double[,] { { 1.0, 2.0 }, { 3.0, 4.0 }, { 5.0, 6.0 } };
        EigenFast.DecomposeSymmetric(symmetric);   // before any scope: not recorded

        using DecompositionLog outer = DecompositionLog.Begin();
        EigenFast.DecomposeSymmetric(symmetric);
        DenseEigen.DecomposeSymmetric(symmetric);   // the dispatcher records once, at the kernel it runs
        using (DecompositionLog inner = DecompositionLog.Begin())
        {
            Eigen.DecomposeSymmetric(symmetric);
            await Task.Run(() => ThinSvd.Decompose(tall));
            Assert.Equal([new Decomposition(DecompositionKernel.ReferenceJacobi, 2, 2), new Decomposition(DecompositionKernel.ThinSvd, 3, 2)],
                         inner.Entries);
        }
        EigenFast.DecomposeSymmetric(symmetric);

#if EIGEN_REFERENCE
        const DecompositionKernel dense = DecompositionKernel.ReferenceJacobi;
#else
        const DecompositionKernel dense = DecompositionKernel.FastJacobi;
#endif
        Assert.Equal(5, outer.Entries.Count);
        Assert.Equal(DecompositionKernel.FastJacobi, outer.Entries[0].Kernel);
        Assert.Equal(dense, outer.Entries[1].Kernel);
        Assert.Equal(1, outer.Count(DecompositionKernel.ThinSvd));
        Assert.Equal(DecompositionKernel.FastJacobi, outer.Entries[4].Kernel);
    }
}
