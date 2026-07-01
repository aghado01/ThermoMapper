#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

public sealed class ZigzagBarcodeNaiveTests
{
    private static void AssertBar(Barcode barcode, int p, double birth, double death, IntervalEnd bEnd, IntervalEnd dEnd)
    {
        var bars = barcode.Bars.Where(b => b.Dimension == p).ToList();
        Assert.Contains(bars, b => 
            b.Birth == birth && 
            b.Death == death && 
            b.BirthEnd == bEnd && 
            b.DeathEnd == dEnd);
    }

    [Fact]
    public void AddThenDeleteOneCell_ProducesClosedOpenInterval()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); // Step 0: Add vertex 0
        f.Delete(0);          // Step 1: Delete vertex 0

        var bc = ZigzagBarcodeNaive.Compute(f, 0);

        // Birth at step 0 (Add -> Closed)
        // Death at step 1 (Delete -> Closed)
        AssertBar(bc, 0, 0, 1, IntervalEnd.Closed, IntervalEnd.Closed);
    }

    [Fact]
    public void TriangleFormsThenBreaks_ProducesClosedClosedDeath_H1()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); // v0
        f.Add(1, new int[0]); // v1
        f.Add(2, new int[0]); // v2
        f.Add(3, new[] { 0, 1 }); // e01
        f.Add(4, new[] { 1, 2 }); // e12
        f.Add(5, new[] { 0, 2 }); // e02 (loop forms at step 5!)
        f.Delete(5); // e02 deleted (loop breaks at step 6!)

        var bc = ZigzagBarcodeNaive.Compute(f, 1);

        // H1 interval born at 5 (Add -> Closed), dies at 6 (Delete -> Closed)
        AssertBar(bc, 1, 5, 6, IntervalEnd.Closed, IntervalEnd.Closed);
    }

    [Fact]
    public void DynamicGraphH0_ComponentsMergeThenSplit()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); // step 0: v0
        f.Add(1, new int[0]); // step 1: v1
        f.Add(2, new int[0]); // step 2: v2
        f.Add(3, new[] { 0, 1 }); // step 3: merge v0, v1 (v1 dies, v0 lives)
        f.Add(4, new[] { 1, 2 }); // step 4: merge v0+1, v2 (v2 dies)
        f.Delete(3); // step 5: split v0 and v1+2 (v1 is born again!)
        
        var bc = ZigzagBarcodeNaive.Compute(f, 0);

        // v0 born at 0, never dies? Wait, we don't delete v0, v1, v2.
        // At the end (step 6), we have v0, and edge (1,2) connecting v1 and v2.
        // Two components at the end: v0, and v1+v2.
        
        // Expected H0 bars:
        // v0: [0, 6] (assuming it never dies, it goes to m=6. Wait, death is 6 with Open end if it survives to the end?
        // Let's just check the counts of intervals.
        Assert.Equal(4, bc.Bars.Count(b => b.Dimension == 0));
        
        // v2 dies at step 4 (Add e12). Birth 2 (Closed), Death 4 (Open since it dies by Add).
        AssertBar(bc, 0, 2, 4, IntervalEnd.Closed, IntervalEnd.Open);
        
        // v1 born at 1 (Closed), dies at 3 (Open, via Add e01).
        AssertBar(bc, 0, 1, 3, IntervalEnd.Closed, IntervalEnd.Open);
        
        // v1 is reborn at step 5 via Delete e01 -> Open birth (born by a backward/Delete arrow).
        // Survives to the end (step 6) -> Open death.
        AssertBar(bc, 0, 5, 6, IntervalEnd.Open, IntervalEnd.Open);
        
        // v0 born at 0 (Closed), lives to end (6) -> Open.
        AssertBar(bc, 0, 0, 6, IntervalEnd.Closed, IntervalEnd.Open);
    }

    [Fact]
    public void MinimalFourTypeTaxonomy_CC_CO_OC_OO()
    {
        // 1. Closed-Closed (cc): Born by Add, Dies by Delete
        var f_cc = new ZigzagFiltration();
        f_cc.Add(0, new int[0]);
        f_cc.Delete(0);
        var bc_cc = ZigzagBarcodeNaive.Compute(f_cc, 0);
        AssertBar(bc_cc, 0, 0, 1, IntervalEnd.Closed, IntervalEnd.Closed);

        // 2. Closed-Open (co): Born by Add, Dies by Add
        var f_co = new ZigzagFiltration();
        f_co.Add(0, new int[0]);
        f_co.Add(1, new int[0]);
        f_co.Add(2, new[] { 0, 1 }); // v1 dies by Add
        var bc_co = ZigzagBarcodeNaive.Compute(f_co, 0);
        AssertBar(bc_co, 0, 1, 2, IntervalEnd.Closed, IntervalEnd.Open);

        // 3. Open-Closed (oc): Born by Delete, Dies by Delete
        var f_oc = new ZigzagFiltration();
        f_oc.Add(0, new int[0]);
        f_oc.Add(1, new int[0]);
        f_oc.Add(2, new[] { 0, 1 });
        f_oc.Delete(2); // v1 born by Delete
        f_oc.Delete(1); // v1 dies by Delete
        var bc_oc = ZigzagBarcodeNaive.Compute(f_oc, 0);
        // Step 0: v0
        // Step 1: v1
        // Step 2: e01 (v1 dies at 2: [1, 2, Closed, Open])
        // Step 3: del e01 (v1 born at 3: Open birth since K_2 -> K_3 is Delete)
        // Step 4: del v1 (v1 dies at 4: Closed death since K_4 <- K_5 is Delete... wait, 4 is a Delete of v1)
        AssertBar(bc_oc, 0, 3, 4, IntervalEnd.Open, IntervalEnd.Closed);

        // 4. Open-Open (oo): Born by Delete, Dies by Add
        var f_oo = new ZigzagFiltration();
        f_oo.Add(0, new int[0]);
        f_oo.Add(1, new int[0]);
        f_oo.Add(2, new int[0]);
        f_oo.Add(3, new[] { 0, 1 });
        f_oo.Add(4, new[] { 1, 2 }); // v1, v2 merged into v0
        f_oo.Delete(3); // e01 deleted. v1 is born again by Delete.
        f_oo.Add(5, new[] { 0, 1 }); // e01 added back. v1 dies by Add.
        var bc_oo = ZigzagBarcodeNaive.Compute(f_oo, 0);
        // Step 5: del e01 -> born by Delete -> Open.
        // Step 6: add e01 -> dies by Add -> Open.
        AssertBar(bc_oo, 0, 5, 6, IntervalEnd.Open, IntervalEnd.Open);
    }
}
