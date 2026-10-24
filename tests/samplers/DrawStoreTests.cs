using System;
using System.Linq;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The shared draw store: it keeps every draw by default, one in k when thinned with the first kept, carries the
/// thinning's phase through a copy so a resumed chain keeps the same draws, and concatenates when chains are pooled.
/// </summary>
public sealed class DrawStoreTests
{
    [Fact]
    public void KeepsEveryDraw_ByDefault()
    {
        var store = new DrawStore<int>();
        for (int i = 0; i < 10; i++) store.Add(i);
        Assert.Equal(Enumerable.Range(0, 10), store);
        Assert.Equal(10, store.Count);
        Assert.Equal(10, store.Offered);
    }

    [Fact]
    public void Thinning_KeepsOneInEvery_TheFirstIncluded()
    {
        var store = new DrawStore<int>(every: 3);
        for (int i = 0; i < 10; i++) store.Add(i);
        Assert.Equal(new[] { 0, 3, 6, 9 }, store.ToArray());
        Assert.Equal(10, store.Offered);
    }

    [Fact]
    public void Copy_CarriesThePhase_SoAResumedChainKeepsTheSameDraws()
    {
        var straight = new DrawStore<int>(every: 4);
        for (int i = 0; i < 23; i++) straight.Add(i);

        var first = new DrawStore<int>(every: 4);
        for (int i = 0; i < 10; i++) first.Add(i);
        DrawStore<int> checkpoint = first.Copy();
        first.Add(99);                                   // the live store moves on; the copy does not
        Assert.Equal(10, checkpoint.Offered);

        DrawStore<int> resumed = checkpoint.Copy();
        for (int i = 10; i < 23; i++) resumed.Add(i);
        Assert.Equal(straight.ToArray(), resumed.ToArray());
        Assert.Equal(new[] { 0, 4, 8 }, checkpoint.ToArray());
    }

    [Fact]
    public void AppendFrom_ConcatenatesRecords()
    {
        var a = new DrawStore<string>();
        var b = new DrawStore<string>();
        a.Add("a1"); a.Add("a2");
        b.Add("b1");
        a.AppendFrom(b);
        Assert.Equal(new[] { "a1", "a2", "b1" }, a.ToArray());
        Assert.Equal(3, a.Offered);
        Assert.Throws<ArgumentOutOfRangeException>(() => new DrawStore<int>(every: 0));
    }
}
