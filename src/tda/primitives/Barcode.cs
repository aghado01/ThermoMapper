#nullable enable
using System;
using System.Collections.Generic;
using Maths.Topology;

using TDA.Ph;
namespace TDA.Primitives;

/// <summary>
/// Computes persistence barcodes from a <see cref="NerveFiltration"/> by walking
/// its <see cref="NerveDiff"/> event stream.
/// </summary>
public static class PersistenceBarcode
{
    /// <summary>
    /// Compute the H0 (connected-component) persistence barcode from a
    /// <see cref="NerveFiltration"/>.
    /// </summary>
    public static Barcode ComputeH0(NerveFiltration filtration)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        if (filtration.Frames.Count == 0)
            return new Barcode(Array.Empty<Bar>(), filtration.ParameterLabel);

        var closed = new List<Bar>();
        var active = new Dictionary<int, (double Birth, int Id)>();
        int nextId = 0;

        int nCcs0 = NerveDiff.CountCcs(filtration.Frames[0].Nerve, out _);
        double param0 = filtration.Frames[0].ParameterValue;
        for (int cc = 0; cc < nCcs0; cc++)
            active[cc] = (param0, nextId++);

        foreach (var diff in filtration.ComputeDiffs())
        {
            double paramTo = diff.ParameterTo;
            var nextActive = new Dictionary<int, (double Birth, int Id)>();

            foreach (var evt in diff.ComponentEvents)
            {
                switch (evt.Kind)
                {
                    case ComponentEventKind.Birth:
                        nextActive[evt.CcsTo[0]] = (paramTo, nextId++);
                        break;

                    case ComponentEventKind.Death:
                    {
                        var comp = active[evt.CcsFrom[0]];
                        closed.Add(new Bar(comp.Birth, paramTo, 0));
                        break;
                    }

                    case ComponentEventKind.Continuation:
                        nextActive[evt.CcsTo[0]] = active[evt.CcsFrom[0]];
                        break;

                    case ComponentEventKind.Merge:
                    {
                        (double Birth, int Id) elder = active[evt.CcsFrom[0]];
                        for (int k = 1; k < evt.CcsFrom.Count; k++)
                        {
                            var cand = active[evt.CcsFrom[k]];
                            if (cand.Birth < elder.Birth || (cand.Birth == elder.Birth && cand.Id < elder.Id))
                                elder = cand;
                        }
                        foreach (int fromCc in evt.CcsFrom)
                        {
                            var comp = active[fromCc];
                            if (comp.Id != elder.Id)
                                closed.Add(new Bar(comp.Birth, paramTo, 0));
                        }
                        nextActive[evt.CcsTo[0]] = elder;
                        break;
                    }

                    case ComponentEventKind.Split:
                    {
                        nextActive[evt.CcsTo[0]] = active[evt.CcsFrom[0]];
                        for (int k = 1; k < evt.CcsTo.Count; k++)
                            nextActive[evt.CcsTo[k]] = (paramTo, nextId++);
                        break;
                    }
                }
            }

            active = nextActive;
        }

        foreach (var (_, comp) in active)
            closed.Add(new Bar(comp.Birth, double.PositiveInfinity, 0));

        return new Barcode(closed, filtration.ParameterLabel);
    }
}
