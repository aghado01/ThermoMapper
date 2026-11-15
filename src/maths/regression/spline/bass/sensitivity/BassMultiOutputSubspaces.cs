using System;
using System.Linq;
using Maths.Regression.Subspaces;

namespace Maths.Regression.Spline.Bass.Sensitivity;

internal static partial class BassActiveSubspaces
{
    internal static CoActiveAnalysis MultiOutput(BassFit fit, BassData data, CoActiveSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        if (spec.OtherDraws is not null) throw new ArgumentException("Two outputs of one fit pair each draw with itself.", nameof(spec));
        double[] a = BassMultiOutput.Functional(data.Response.Outputs, spec.Output ?? throw new ArgumentException("The first output functional is required.", nameof(spec)));
        double[] b = BassMultiOutput.Functional(data.Response.Outputs, spec.OtherOutput ?? throw new ArgumentException("The second output functional is required.", nameof(spec)));
        var inputs = BassInputs.Prepare(data, fit.Transform);
        Frame frame = FrameOf(inputs, fit.Transform, fit.Config.Structure.Degree, fit.Config.Structure.ResolvedHingeNormalization, spec.Measures, spec.Sample, spec.Metric, spec.Route, spec.Points, 1,
                              Separable(inputs, fit.Config.Structure.Degree, fit.Config.Structure.ResolvedHingeNormalization, [fit.Terms]));
        var pooled = fit.Draws.SelectMany(c => c).ToArray();
        var (selected, weights) = Draws3.Select(spec.Draws, pooled.Length, fit.DrawWeights, nameof(spec));
        var pairs = new CoActivePair[selected.Length]; CoActiveError[]? errors = null;
        if (frame.Route == ActiveSubspaceRoute.Exact)
        {
            var moments = BassMoments.On(inputs, frame.Degree, frame.Measure!, frame.Hinges);
            for (int s = 0; s < selected.Length; s++)
            {
                var draw = pooled[selected[s]];
                var terms = draw.Structure.TermIds.Select(id => fit.Terms[id]).ToArray();
                pairs[s] = Pair(moments, frame, terms, BassMultiOutput.Project(draw, a)[1..], terms, BassMultiOutput.Project(draw, b)[1..]);
            }
        }
        else
        {
            var points = PointsOf(frame, spec.Points, spec.Seed); errors = new CoActiveError[selected.Length];
            // The gradient in draw batches of about BatchValues partials, as E8's functional readouts batch theirs; each draw's pair is the same bits.
            int batch = (int)Math.Max(1, BatchValues / ((long)points.N * data.Response.Outputs * frame.P));
            for (int start = 0; start < selected.Length; start += batch)
            {
                int end = Math.Min(selected.Length, start + batch);
                var gradient = BassMultiOutput.Gradient(fit, data, points.Cases, selected[start..end]);
                for (int s = 0; s < end - start; s++)
                    (pairs[start + s], errors[start + s]) = SampledPair(frame, points, Directional(gradient, s, a), Directional(gradient, s, b));
            }
        }
        return new CoActiveAnalysis(frame.Names, frame.Excluded, frame.Route, spec.Metric, frame.Map, frame.Measure, PointCount(frame, spec.Points),
            CoActivePosterior.Of(pairs, weights), errors, selected, null);
    }
}
