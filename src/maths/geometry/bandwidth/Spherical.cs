using System;

namespace Maths.Geometry.Bandwidth
{
    /// <summary>
    /// Intrinsic Gaussian heat-kernel bandwidth calibration on S^m.
    /// <para>
    /// <paramref name="intrinsicDimension"/> is m (manifold dimension), not ambient n.
    /// For points in R^n on S^(n-1), pass m = n − 1.
    /// </para>
    /// </summary>
    public static class Spherical
    {
        private const double ReferenceLocalPopulation = 32.0;
        private const double LogScaleFloor = 1e-12;
        private const double BetaFloor = 1e-6;

        public static double MatchBeta(double observedSecondMoment, int intrinsicDimension)
        {
            if (!double.IsFinite(observedSecondMoment) || observedSecondMoment <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(observedSecondMoment), observedSecondMoment, "Observed second moment must be finite and positive.");

            double ceilingMoment = SecondMoment(intrinsicDimension, BetaFloor);
            if (observedSecondMoment >= ceilingMoment)
                return BetaFloor;

            double guess = Math.Max(BetaFloor, intrinsicDimension / Math.Max(2.0 * observedSecondMoment, 1e-6));
            double betaLow = guess / 16.0;
            double betaHigh = guess * 16.0;

            double momentLow = SecondMoment(intrinsicDimension, betaLow);
            int expansions = 0;
            while (momentLow < observedSecondMoment && betaLow > BetaFloor && expansions++ < 32)
            {
                betaLow = Math.Max(BetaFloor, betaLow * 0.25);
                momentLow = SecondMoment(intrinsicDimension, betaLow);
            }

            double momentHigh = SecondMoment(intrinsicDimension, betaHigh);
            expansions = 0;
            while (momentHigh > observedSecondMoment && expansions++ < 16)
            {
                betaHigh *= 4.0;
                momentHigh = SecondMoment(intrinsicDimension, betaHigh);
            }

            for (int iteration = 0; iteration < 48; iteration++)
            {
                double betaMid = 0.5 * (betaLow + betaHigh);
                double momentMid = SecondMoment(intrinsicDimension, betaMid);
                if (momentMid > observedSecondMoment)
                    betaLow = betaMid;
                else
                    betaHigh = betaMid;
            }

            return 0.5 * (betaLow + betaHigh);
        }

        public static double SecondMoment(int intrinsicDimension, double beta)
        {
            if (intrinsicDimension < 1)
                throw new ArgumentOutOfRangeException(nameof(intrinsicDimension), intrinsicDimension, "Intrinsic dimension must be positive.");
            if (!double.IsFinite(beta) || beta <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(beta), beta, "Beta must be finite and positive.");

            const int gridCount = 1024;
            double radiusMax = ResolveRadiusMax();
            double maxLogWeight = double.NegativeInfinity;

            for (int index = 1; index <= gridCount; index++)
            {
                double radius = radiusMax * index / gridCount;
                double logWeight = LogGaussianWeight(intrinsicDimension, beta, radius);
                if (logWeight > maxLogWeight)
                    maxLogWeight = logWeight;
            }

            var radii = new double[gridCount + 1];
            var cumulativeMass = new double[gridCount + 1];
            var weights = new double[gridCount + 1];
            double previousRadius = 0.0;
            double previousWeight = 0.0;

            for (int index = 1; index <= gridCount; index++)
            {
                double radius = radiusMax * index / gridCount;
                double weight = Math.Exp(LogGaussianWeight(intrinsicDimension, beta, radius) - maxLogWeight);
                double step = radius - previousRadius;

                radii[index] = radius;
                weights[index] = weight;
                cumulativeMass[index] = cumulativeMass[index - 1] + 0.5 * (weight + previousWeight) * step;

                previousRadius = radius;
                previousWeight = weight;
            }

            double totalMass = cumulativeMass[gridCount];
            double normalization = 1.0 - Math.Exp(-ReferenceLocalPopulation);
            double weightedMoment = 0.0;

            for (int index = 1; index <= gridCount; index++)
            {
                double leftRadius = radii[index - 1];
                double rightRadius = radii[index];
                double step = rightRadius - leftRadius;
                double leftPdf = NearestNeighborPdf(weights[index - 1], cumulativeMass[index - 1], totalMass, normalization);
                double rightPdf = NearestNeighborPdf(weights[index], cumulativeMass[index], totalMass, normalization);
                weightedMoment += 0.5 * ((leftRadius * leftRadius * leftPdf) + (rightRadius * rightRadius * rightPdf)) * step;
            }

            return weightedMoment;
        }

        private static double LogGaussianWeight(int intrinsicDimension, double beta, double radius)
        {
            if (radius <= 0.0 || radius >= Math.PI)
                return double.NegativeInfinity;

            double logWeight = -beta * radius * radius;
            double exponent = (intrinsicDimension - 1) / 2.0;
            if (exponent <= 0.0)
                return logWeight;

            return logWeight + exponent * (Math.Log(radius) + LogSin(radius));
        }

        private static double LogSin(double radius)
        {
            if (radius < 1e-8)
                return Math.Log(Math.Max(radius, LogScaleFloor));

            return Math.Log(Math.Max(Math.Sin(radius), LogScaleFloor));
        }

        private static double ResolveRadiusMax() => Math.PI - 1e-9;

        private static double NearestNeighborPdf(double weight, double cumulativeMass, double totalMass, double normalization)
        {
            if (weight <= 0.0 || totalMass <= 0.0)
                return 0.0;

            double intensity = ReferenceLocalPopulation * cumulativeMass / totalMass;
            return (ReferenceLocalPopulation / totalMass)
                * weight
                * Math.Exp(-intensity)
                / Math.Max(normalization, 1e-12);
        }
    }
}
