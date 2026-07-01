using System;

namespace Maths.Geometry.Bandwidth
{
    /// <summary>
    /// Intrinsic Gaussian heat-kernel bandwidth calibration on H^d (Poincaré ball).
    /// Moment-matches observed k-NN distances against the Van Vleck–corrected radial reference.
    /// </summary>
    public static class Hyperbolic
    {
        private const double ReferenceLocalPopulation = 32.0;
        private const double LogScaleFloor = 1e-12;

        public static double MatchBeta(double observedSecondMoment, int ambientDimension)
        {
            if (!double.IsFinite(observedSecondMoment) || observedSecondMoment <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(observedSecondMoment), observedSecondMoment, "Observed second moment must be finite and positive.");

            double guess = Math.Max(1e-6, ambientDimension / Math.Max(2.0 * observedSecondMoment, 1e-6));
            double betaLow = guess / 16.0;
            double betaHigh = guess * 16.0;

            double momentLow = SecondMoment(ambientDimension, betaLow);
            int expansions = 0;
            while (momentLow < observedSecondMoment && expansions++ < 16)
            {
                betaLow *= 0.25;
                momentLow = SecondMoment(ambientDimension, betaLow);
            }

            double momentHigh = SecondMoment(ambientDimension, betaHigh);
            expansions = 0;
            while (momentHigh > observedSecondMoment && expansions++ < 16)
            {
                betaHigh *= 4.0;
                momentHigh = SecondMoment(ambientDimension, betaHigh);
            }

            for (int iteration = 0; iteration < 48; iteration++)
            {
                double betaMid = 0.5 * (betaLow + betaHigh);
                double momentMid = SecondMoment(ambientDimension, betaMid);
                if (momentMid > observedSecondMoment)
                    betaLow = betaMid;
                else
                    betaHigh = betaMid;
            }

            return 0.5 * (betaLow + betaHigh);
        }

        public static double SecondMoment(int ambientDimension, double beta)
        {
            if (ambientDimension < 1)
                throw new ArgumentOutOfRangeException(nameof(ambientDimension), ambientDimension, "Ambient dimension must be positive.");
            if (!double.IsFinite(beta) || beta <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(beta), beta, "Beta must be finite and positive.");

            const int gridCount = 1024;
            double radiusMax = ResolveRadiusMax(ambientDimension, beta);
            double maxLogWeight = double.NegativeInfinity;

            for (int index = 1; index <= gridCount; index++)
            {
                double radius = radiusMax * index / gridCount;
                double logWeight = LogGaussianWeight(ambientDimension, beta, radius);
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
                double weight = Math.Exp(LogGaussianWeight(ambientDimension, beta, radius) - maxLogWeight);
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

        private static double LogGaussianWeight(int ambientDimension, double beta, double radius)
        {
            if (radius <= 0.0)
                return double.NegativeInfinity;

            double logWeight = -beta * radius * radius;
            double exponent = (ambientDimension - 1) / 2.0;
            if (exponent <= 0.0)
                return logWeight;

            return logWeight + exponent * (Math.Log(radius) + LogSinh(radius));
        }

        private static double LogSinh(double radius)
        {
            if (radius < 1e-8)
                return Math.Log(Math.Max(radius, LogScaleFloor));

            if (radius < 20.0)
                return Math.Log(Math.Sinh(radius));

            return radius - Math.Log(2.0);
        }

        private static double ResolveRadiusMax(int ambientDimension, double beta)
        {
            double bandwidth = HeatKernel.BandwidthFromBeta(beta);
            double bulkMode = ambientDimension <= 1
                ? 0.0
                : ((ambientDimension - 1) * bandwidth * bandwidth) / 2.0;

            return Math.Max(8.0, bulkMode + 12.0 * bandwidth + ambientDimension);
        }

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
