namespace Graphs.Coupling
{
    public readonly record struct MixtureBandwidth(
        double Gaussian,
        double Cauchy,
        double Laplacian);

    public readonly record struct MixtureWeights(
        double Gaussian,
        double Cauchy,
        double Laplacian);
}
