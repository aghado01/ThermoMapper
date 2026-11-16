using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace GenerativePrior;

internal sealed class Matrix
{
    public double[][] Rows { get; }
    private readonly double[][] columns;

    public Matrix(int height, int width, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double sd = 1.0 / Math.Sqrt(height);
        Rows = new double[height][];
        columns = new double[width][];
        for (int j = 0; j < width; j++) columns[j] = new double[height];
        for (int i = 0; i < height; i++)
        {
            Rows[i] = new double[width];
            for (int j = 0; j < width; j++)
                columns[j][i] = Rows[i][j] = sd * Normal.SampleStandard(rng);
        }
    }

    public void Apply(double[] x, double[] result)
    {
        for (int i = 0; i < result.Length; i++)
            result[i] = PairwiseSum.Dot(Rows[i], x);
    }

    public void ApplyTranspose(double[] x, double[] result)
    {
        for (int i = 0; i < result.Length; i++)
            result[i] = PairwiseSum.Dot(columns[i], x);
    }
}

// Storage for one of the two signs. Count exact zeros before overwriting preactivations.
internal sealed class Forward(Settings s)
{
    public readonly double[] Hidden = new double[s.HiddenWidth];
    public readonly double[] Output = new double[s.OutputWidth];
    public readonly double[] Residual = new double[s.Measurements];
    public int Zeros;
    public double Objective;
}

internal sealed class Network(Settings s, int[] seeds)
{
    public readonly Matrix W1 = new(s.HiddenWidth, s.CodeDimension, seeds[0]);
    public readonly Matrix W2 = new(s.OutputWidth, s.HiddenWidth, seeds[1]);
    public readonly Matrix A = new(s.Measurements, s.OutputWidth, seeds[2]);
    private readonly double[] hiddenAdjoint = new double[s.HiddenWidth];
    private readonly double[] outputAdjoint = new double[s.OutputWidth];

    private static int Relu(double[] values)
    {
        int zeros = 0;
        for (int i = 0; i < values.Length; i++)
        {
            if (values[i] == 0) zeros++;
            values[i] = Math.Max(0, values[i]);
        }
        return zeros;
    }

    public void Predict(double[] code, Forward state)
    {
        W1.Apply(code, state.Hidden);
        state.Zeros = Relu(state.Hidden);
        W2.Apply(state.Hidden, state.Output);
        state.Zeros += Relu(state.Output);
        A.Apply(state.Output, state.Residual);
    }

    public void Evaluate(double[] code, double[] y, Forward state)
    {
        Predict(code, state);
        for (int i = 0; i < y.Length; i++) state.Residual[i] -= y[i];
        state.Objective = 0.5 * PairwiseSum.Dot(state.Residual, state.Residual);
    }

    public void Gradient(Forward state, double[] result)
    {
        A.ApplyTranspose(state.Residual, outputAdjoint);
        for (int i = 0; i < outputAdjoint.Length; i++)
            if (state.Output[i] == 0) outputAdjoint[i] = 0;
        W2.ApplyTranspose(outputAdjoint, hiddenAdjoint);
        for (int i = 0; i < hiddenAdjoint.Length; i++)
            if (state.Hidden[i] == 0) hiddenAdjoint[i] = 0;
        W1.ApplyTranspose(hiddenAdjoint, result);
    }
}
