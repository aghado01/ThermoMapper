using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass.Calibration;

namespace LatentRecoveryProbe;

/// <summary>
/// The synthesizer itself behind the emulator contract, at a set of kept frequencies: one draw, no emulator noise. It is
/// the pilot's reference inversion, the same calibration with a perfect generator, so what the BASS emulator adds to the
/// error reads as the difference. Thread-safe: every call allocates its own waveform.
/// </summary>
public sealed class SynthesizerEmulator : ICalibrationEmulator
{
    private readonly Synthesizer _synth;
    private readonly double[,] _cos, _sin;

    public SynthesizerEmulator(Synthesizer synth, IReadOnlyList<CalibrationParameter> parameters, IReadOnlyList<double> keptFrequencies)
    {
        _synth = synth;
        Parameters = parameters;
        (_cos, _sin) = synth.Table(keptFrequencies);
        Outputs = keptFrequencies.Count;
    }

    public IReadOnlyList<CalibrationParameter> Parameters { get; }
    public int DrawCount => 1;
    public IReadOnlyList<double>? DrawWeights => null;
    public int Outputs { get; }
    public EmulatorNoise Noise(int draw) => new IndependentNoise(0.0);
    public double[,]? TruncationResiduals() => null;

    public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output)
    {
        double[] x = _synth.Waveform(continuous);
        Synthesizer.Power(x, _cos, _sin, output);
        for (int l = 0; l < output.Length; l++) output[l] = Math.Log(output[l] + _synth.Floor);
    }
}
