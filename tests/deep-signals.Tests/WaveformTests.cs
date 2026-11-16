using System;
using System.IO;
using NAudio.Wave;
using TorchSharp;
using Xunit;

namespace DeepSignals.Tests;

public class WaveformTests
{
    [Fact]
    public void StereoWindowPreservesRateChannelsAndFinalFrame()
    {
        using var bytes = new MemoryStream();
        using (var writer = new WaveFileWriter(bytes, WaveFormat.CreateIeeeFloatWaveFormat(48000, 2)))
            writer.WriteSamples([0f, 0f, -0.75f, 0.25f, 0.5f, -0.125f], 0, 6);
        using var source = new MemoryStream(bytes.ToArray());
        using (var window = Waveform.ReadFrames(source, 1, 50))
        {
            Assert.Equal(48000, window.SampleRate);
            Assert.Equal(new long[] { 2, 2 }, window.Samples.shape);
            Assert.Equal(torch.ScalarType.Float32, window.Samples.dtype);
            using var values = window.Samples.data<float>();
            Assert.Equal(new[] { -0.75f, 0.5f, 0.25f, -0.125f }, values.ToArray());
        }
        Assert.True(source.CanRead);
        using var end = Waveform.ReadFrames(source, 3, 10);
        Assert.Equal(new long[] { 2, 0 }, end.Samples.shape);
        Assert.Throws<ArgumentOutOfRangeException>(() => Waveform.ReadFrames(source, 4, 1));
    }

    [Fact]
    public void MonoPcmFeedsAKnownGradientAndOptimizerStep()
    {
        using var bytes = new MemoryStream();
        using (var writer = new WaveFileWriter(bytes, new WaveFormat(16000, 16, 1)))
            writer.WriteSamples(new short[] { 16384, -16384 }, 0, 2);
        using var source = new MemoryStream(bytes.ToArray());
        using var clip = Waveform.ReadFrames(source, 0, 2);
        Assert.Equal(16000, clip.SampleRate);
        Assert.Equal(new long[] { 1, 2 }, clip.Samples.shape);
        using var scope = torch.NewDisposeScope();
        using var x = clip.Samples.transpose(0, 1);
        using var layer = torch.nn.Linear(1, 1, hasBias: false);
        var weight = layer.weight!;
        using (torch.no_grad()) weight.zero_();
        using var optimizer = torch.optim.SGD(layer.parameters(), 1);
        using var predicted = layer.forward(x);
        using var loss = torch.nn.functional.mse_loss(predicted, x);
        loss.backward();
        using var gradient = weight.grad!;
        Assert.Equal(-0.5f, gradient.item<float>());
        optimizer.step();
        Assert.Equal(0.5f, weight.item<float>());
    }
}
