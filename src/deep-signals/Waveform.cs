using System;
using System.IO;
using NAudio.Wave;
using TorchSharp;

namespace DeepSignals;

/// <summary>A segment of a recording at its source sample rate. Dispose it to release the tensor.</summary>
public sealed class Waveform : IDisposable
{
    private Waveform(int sampleRate, torch.Tensor samples)
    {
        SampleRate = sampleRate;
        Samples = samples;
    }

    public int SampleRate { get; }

    /// <summary>CPU float32 samples indexed by [channel, frame].</summary>
    public torch.Tensor Samples { get; }

    /// <summary>Decodes up to <paramref name="maxFrames"/> frames of a seekable PCM or IEEE-float WAV.
    /// Frame indices are relative to the start of the WAV file. The caller owns the stream.
    /// Starting at the end returns an empty frame axis.</summary>
    public static Waveform ReadFrames(Stream source, long startFrame, int maxFrames)
    {
        ArgumentNullException.ThrowIfNull(source);
        ArgumentOutOfRangeException.ThrowIfNegative(startFrame);
        ArgumentOutOfRangeException.ThrowIfNegative(maxFrames);
        if (!source.CanRead || !source.CanSeek)
            throw new ArgumentException("The WAV stream must be readable and seekable.", nameof(source));
        source.Position = 0;
        using var reader = new WaveFileReader(source);
        var provider = reader.ToSampleProvider();
        int channels = provider.WaveFormat.Channels;
        int sampleRate = provider.WaveFormat.SampleRate;
        int blockSize = reader.WaveFormat.BlockAlign;
        if (channels <= 0 || sampleRate <= 0 || blockSize <= 0 || reader.Length % blockSize != 0)
            throw new InvalidDataException("Invalid WAV format or incomplete audio frame.");

        long totalFrames = reader.Length / blockSize;
        ArgumentOutOfRangeException.ThrowIfGreaterThan(startFrame, totalFrames);
        int frames = (int)Math.Min(maxFrames, totalFrames - startFrame);
        if (frames == 0)
            return new Waveform(sampleRate, torch.empty(channels, 0, dtype: torch.ScalarType.Float32, device: torch.CPU));

        reader.Position = startFrame * blockSize;
        var interleaved = new float[checked(frames * channels)];
        int offset = 0;
        while (offset < interleaved.Length)
        {
            int read = provider.Read(interleaved.AsSpan(offset));
            if (read == 0) throw new EndOfStreamException("WAV ended before the requested frames were decoded.");
            offset += read;
        }

        using var flat = torch.tensor(interleaved, device: torch.CPU);
        using var frameMajor = flat.reshape(frames, channels);
        using var channelMajor = frameMajor.transpose(0, 1);
        return new Waveform(sampleRate, channelMajor.contiguous());
    }

    public void Dispose() => Samples.Dispose();
}
