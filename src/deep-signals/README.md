# deep-signals

This .NET 10 library provides audio input for models trained with TorchSharp. `Waveform.ReadFrames` decodes a bounded window of a PCM or IEEE-float WAV with NAudio.Core. It preserves the sample rate and channel order and returns a CPU float32 tensor shaped `[channels, frames]`. A final window can be shorter than requested.

```csharp
using DeepSignals;
using System.IO;

using var recording = File.OpenRead(path);
using var clip = Waveform.ReadFrames(recording, startFrame: 0, maxFrames: 16_000);
using var batch = clip.Samples.unsqueeze(0); // [1, channels, frames]
```

The host supplies a native LibTorch runtime compatible with the pinned TorchSharp package. The test project selects LibTorch 2.10.0 CPU for the SDK host and checks WAV decoding, the gradient, and one optimizer update on a known PCM waveform. Run it from the repository root.

```powershell
dotnet test tests/deep-signals.Tests/deep-signals.Tests.csproj -c Release
```
