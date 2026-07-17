using System;
using System.Threading;
using System.Threading.Tasks;

namespace Viz.Scene;

/// <summary>
/// Lifecycle boundary implemented by a concrete spatial renderer. Scientific and
/// workbench state are compiled before they reach this interface.
/// </summary>
public interface ISceneRenderer : IAsyncDisposable
{
    ValueTask ApplyAsync(
        SceneSnapshot snapshot,
        CancellationToken cancellationToken = default);
}

/// <summary>Receives renderer interaction as typed view actions.</summary>
public interface IViewActionSink
{
    ValueTask DispatchAsync(
        ViewAction action,
        CancellationToken cancellationToken = default);
}
