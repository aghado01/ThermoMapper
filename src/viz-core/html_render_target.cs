using System;
using System.IO;
using System.Reflection;
using System.Text;
using Viz;

namespace Viz.Renderers
{
    /// <summary>
    /// Renders a ScenePackage as a self-contained Three.js HTML file.
    ///
    /// All serialization is delegated to <see cref="JsonExportRenderTarget"/>.
    /// The compact JSON payload is injected into the static viewer template
    /// (viewer.html, embedded in this assembly) at a single __SCENE_DATA__
    /// placeholder. The viewer template owns all HTML, CSS, and JS — this
    /// class only performs the injection.
    /// </summary>
    public sealed class ThreeJsHtmlRenderTarget : IRenderTarget
    {
        private static readonly JsonExportRenderTarget _json = new(compact: true);
        private static readonly string _template = LoadTemplate();

        private static string LoadTemplate()
        {
            var asm = Assembly.GetExecutingAssembly();
            using var stream = asm.GetManifestResourceStream("Viz.viewer.html")
                ?? throw new InvalidOperationException(
                    "Embedded resource 'Viz.viewer.html' not found. " +
                    "Ensure viewer.html is marked EmbeddedResource in VizCore.csproj.");
            using var reader = new StreamReader(stream, Encoding.UTF8);
            return reader.ReadToEnd();
        }

        public void Render(ScenePackage scene, Stream output)
        {
            if (scene is null) throw new ArgumentNullException(nameof(scene));
            if (output is null) throw new ArgumentNullException(nameof(output));

            // Serialize via the canonical JSON renderer (compact — no whitespace bloat)
            using var jsonStream = new MemoryStream();
            _json.Render(scene, jsonStream);
            string json = Encoding.UTF8.GetString(jsonStream.ToArray());

            // Inject at the single defined placeholder in the viewer template
            string html = _template.Replace("__SCENE_DATA__", json);
            byte[] bytes = Encoding.UTF8.GetBytes(html);
            output.Write(bytes, 0, bytes.Length);
        }
    }
}
