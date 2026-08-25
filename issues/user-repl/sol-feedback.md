The module boundary is good; the interior has accumulated coupling.
SpcCommand.cs is roughly 45.6K characters—26% of the entire module—and combines argument parsing, policy validation, graph caching, engine construction, execution, artifact writing, resolver logic, and console presentation.
Argument parsers and CSV/config handling are duplicated across commands.
SpcPreset.ApplyTo(SpcCommand.Options) makes configuration mutate a private CLI-options bag. It is configuration shaped around parser internals rather than a stable declarative request.
SpcUserSession.Run both invokes the engine and writes a fixed artifact bundle. Computation and materialization are therefore difficult to compose separately.
Manifests serialize concrete engine objects and interfaces such as GraphCompilerConfig and IEdgeProjection, coupling the wire schema to implementation types.
SpcUserDataset presents mutable jagged arrays through a nominally immutable record and does not snapshot them.
Graph fingerprints use BitConverter.GetBytes(double) without canonical byte order, framing, or an explicit algorithm/version identity.
Cache loading collapses absent, stale, incompatible, and corrupt artifacts into a single false.
Broad exception catches flatten all failures into exit code 1 and a message.
