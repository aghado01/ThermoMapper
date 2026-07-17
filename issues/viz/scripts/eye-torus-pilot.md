# EyeTorus pilot — The crime of proximity

## Header

- **Film ID:** `thermomapper.eye.proximity.v1`
- **Status:** treatment; first v2 vertical slice
- **Purpose:** demonstrative, diagnostic, exploratory
- **Question:** when does ambient proximity misrepresent the intrinsic eye, and
  how do graph construction, coupling, and geometry change the answer?
- **Audience:** technically curious; no SPC knowledge required for the first cut
- **Runtime modes:** live execution and static replay
- **Canonical fixture:** `EyeTorusToy`, three-dimensional Euclidean realization,
  seed 42
- **Related takes:** Poincaré realization and four-dimensional unlinked
  realization

The pilot ends at graph diagnosis. Empirical flow and SPC are sequel sequences
which reuse the study once their adapters are available.

## Scientific contract

The canonical graph producers receive only observed sample features and a pinned
graph recipe. They cannot consume the skeleton, generative labels, latent stroke
parameters, or cross-realization correspondence.

The study may additionally contain oracle artifacts:

- coarse hierarchy: structure versus background;
- fine hierarchy: iris, upper arc, lower arc, optional pupil, background;
- clean skeleton strokes and their parameterization;
- generating tangent/frame information when the Synthetic producer exposes it;
- exact correspondence between structural samples shared by Euclidean and
  Poincaré takes.

Oracle cross-group edge masks score an observed graph. They never add, remove, or
reweight its edges.

## Takes

### Take A — Euclidean eye

- `EyeTorusToy.EyeTorusToyConfig`, dimension 3, seed 42.
- Canonical release take uses the generator's pinned production configuration.
- Development may use reduced point counts, but the UI must disclose that it is a
  development take and the package must preserve the exact config.
- Feature and initial display coordinates happen to agree, but remain separate
  coordinate artifacts.

### Take B — graph alternative

- Reuses Take A's exact sample entity set and feature coordinates.
- Changes only the declared graph recipe dimension under investigation: graph
  rule, metric, connectivity policy, kernel, or bandwidth.
- The first implementation should compare one adjacency-changing alternative;
  later scripts can freeze adjacency and isolate coupling.

### Take C — hyperbolic eye

- `HyperbolicEyeTorus` with the same dimension and seed.
- Structural samples correspond to Take A through their shared local sample;
  the correspondence must be explicit rather than inferred from array order.
- Background samples do not currently correspond: Euclidean background is a
  colored spectral field and hyperbolic background is sampled from native ball
  volume. They dissolve/reappear during a morph instead of being paired falsely.
- Poincaré display uses a visible ball boundary and geodesic graph arcs.

### Take D — unlinked 4D eye

- Dimension 4 with stroke planes fanned according to the fixture.
- Requires an explicit 4D-to-3D projection artifact and projection controls.
- No claim of matched pointwise noise with Take A is made unless the Synthetic
  case contract later guarantees it.
- Hyperbolic 4D background is blocked as evidence until its radial volume density
  uses the dimension-correct `sinh^(D-1)(rho)` factor.

## Required artifacts

- sample entity set with stable IDs and source-role relation;
- feature coordinate set;
- one or more display coordinate sets with projection provenance;
- fixture and realizer manifests;
- oracle label fields at both hierarchy levels;
- optional oracle skeleton and tangent artifacts;
- two complete graph build artifacts, preserving CSR row offsets, columns,
  weights, edge slots, weight kind, metric properties, repair provenance, health
  report, and diagnostics log;
- diagnostic edge fields for common/only-A/only-B and oracle cross-group status;
- graph comparison report;
- two spatial panel descriptors and one report/chart region;
- linked selection and linked camera policy;
- film shots and cues.

## Panel plan

The opening uses one spatial panel. Graph comparison expands into two linked
spatial panels with a shared inspector. Diagnostics occupy a collapsible lower or
right panel rather than being squeezed into 3D labels.

The inspector presents, as applicable:

- sample ID, source role, coordinates, and labels when oracle disclosure permits;
- edge CSR slot, endpoints, metric distance, weight, construction source, repair
  status, and diagnostic classifications;
- producer, inputs, requested/resolved configuration, fingerprint, and warnings.

## Sequence

### Shot 01 — The plate

**Beat:** recognition.

Begin with a face-on orthographic rendering of the clean, flattened generating
strokes: iris, brow, and lower arc. This is explicitly labeled **oracle anatomy**.
The visual quotes the familiar two-dimensional eye without yet claiming that the
method observes it.

Entry cue: cinematic fade. No scientific state changes.

Allowed interaction: pause, orbit disabled until the lift begins, inspect the
fixture reference and take provenance.

### Shot 02 — The lift

**Beat:** the plate becomes a spatial object.

Interpolate from the flattened presentation coordinates to the Euclidean 3D
realization while the camera moves to a shallow oblique angle. Reveal tube
cross-sections, taper, density bias, and optional pupil in separate cues.

Truth note: interpolation between coordinate artifacts is
`PresentationDerived`; the endpoint coordinates come from the fixture producer.

### Shot 03 — How the world was made

**Beat:** distinguish latent cause from observation.

Accumulate structural samples on the skeleton, introduce jitter, then add the
colored background blanket. Construction layers use oracle styling and remain
clearly labeled.

Truth note: the accumulation is a cinematic reveal unless the Synthetic producer
later emits an actual sampling trace.

### Shot 04 — Epistemic cut

**Beat:** everything privileged disappears.

Fade skeleton, source colors, labels, latent parameters, and clean positions.
Leave one neutrally styled observed point cloud. The inspector switches to its
canonical-method admissible view.

Narration beat: “Everything before this moment is how the world was made.
Everything after it is what the method is allowed to know.”

The resulting shot is the default workbench entry state for users who skip the
opening.

### Shot 05 — One neighborhood

**Beat:** a neighborhood is a hypothesis.

Select a sample in a geometrically ambiguous region. Show metric distances to a
small candidate set, then reveal which candidates the pinned graph rule selects.
The inspector displays real producer values.

No full-graph oracle classification is visible yet.

### Shot 06 — The first graph

**Beat:** local choices become global structure.

Reveal graph A progressively for legibility, then settle on a level-of-detail
representation. Display graph health and construction provenance outside the
spatial view.

Truth note: progressive edge drawing is cinematic. The graph was produced as one
artifact. Edge thinning is `PresentationDerived` and the shown/full counts remain
visible.

### Shot 07 — The crime

**Beat:** apparent closeness can cross intrinsic anatomy.

Enable the oracle fine-level cross-group edge field. Cross-group edges become a
diagnostic accent while all other graph geometry remains unchanged. Selecting
one reveals its CSR slot, endpoints, distance, weight, and oracle classification.

Allow switching to the coarse signal/background hierarchy. The classification
may change because the scientific question changed; the graph does not.

### Shot 08 — A/B construction

**Beat:** state precisely what changed.

Split into graph A and graph B with linked camera and selection. Common edges are
neutral; A-only and B-only edges receive symmetric contrasting encodings. The
comparison header lists the exact differing recipe fields and confirms all held
fixed inputs.

The user can unlock cameras, but the film's reset restores the linked preset.

### Shot 09 — Freeze adjacency, change weight

**Beat:** neighbor selection and coupling are different decisions.

Return to one adjacency artifact and compare two weight fields over it. A selected
edge is linked to a small kernel chart showing its distance and resulting weight.
Only appearance changes between weight fields; topology is constant and stated.

This shot may be deferred from the first executable cut if the graph adapter does
not yet expose a reusable adjacency/weight relation cleanly.

### Shot 10 — Geometry changes the reading

**Beat:** the same intrinsic eye may inhabit a different world.

Open Euclidean and Poincaré takes with explicitly corresponding structural
samples selected in both. Show the Poincaré boundary and native geodesic edges.
Backgrounds crossfade rather than morph. Metric and coordinate-space provenance
stay visible.

This is the final shot of the full pilot and may follow the minimal graph-lab
release as a second milestone.

### Shot 11 — Open laboratory

**Beat:** the film becomes a research state.

Stop autoplay with all artifacts available in the outline. The user may probe,
change hierarchy level, inspect diagnostics, select an existing take, or fork a
new graph recipe. Forking exits the canonical sequence but retains a return point.

## Interaction exits

- Selecting or probing preserves the current shot.
- Ordinary camera movement marks the camera as diverged but preserves playback.
- Changing a view encoding preserves the take and leaves the saved shot.
- Scrubbing a declared scientific axis preserves the study but leaves the shot.
- Editing a recipe creates a new take and exits the canonical film.
- `ReturnToFilm` restores view state only; it does not delete the experimental
  branch.

## Minimal executable cut

The first shippable slice comprises Shots 04–08 and 11:

1. load or generate one Euclidean EyeTorus take;
2. display the observed points;
3. adapt two complete graph results;
4. compare them in linked spatial panels;
5. disclose an oracle edge classification;
6. inspect samples and edges through stable identities;
7. save and reopen the exact study statically.

The anatomy opening, weight-isolation shot, Poincaré act, empirical-flow sequel,
and SPC sequel remain in the script but cannot delay proof of the core seam.

## Acceptance criteria

- Static and live modes expose the same study fingerprint and artifact inventory.
- Selecting an edge resolves to its original CSR slot after save/reload.
- A/B comparison declares held-fixed and changed recipe fields.
- Oracle toggling never changes graph fingerprints.
- Feature and display coordinates have distinct IDs even when numerically equal.
- Presentation thinning reports both displayed and full counts.
- Unsupported package versions fail visibly.
- No graph, metric, tangent, or oracle classification is implemented in
  JavaScript.
- The final workbench state remains interactive after playback.

## Sequels

- **EyeTorus II — The directions hidden in the eye:** empirical versus oracle
  tangent/line fields and their disagreement.
- **EyeTorus III — Heating the eye:** SPC temperature, edge currencies,
  assignments, curves, hierarchy, and lineage.
- **EyeTorus IV — Unlinked:** three-dimensional ambient false proximity versus
  the four-dimensional realization under explicit projections.
