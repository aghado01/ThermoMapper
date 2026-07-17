# Visualization screenplay format

A screenplay is a human-reviewable specification for a reproducible film over a
`VizStudy`. It directs computation, view state, disclosure, and narration without
placing scientific state in the renderer.

The first scripts are Markdown. A typed/serializable film contract should emerge
only after two or more scripts demonstrate the stable grammar.

## Vocabulary boundary

The production vocabulary is an authoring-layer device for scripts and planning
prose. Typed contracts use idiomatic technical names, and this section records
the mapping as each surface is minted:

- **take** → `RunDescriptor` / `RunId` — a pinned recipe execution or branch;
- **film**, **shot**, **cue** — unminted; technical names are chosen when the
  typed contract emerges from stable scripts.

The metaphor may also surface in the front end as user-facing strings and
visual affordances — labels, tooltips, playback chrome — where it earns its
keep as intuitive UI. It stays out of typed identifiers there too: a component
named `SequencePlayer` may render the label "Film".

## Required header

Every film declares:

```text
Film ID          stable semantic ID
Title            displayed title
Status           treatment | blocked | executable | canonical
Purpose          diagnostic | exploratory | demonstrative (one or more)
Question         the single scientific question organizing the film
Audience         assumed background
Recipe           canonical or intervened operation graph
Takes            pinned inputs, seeds, configs, and branch identities
Evidence policy  which oracle artifacts exist and when they may appear
Panels           required spatial and non-spatial views
Axes             scientific axes used by the film
Runtime modes    static replay, live execution, or both
```

## Production vocabulary

- **Study:** durable collection of entities, artifacts, relations, panels,
  provenance, presets, and films.
- **Recipe:** typed operation graph that produces scientific artifacts.
- **Take:** one pinned recipe execution or branch with immutable provenance.
- **Shot:** a saved arrangement of panels and durable view state.
- **Cue:** an ordered transition or action connecting shots.
- **Beat:** the scientific idea the audience should understand.
- **Sequence:** ordered shots and cues, possibly spanning related takes.
- **Probe:** an interactive inspection that does not alter scientific artifacts.
- **Intervention:** a recorded recipe change that can alter scientific artifacts.

The renderer may know scenes and frames internally. Those are implementation
terms and must not be confused with a study, take, or scientific-axis frame.

## The three clocks

Each cue identifies exactly one primary clock:

1. **Pipeline clock** — a real stage boundary or exposed algorithm event.
2. **Scientific clock** — movement along temperature, filtration, sweep, or
   another domain axis.
3. **Cinematic clock** — presentation-only camera, fade, reveal, or interpolation.

A cue may coordinate clocks but cannot silently substitute one for another. If a
graph is progressively revealed after being computed in one batch, the cue is
cinematic, not an execution trace. If displayed intermediate states are sampled
from a solver, the producer must have emitted them as artifacts.

## Shot specification

Each shot records:

```text
Shot ID
Beat
Take/artifact bindings
Panel layout
Camera and projection
Visible layers and encodings
Selection and scientific-axis cursor
Annotations and narration
Evidence disclosure
Entry cue / exit cue
Allowed interactions
Reset behavior
Performance policy
```

`Performance policy` records presentation-derived thinning, level of detail, or
sampling. A view may draw fewer edges than the artifact contains, but it must say
so and diagnostics must continue to describe the full artifact.

## Cue specification

Each cue records:

```text
Cue ID
Primary clock
Trigger          autoplay | user action | axis change | artifact arrival
From / to
Duration/easing  cinematic cues only
Action           typed view action or typed recipe intervention
Truth note       what changes scientifically, and what merely changes visually
Interruptible
```

Canonical typed view actions mirror `Viz.Scene.ViewAction`; the code is
authoritative:

- `SelectEntities`
- `SetLayerVisibility`
- `SetLayerEncodings`
- `SetCamera`
- `SetDisplayCoordinates`
- `SetScientificAxisCursor`
- `FrameSelection`

Actions scripts anticipate that are not yet typed (`SetProjection`,
`OpenInspector`) are minted in `Viz.Scene` when a script first requires them.

Recipe interventions are never encoded as view actions.

## Evidence disclosure

Every visible artifact declares one of the study evidence roles. Scripts state:

- whether the canonical recipe could consume it;
- whether it is visible initially;
- which cue first discloses it;
- how it is visually distinguished;
- whether a derived comparison changes the result or only scores it.

Oracle content should generally arrive after the observed result in a
demonstration. Diagnostic films may start with it visible if the title and UI
make that posture unmistakable.

## Interaction contract

A film is not a video. At any shot, the user may be allowed to:

- orbit, zoom, and pan;
- hover or select scientific entities;
- inspect provenance and raw values;
- toggle declared layers;
- scrub a scientific axis;
- compare takes;
- fork a new experimental branch.

The script specifies which actions preserve the sequence and which leave it. On
leaving, the current research state is retained; `ReturnToFilm` restores the
saved shot without discarding newly created branches.

## Probe admissibility

A probe reads; it does not compute. In static replay a probe may only read
stored artifacts — values, relations, provenance, and diagnostics the producers
already emitted. Computation-backed probes (ad-hoc distances, neighborhood
re-evaluation) are live-mode capabilities served by the host's job protocol, or
are pinned at production time: a script that needs a computed probe in static
mode pins the probed selection and emits the probe result as a study artifact.

## Truth and honesty checklist

Before a film becomes canonical:

- The narration does not call a cinematic reveal an algorithm trace.
- Every metric, graph, field, and result comes from its owning public producer.
- Oracle artifacts are typed and disclosed.
- Cross-take point or entity correspondence is explicit; visual proximity is not
  treated as identity.
- Coordinate projections and geometry models are visible in provenance.
- Native geometry is used when straight Euclidean rendering would mislead.
- Thinning, interpolation, smoothing, and aggregation are marked
  `PresentationDerived`.
- A/B comparisons state what was held fixed.
- Scientific-axis direction and units are visible.
- Static playback uses stored artifacts rather than silently recomputing.
- The final shot can be inspected as a normal workbench state.

## Suggested Markdown shape

```markdown
# film title

## Header
## Scientific contract
## Cast and evidence
## Takes
## Panel plan
## Sequence
### Shot ...
### Shot ...
## Interaction exits
## Required producers and adapters
## Acceptance criteria
## Open scientific questions
```

Scripts should read well as prose while remaining precise enough to implement.
Tables are useful for short shot lists; complex sequences should give each shot
its own section rather than compressing the scientific argument into cells.
