// Presentation-only transforms and lighting for the Eye fixture. Source samples stay unchanged.
function advanceEyeMotion(state, dt) {
  // Integrate constant angular jerk exactly over this frame (degrees and seconds).
  state.angle = (state.angle + state.velocity * dt + state.currentAcceleration * dt * dt / 2
    + state.jerk * dt * dt * dt / 6) % 360;
  state.velocity += state.currentAcceleration * dt + state.jerk * dt * dt / 2;
  state.currentAcceleration += state.jerk * dt;
}

function eyeCoastStep(speed, damping, dt) {
  // Exact exponential drag: subdivision into animation frames does not change the travel.
  const decay = Math.exp(-damping * dt);
  return { speed: speed * decay, angle: damping > 0 ? speed * -Math.expm1(-damping * dt) / damping : speed * dt };
}

function createEyeViewer({ THREE, scene, camera, renderer, controls, root, appearance, frame, reset, recolor, photo, initialPointSize = 0.07 }) {
  const definitions = [
    { key: "central", name: "Iris", axis: "X", speed: 12 },
    { key: "upper", name: "Upper lid", axis: "Y", speed: -16 },
    { key: "lower", name: "Lower lid", axis: "Diagonal", speed: 9 },
    { key: "pupil", name: "Pupil", axis: "Y", speed: 20 },
    { key: "background", name: "Background", axis: "Z", speed: 0 },
  ];
  const fallbackColors = ["#00d2df", "#ffab00", "#8195ad", "#a78bfa", "#48505f"];
  const group = new THREE.Group();
  scene.add(group);
  const objects = new Map();
  let active = false, selected = "central", playing = !photo, collected = false, planeZ = 0;
  let lighting = { fill: 0.45, azimuth: -35, elevation: 45, background: "#0e0e12" };
  const light = new THREE.Vector3();
  const axis = new THREE.Vector3();
  const orientation = new THREE.Quaternion();
  const dragAxis = new THREE.Vector3();
  const dragRotation = new THREE.Quaternion();
  const inverseGroupRotation = new THREE.Quaternion();
  const coastRotation = new THREE.Quaternion();
  const motionDelta = new THREE.Quaternion();
  const sampleVelocity = new THREE.Vector3();
  const sampleAcceleration = new THREE.Vector3();
  const sampleJerk = new THREE.Vector3();
  const origin = new THREE.Vector3();
  const axes = { X: [1, 0, 0], Y: [0, 1, 0], Z: [0, 0, 1], Diagonal: [1, 1, 0.4] };

  const range = (id, label, min, max, step, value) => `<label class="eye-range"><span>${label}</span><input id="${id}" aria-label="${label}" type="range" min="${min}" max="${max}" step="${step}" value="${value}"><output for="${id}">${value}</output></label>`;
  root.innerHTML = `
    <div class="eye-toolbar"><button id="eye-front">Front</button><button id="eye-oblique">Oblique</button><button id="eye-fit">Fit</button></div>
    <label title="Choose whether dragging turns one component or the whole eye.">Manipulation <select id="eye-manipulation" aria-label="Manipulation"><option value="individual">Individual</option><option value="assembly">Assembly</option></select></label>
    <label title="Carry the motion of your hand after releasing an object."><input id="eye-inertia" type="checkbox" aria-label="Inertia" checked> Inertia</label>
    ${range("eye-damping", "Damping", 0, 5, 0.1, 1.5)}
    <p id="eye-hand-motion" class="eye-help">Drag and release to impart momentum.</p>
    <label class="eye-object-label">Object <select id="eye-object" aria-label="Object"></select></label>
    <div class="eye-identity"><label><input id="eye-visible" type="checkbox" aria-label="Show selected object"> Visible</label><label>Color <input id="eye-color" aria-label="Object color" type="color"></label></div>
    <details class="eye-section"><summary>Shape</summary><div id="eye-shape"></div></details>
    <details class="eye-section" open><summary>Appearance</summary>
      ${range("eye-size", "Point size", 0.015, 0.18, 0.005, 0.07)}
      ${range("eye-brightness", "Brightness", 0.1, 3, 0.05, 1)}
      ${range("eye-glow", "Glow", 0, 1.5, 0.05, 0.3)}
      ${range("eye-shine", "Shine", 0, 1, 0.05, 0.65)}
    </details>
    <details class="eye-section" open><summary>Motion</summary>
      <label><input id="eye-spin" type="checkbox" aria-label="Spin"> Spin</label>
      <label>Drag changes <select id="eye-drag-mode" aria-label="Drag changes"><option value="pose" selected>Pose</option><option value="speed">Speed</option><option value="acceleration">Acceleration</option><option value="jerk">Jerk</option></select></label>
      <div class="eye-select-row"><label>Axis <select id="eye-axis" aria-label="Spin axis"><option>X</option><option>Y</option><option>Z</option><option>Diagonal</option></select></label>
      <label>Pivot <select id="eye-pivot" aria-label="Rotation pivot"><option value="eye">Eye center</option><option value="object">Object center</option></select></label></div>
      ${range("eye-speed", "Speed °/s", -60, 60, 1, 12)}
      ${range("eye-acceleration", "Accel °/s²", -20, 20, 0.5, 0)}
      ${range("eye-jerk", "Jerk °/s³", -5, 5, 0.1, 0)}
      <p id="eye-motion-now" class="eye-help"></p>
      ${range("eye-phase", "Phase °", 0, 360, 1, 0)}
    </details>
    <details class="eye-section"><summary>Lighting & background</summary>
      ${range("eye-fill", "Fill light", 0, 1, 0.05, 0.45)}
      ${range("eye-light-azimuth", "Light direction °", -180, 180, 1, -35)}
      ${range("eye-light-elevation", "Light elevation °", -80, 80, 1, 45)}
      <label>Backdrop <input id="eye-backdrop" aria-label="Backdrop color" type="color" value="#0e0e12"></label>
    </details>
    <details class="eye-section"><summary>Sampling</summary><div id="eye-sampling"></div></details>
    <p class="eye-help">Drag naturally in Pose mode. Inertia carries your motion after release; turn it off to hold the pose.<br>Assembly turns the whole eye. Drag empty space to orbit · scroll to zoom.</p>
    <button id="eye-reset">Reset eye</button><div id="regen-status" role="status" aria-live="polite"></div>`;

  const el = id => document.getElementById(id);
  const playButton = el("eye-play");
  const collectButton = el("eye-collect");
  function stateFor(def) {
    let state = appearance.get(def.name);
    if (!state) { state = {}; appearance.set(def.name, state); }
    const defaults = { visible: true, size: initialPointSize, brightness: 1, glow: 0.3, shine: 0.65,
      axis: def.axis, speed: def.speed, acceleration: 0, jerk: 0,
      phase: 0, angle: 0, orientation: [0, 0, 0, 1], depthShift: 0,
      spin: def.speed !== 0, pivot: "eye" };
    for (const [key, value] of Object.entries(defaults)) state[key] ??= value;
    state.velocity ??= state.speed;
    state.currentAcceleration ??= state.acceleration;
    return state;
  }
  const current = () => stateFor(definitions.find(d => d.key === selected));
  const dragScales = { speed: 0.4, acceleration: 0.15, jerk: 0.03 };
  let dragMode = "pose";
  let moveTogether = false;
  let inertia = true, damping = 1.5;
  const makeMotion = () => ({ velocity: new THREE.Vector3(), acceleration: new THREE.Vector3(), jerk: new THREE.Vector3() });
  const assemblyMotion = makeMotion();
  const manualMotions = new Map(definitions.map(def => [def.key, makeMotion()]));
  const motionEpsilon = 0.002;
  let drag = null;
  const projected = new THREE.Vector3();
  const canvas = renderer.domElement;
  function clearMotion(motion) {
    motion.velocity.set(0, 0, 0); motion.acceleration.set(0, 0, 0); motion.jerk.set(0, 0, 0);
  }
  function stopManualMotion() {
    clearMotion(assemblyMotion);
    for (const motion of manualMotions.values()) clearMotion(motion);
  }
  function hasMomentum() {
    return assemblyMotion.velocity.length() > motionEpsilon ||
      [...manualMotions.values()].some(motion => motion.velocity.length() > motionEpsilon);
  }
  function refreshPlaybackButton() {
    const moving = playing || hasMomentum();
    const label = moving ? "Pause motion" : "Play motion";
    if (playButton.textContent !== label) playButton.textContent = label;
    playButton.setAttribute("aria-pressed", String(moving));
  }
  function recordDragMotion(quaternion, timestamp) {
    const dt = (timestamp - drag.lastTime) / 1000;
    const motion = drag.motion;
    motionDelta.copy(drag.lastOrientation).invert().premultiply(quaternion).normalize();
    // q and -q are the same pose; choose the shortest angular displacement.
    if (motionDelta.w < 0) motionDelta.set(-motionDelta.x, -motionDelta.y, -motionDelta.z, -motionDelta.w);
    const sine = Math.hypot(motionDelta.x, motionDelta.y, motionDelta.z);
    const angle = 2 * Math.atan2(sine, Math.max(0, motionDelta.w));
    if (dt > 0) {
      sampleVelocity.set(motionDelta.x, motionDelta.y, motionDelta.z);
      if (sine > 1e-8) sampleVelocity.multiplyScalar(Math.min(angle / dt, 12) / sine);
      else sampleVelocity.set(0, 0, 0);
      const blend = -Math.expm1(-dt / 0.03);
      sampleVelocity.lerpVectors(motion.velocity, sampleVelocity, blend);
      sampleAcceleration.copy(sampleVelocity).sub(motion.velocity).divideScalar(dt);
      sampleAcceleration.lerpVectors(motion.acceleration, sampleAcceleration, blend);
      sampleJerk.copy(sampleAcceleration).sub(motion.acceleration).divideScalar(dt);
      motion.jerk.lerp(sampleJerk, blend);
      motion.acceleration.copy(sampleAcceleration);
      motion.velocity.copy(sampleVelocity);
    }
    if (angle > 1e-6) drag.lastMovement = timestamp;
    drag.lastTime = timestamp;
    drag.lastOrientation.copy(quaternion);
  }
  // Pick the displayed points, including the shader's pivot and Collect translation.
  function pickObject(event) {
    const rect = canvas.getBoundingClientRect();
    const x = event.clientX - rect.left, y = event.clientY - rect.top;
    let best = null, nearest = Infinity;
    scene.updateMatrixWorld(true);
    for (const def of definitions) {
      const obj = objects.get(def.key);
      if (!obj?.core.visible) continue;
      const positions = obj.geometry.attributes.position;
      for (let i = 0; i < positions.count; i++) {
        projected.fromBufferAttribute(positions, i).sub(obj.core.material.uniforms.pivot.value)
          .applyMatrix4(obj.core.matrixWorld).project(camera);
        if (projected.z < -1 || projected.z > 1 || projected.z >= nearest) continue;
        const dx = (projected.x + 1) * rect.width / 2 - x;
        const dy = (1 - projected.y) * rect.height / 2 - y;
        if (dx * dx + dy * dy <= 144) { nearest = projected.z; best = def; }
      }
    }
    return best;
  }
  function setValue(key, value) {
    if (key === "phase") setCollected(false, true);
    const state = current();
    state[key] = value;
    if (key === "speed") state.velocity = value;
    if (key === "acceleration") state.currentAcceleration = value;
    syncRange("eye-" + key, value);
    showMotion(state);
  }
  function endDrag(event) {
    if (!drag || (event && event.pointerId !== drag.pointerId)) return;
    const pointerId = drag.pointerId;
    if (drag.mode === "pose") {
      const idle = event ? Math.max(0, (event.timeStamp - drag.lastMovement) / 1000) : Infinity;
      if (!inertia || event?.type !== "pointerup" || idle > 0.15) clearMotion(drag.motion);
      else drag.motion.velocity.multiplyScalar(Math.exp(-idle / 0.08));
    }
    controls.enabled = drag.orbitEnabled;
    drag = null;
    canvas.style.cursor = "";
    if (canvas.hasPointerCapture(pointerId)) canvas.releasePointerCapture(pointerId);
    refreshPlaybackButton();
  }
  canvas.addEventListener("pointerdown", event => {
    if (!active || event.button !== 0 || !event.isPrimary || drag) return;
    const def = pickObject(event);
    if (!def) return;
    if (!moveTogether) {
      el("eye-object").value = def.key;
      selectObject();
    }
    const state = current();
    const mode = moveTogether || !playing || !state.spin ? "pose" : dragMode;
    const motion = moveTogether ? assemblyMotion : manualMotions.get(def.key);
    if (mode === "pose") {
      if (moveTogether) stopManualMotion();
      else { clearMotion(assemblyMotion); clearMotion(motion); }
    }
    drag = { pointerId: event.pointerId, x: event.clientX, y: event.clientY, mode, key: def.key,
      together: moveTogether,
      orientation: moveTogether ? group.quaternion.clone() : new THREE.Quaternion().fromArray(state.orientation),
      lastOrientation: moveTogether ? group.quaternion.clone() : new THREE.Quaternion().fromArray(state.orientation),
      motion, lastTime: event.timeStamp, lastMovement: event.timeStamp,
      start: mode === "speed" ? state.velocity : mode === "acceleration" ? state.currentAcceleration : state.jerk,
      orbitEnabled: controls.enabled };
    controls.enabled = false;
    canvas.style.cursor = "grabbing";
    canvas.setPointerCapture(event.pointerId);
    event.preventDefault();
    event.stopImmediatePropagation();
  }, true);
  canvas.addEventListener("pointermove", event => {
    if (!drag || event.pointerId !== drag.pointerId) return;
    const dx = event.clientX - drag.x, dy = event.clientY - drag.y;
    if (drag.mode === "pose") {
      const distance = Math.hypot(dx, dy);
      if (distance > 0) {
        // Leave Collect without moving any component out of its held position.
        setCollected(false, true);
        dragAxis.set(dy, dx, 0).normalize().applyQuaternion(camera.quaternion);
        if (!drag.together) dragAxis.applyQuaternion(inverseGroupRotation.copy(group.quaternion).invert());
        dragRotation.setFromAxisAngle(dragAxis, distance * 0.008).multiply(drag.orientation).normalize();
      } else dragRotation.copy(drag.orientation);
      recordDragMotion(dragRotation, event.timeStamp);
      if (drag.together) group.quaternion.copy(dragRotation);
      else dragRotation.toArray(current().orientation);
    } else {
      const input = el("eye-" + drag.mode), step = Number(input.step);
      const value = Math.max(Number(input.min), Math.min(Number(input.max),
        Math.round((drag.start + (dx - dy) * dragScales[drag.mode]) / step) * step));
      setValue(drag.mode, value);
    }
    event.preventDefault();
    event.stopImmediatePropagation();
  }, true);
  canvas.addEventListener("pointerup", endDrag);
  canvas.addEventListener("pointercancel", endDrag);
  canvas.addEventListener("lostpointercapture", endDrag);
  function syncRange(id, value) {
    const input = el(id);
    input.value = value;
    input.nextElementSibling.value = Number(value).toFixed(Number(input.step) < 1 ? 2 : 0);
  }
  function selectObject() {
    selected = el("eye-object").value;
    const state = current();
    el("eye-visible").checked = state.visible;
    el("eye-color").value = state.color;
    el("eye-spin").checked = state.spin;
    el("eye-axis").value = state.axis;
    el("eye-pivot").value = state.pivot;
    for (const key of ["size", "brightness", "glow", "shine", "speed", "acceleration", "jerk", "phase"])
      syncRange("eye-" + key, state[key]);
    showMotion(state);
    syncDragMode();
    for (const section of el("eye-shape").children) section.hidden = section.dataset.eyeComponent !== selected;
  }
  function attachGeometry(body) {
    el("eye-shape").replaceChildren();
    el("eye-sampling").replaceChildren();
    for (const section of Array.from(body.querySelectorAll("details[data-eye-component]"))) {
      section.open = true;
      (section.dataset.eyeComponent === "sampling" ? el("eye-sampling") : el("eye-shape")).appendChild(section);
    }
    const regen = body.querySelector("#btn-regen");
    if (regen) el("eye-sampling").appendChild(regen);
    for (const section of el("eye-shape").children) section.hidden = section.dataset.eyeComponent !== selected;
  }
  function syncDragMode() {
    const frozen = !playing || !current().spin;
    const input = el("eye-drag-mode");
    input.value = moveTogether || frozen ? "pose" : dragMode;
    input.disabled = moveTogether || frozen;
    input.title = moveTogether ? "Assembly rotates the whole eye. Object settings still apply to the selected object."
      : frozen ? "Paused objects drag their pose. Play motion and enable Spin to adjust motion by dragging."
      : "Choose what dragging an object changes.";
  }
  function setCollected(value, preservePose = false) {
    if (collected && !value) {
      for (const def of definitions) {
        const obj = objects.get(def.key);
        stateFor(def).depthShift = preservePose && obj && def.key !== "background" ? planeZ - obj.offset.z : 0;
      }
    }
    collected = value;
    collectButton.setAttribute("aria-pressed", String(collected));
    collectButton.textContent = collected ? "Release" : "Collect";
    collectButton.title = collected
      ? "Release the shared plane; motion stays paused."
      : "Align the eye components with the iris plane and hold them there.";
  }
  function setPlaying(value) {
    if (value) setCollected(false, true);
    else { endDrag(); stopManualMotion(); }
    playing = value;
    refreshPlaybackButton();
    syncDragMode();
  }
  function resetPose() {
    endDrag();
    setCollected(false);
    group.quaternion.identity();
    for (const def of definitions) {
      const s = stateFor(def);
      s.angle = 0; s.phase = 0;
      s.orientation = [0, 0, 0, 1]; s.depthShift = 0;
      s.velocity = s.speed; s.currentAcceleration = s.acceleration;
    }
    setPlaying(false);
    syncRange("eye-phase", 0);
  }
  el("eye-object").addEventListener("change", selectObject);
  el("eye-manipulation").addEventListener("change", event => {
    endDrag();
    moveTogether = event.target.value === "assembly";
    syncDragMode();
  });
  el("eye-inertia").addEventListener("change", event => {
    inertia = event.target.checked;
    if (!inertia) stopManualMotion();
    refreshPlaybackButton();
  });
  el("eye-damping").title = "How quickly a released drag slows down. Zero keeps its momentum.";
  el("eye-damping").addEventListener("input", event => {
    damping = Number(event.target.value);
    syncRange("eye-damping", damping);
  });
  el("eye-visible").addEventListener("change", event => { current().visible = event.target.checked; });
  el("eye-color").addEventListener("input", event => { current().color = event.target.value; recolor(); });
  el("eye-spin").addEventListener("change", event => {
    current().spin = event.target.checked;
    if (!event.target.checked) { clearMotion(manualMotions.get(selected)); clearMotion(assemblyMotion); }
    syncDragMode();
    refreshPlaybackButton();
  });
  el("eye-drag-mode").addEventListener("change", event => { dragMode = event.target.value; });
  el("eye-axis").addEventListener("change", event => { current().axis = event.target.value; });
  el("eye-pivot").addEventListener("change", event => { current().pivot = event.target.value; });
  for (const key of ["size", "brightness", "glow", "shine", "speed", "acceleration", "jerk", "phase"])
    el("eye-" + key).addEventListener("input", event => {
      setValue(key, Number(event.target.value));
    });
  el("eye-speed").title = "Set the current angular speed. Collect restarts from this value.";
  el("eye-acceleration").title = "Set the current angular acceleration. Collect restarts from this value.";
  el("eye-jerk").title = "Rate of change of acceleration. Zero holds the current acceleration.";
  for (const [id, key] of [["eye-fill", "fill"], ["eye-light-azimuth", "azimuth"], ["eye-light-elevation", "elevation"]])
    el(id).addEventListener("input", event => { lighting[key] = Number(event.target.value); syncRange(id, lighting[key]); });
  el("eye-backdrop").addEventListener("input", event => { lighting.background = event.target.value; renderer.setClearColor(lighting.background); });
  el("eye-front").addEventListener("click", () => frame(0, 0));
  el("eye-oblique").addEventListener("click", () => frame(0.5, 0.35));
  el("eye-fit").addEventListener("click", () => frame(null, null));
  el("eye-reset").addEventListener("click", () => { resetPose(); reset(); });
  playButton.addEventListener("click", () => setPlaying(!(playing || hasMomentum())));
  collectButton.addEventListener("click", () => {
    if (collected) { setCollected(false); return; }
    resetPose();
    setCollected(true);
  });
  setCollected(false);
  setPlaying(playing);

  const vertexShader = `
    uniform float pointSize;
    uniform float viewportScale;
    uniform float halo;
    uniform vec3 pivot;
    varying vec3 surfaceNormal;
    void main() {
      vec4 mvPosition = modelViewMatrix * vec4(position - pivot, 1.0);
      surfaceNormal = normalize(normalMatrix * normal);
      gl_Position = projectionMatrix * mvPosition;
      gl_PointSize = clamp(pointSize * viewportScale * projectionMatrix[1][1] / max(0.01, -mvPosition.z), 1.0, 96.0) * (1.0 + 2.0 * halo);
    }`;
  const fragmentShader = `
    uniform vec3 tint;
    uniform vec3 lightDirection;
    uniform float brightness;
    uniform float glow;
    uniform float shine;
    uniform float fillLight;
    uniform float halo;
    varying vec3 surfaceNormal;
    void main() {
      vec2 p = gl_PointCoord * 2.0 - 1.0;
      p.y = -p.y;
      float rr = dot(p, p);
      if (rr > 1.0) discard;
      if (halo > 0.5) {
        float alpha = glow * 0.18 * exp(-4.0 * rr) * (1.0 - smoothstep(0.65, 1.0, rr));
        gl_FragColor = vec4(tint * brightness, alpha);
      } else {
        vec3 beadNormal = vec3(p, sqrt(max(0.0, 1.0 - rr)));
        vec3 halfVector = normalize(lightDirection + vec3(0.0, 0.0, 1.0));
        float bead = 0.35 + 0.65 * max(dot(beadNormal, lightDirection), 0.0);
        float form = 0.35 + 0.65 * max(dot(normalize(surfaceNormal), lightDirection), 0.0);
        float specular = pow(max(dot(beadNormal, halfVector), 0.0), 24.0);
        specular += 0.35 * pow(max(dot(normalize(surfaceNormal), halfVector), 0.0), 18.0);
        vec3 lit = tint * brightness * (fillLight + 1.25 * bead * form);
        lit += vec3(shine * specular * 1.3);
        gl_FragColor = vec4(lit, 1.0);
      }
      #include <tonemapping_fragment>
      #include <colorspace_fragment>
    }`;
  function material(halo) {
    return new THREE.ShaderMaterial({
      uniforms: {
        tint: { value: new THREE.Color() }, pointSize: { value: 0.07 }, viewportScale: { value: 360 },
        brightness: { value: 1 }, glow: { value: 0.3 }, shine: { value: 0.65 }, fillLight: { value: 0.45 },
        lightDirection: { value: new THREE.Vector3() }, pivot: { value: new THREE.Vector3() }, halo: { value: halo ? 1 : 0 },
      }, vertexShader, fragmentShader,
      transparent: halo, depthWrite: !halo,
      blending: halo ? THREE.AdditiveBlending : THREE.NormalBlending,
    });
  }
  function clearObjects() {
    for (const obj of objects.values()) {
      group.remove(obj.core, obj.halo);
      obj.geometry.dispose(); obj.core.material.dispose(); obj.halo.material.dispose();
    }
    objects.clear();
  }
  function setScene(data) {
    endDrag();
    active = data.generator_params?.generator === "EyeTorusToy";
    clearObjects();
    group.visible = active;
    root.hidden = !active;
    playButton.hidden = collectButton.hidden = !active;
    if (!active) {
      el("eye-shape").replaceChildren(); el("eye-sampling").replaceChildren();
      renderer.setClearColor(0x0e0e12); return;
    }
    renderer.setClearColor(lighting.background);
    const gp = data.generator_params;
    planeZ = gp.centralOffset?.[2] ?? 0;
    const labels = (data.label_layers.find(layer => layer.kind === "GroundTruth") ?? data.label_layers[0])?.labels ?? [];
    const names = gp.componentNames ?? { 0: "Iris", 1: "Upper lid", 2: "Lower lid", ...(gp.pupilPoints > 0 ? { 3: "Pupil" } : {}), [gp.backgroundLabel]: "Background" };
    const labelNames = new Map(Object.entries(names).map(([id, name]) => [Number(id), name]));
    const buckets = new Map(definitions.map(def => [def.name, []]));
    for (let i = 0; i < labels.length; i++) buckets.get(labelNames.get(labels[i]))?.push(i);
    const raw = data.points.features, d = data.points.d;
    el("eye-object").replaceChildren();
    definitions.forEach((def, index) => {
      const state = stateFor(def);
      const id = Object.keys(names).find(key => names[key] === def.name);
      state.color ??= gp.componentColors?.[id] ?? fallbackColors[index];
      const indices = buckets.get(def.name);
      const option = document.createElement("option");
      option.value = def.key; option.textContent = `${def.name} · ${indices.length.toLocaleString()}`;
      el("eye-object").appendChild(option);
      if (!indices.length) return;
      const offset = new THREE.Vector3().fromArray(gp[def.key + "Offset"] ?? [0, 0, 0]);
      const positions = new Float32Array(indices.length * 3), normals = new Float32Array(indices.length * 3);
      const normal = new THREE.Vector3();
      indices.forEach((sample, row) => {
        const x = raw[sample * d], y = d > 1 ? raw[sample * d + 1] : 0, z = d > 2 ? raw[sample * d + 2] : 0;
        positions.set([x, y, z], row * 3);
        const dx = x - offset.x, dy = y - offset.y, dz = z - offset.z;
        if (index < 3) {
          const radial = Math.hypot(dx, dy) || 1;
          const major = gp[def.key + "MajorR"] ?? radial;
          normal.set(dx / radial * (radial - major), dy / radial * (radial - major), dz / Math.max(0.01, (gp.zThickness ?? 1) ** 2));
        } else normal.set(dx, dy, dz);
        if (normal.lengthSq() < 1e-12) normal.set(0, 0, 1);
        normal.normalize().toArray(normals, row * 3);
      });
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
      geometry.setAttribute("normal", new THREE.BufferAttribute(normals, 3));
      const core = new THREE.Points(geometry, material(false));
      const halo = new THREE.Points(geometry, material(true));
      // GPU pivot transforms change bounds; keep these five modest point batches unculled.
      core.frustumCulled = halo.frustumCulled = false;
      group.add(core, halo);
      objects.set(def.key, { core, halo, geometry, offset });
    });
    el("eye-object").value = selected;
    selectObject();
  }
  function showMotion(state) {
    const status = !playing ? "Auto paused · resume" : !state.spin ? "Spin off · resume" : "Auto";
    const text = `${status}: ${state.velocity.toFixed(1)} °/s · ${state.currentAcceleration.toFixed(1)} °/s²`;
    const output = el("eye-motion-now");
    if (output.textContent !== text) output.textContent = text;
  }
  function coast(motion, quaternion, dt) {
    const speed = motion.velocity.length();
    if (speed <= motionEpsilon) { clearMotion(motion); return; }
    const next = eyeCoastStep(speed, damping, dt);
    coastRotation.setFromAxisAngle(axis.copy(motion.velocity).divideScalar(speed), next.angle);
    quaternion.premultiply(coastRotation).normalize();
    motion.velocity.multiplyScalar(next.speed / speed);
    motion.acceleration.copy(motion.velocity).multiplyScalar(-damping);
    motion.jerk.copy(motion.velocity).multiplyScalar(damping * damping);
  }
  function showHandMotion() {
    const motion = moveTogether ? assemblyMotion : manualMotions.get(selected);
    const held = drag?.mode === "pose" && (drag.together === moveTogether) && (moveTogether || drag.key === selected);
    const speed = motion.velocity.length() * 180 / Math.PI;
    const status = held ? "Drag" : speed > motionEpsilon * 180 / Math.PI ? "Coast" : "At rest";
    const text = status === "At rest" ? "At rest · drag to rotate" :
      `${status}: ${speed.toFixed(1)} °/s · Accel ${(motion.acceleration.length() * 180 / Math.PI).toFixed(0)} °/s² · Jerk ${(motion.jerk.length() * 180 / Math.PI).toFixed(0)} °/s³`;
    const output = el("eye-hand-motion");
    if (output.textContent !== text) output.textContent = text;
  }
  function tick(dt) {
    if (!active) return;
    if (drag?.mode === "pose" && performance.now() - drag.lastMovement > 150) clearMotion(drag.motion);
    if (inertia && !(drag?.mode === "pose" && drag.together) && dt > 0) coast(assemblyMotion, group.quaternion, dt);
    const az = lighting.azimuth * Math.PI / 180, elevation = lighting.elevation * Math.PI / 180;
    light.set(Math.sin(az) * Math.cos(elevation), Math.sin(elevation), Math.cos(az) * Math.cos(elevation));
    light.transformDirection(camera.matrixWorldInverse);
    for (const def of definitions) {
      const obj = objects.get(def.key);
      if (!obj) continue;
      const state = stateFor(def);
      const held = drag?.mode === "pose" && (drag.together || drag.key === def.key);
      if (playing && state.spin && !held && dt > 0) advanceEyeMotion(state, dt);
      if (inertia && !held && dt > 0) {
        orientation.fromArray(state.orientation);
        coast(manualMotions.get(def.key), orientation, dt);
        orientation.toArray(state.orientation);
      }
      axis.fromArray(axes[state.axis]).normalize();
      for (const points of [obj.core, obj.halo]) {
        points.visible = state.visible && (points === obj.core || state.glow > 0);
        points.quaternion.setFromAxisAngle(axis, (state.angle + state.phase) * Math.PI / 180);
        points.quaternion.premultiply(orientation.fromArray(state.orientation));
        points.position.copy(state.pivot === "object" ? obj.offset : origin);
        const u = points.material.uniforms;
        u.pivot.value.copy(points.position);
        // Apply the collected depth before rotation, including the object's pivot.
        const depthShift = collected && def.key !== "background" ? planeZ - obj.offset.z : state.depthShift;
        if (state.pivot === "object") points.position.z += depthShift;
        else u.pivot.value.z -= depthShift;
        u.tint.value.set(state.color ?? fallbackColors[definitions.indexOf(def)]);
        u.pointSize.value = state.size;
        u.viewportScale.value = renderer.domElement.height * 0.5;
        u.brightness.value = state.brightness; u.glow.value = state.glow; u.shine.value = state.shine;
        u.fillLight.value = lighting.fill; u.lightDirection.value.copy(light);
      }
    }
    showMotion(current());
    showHandMotion();
    refreshPlaybackButton();
  }
  return { setScene, attachGeometry, tick, setPlaying };
}
