## Summary: OpenDrift vertical behavior + LarvalFish and HAB minimal feature decisions

### Scope and constraints

* Target forcing: **3D ROMS** (Cook Inlet), available fields: **u, v, w, temp, salinity, zeta**.
* **No biology** in HAB model (transport-focused only).
* **No viability limits / no killing**, **no cyst/germination**, **no settlement**.
* **One species per run**.
* ROMS output **does not include vertical diffusivity for tracers**.

---

## Vertical behavior design shared by HAB + larvae

### Behavior modes (final)

Only two user-selectable vertical behavior modes:

1. **Depth**: user provides a single preferred depth `z_pref`.
2. **DVM**: user provides two preferred depths `z_day` and `z_night` (day/night schedule).

No explicit user-input depth bands; bands are derived internally.

### Unified vertical effort parameter

* Use a single parameter name across all taxa and both modes: **`w_active`** (m/s).
* Interpretation:

  * literal swimming for larvae and motile dinoflagellates
  * **effective vertical positioning tendency** for weakly/non-motile taxa (e.g., Pseudo-nitzschia), not necessarily literal swimming.

### Internal depth-band expansion from a single depth (updated)

* Any user depth (`z_pref`, `z_day`, `z_night`) is expanded internally into a band using a global rule:

  * `dz = clamp(dz_min, dz_rel * z, dz_max)`
  * chosen values: **`dz_min = 1 m`**, **`dz_rel = 0.1`**, **`dz_max = 15 m`**
  * internal band: `[z - dz, z + dz]`
* When constructing the band, **validate and clip only to the physical water column**:

  * top = sea surface (0 or `zeta`-referenced surface depending on internal convention)
  * bottom = seabed depth at that location
  * No additional surface/bottom buffers are used; just ensure requested depths/bands fall within `[top, bottom]`.

### DVM mechanics (agreed approach)

* DVM selects the active band based on day/night (use what was written into harmfulalgalbloom.py in models).
* Behavioral vertical velocity term:

  * If particle is inside the current band: `w_behavior = 0`
  * If particle is shallower than band: move downward with magnitude up to `w_active`
  * If particle is deeper than band: move upward with magnitude up to `w_active`
* `w_active` used symmetrically up/down (no separate ascent/descent speeds).

---

## LarvalFish lifecycle decisions

### Egg stage: hatching only 

* Egg stage is included only where appropriate (pollock, cod; **no egg stage** for herring, rockfish, razor clams; salmon omitted).
* Eggs are used primarily as a **time delay until hatch**.


### Hatch timing 

* For v1, use **fixed time-to-hatch only**:

  * Parameter: `hatch_time_hours` (or days, but store internally as seconds/hours)


### Larval stages

* After hatch, treat larvae using the same **Depth** or **DVM** vertical behavior framework as phytoplankton/HAB particles.

---

## HAB model decisions

* HAB model is transport-only; no biology.
* Use the same vertical behavior system:

  * Motile dinoflagellates: typically **DVM + `w_active`**
  * Pseudo-nitzschia (PN): **Depth + `w_active`** (effective positioning, not literal swimming)

---

## Species list and recommended behavior mapping

Included:

* HAB: **Alexandrium**, **Dinophysis**, **Pseudo-nitzschia (PN)**
* Larvae: **Walleye pollock**, **Pacific cod**, **Pacific herring**, **Rockfish**, **Razor clams**
* Skipped: salmon; capelin/sand lance/flatfishes removed.

Behavior choices:

* Alexandrium: DVM
* Dinophysis: DVM
* PN: Depth (with `w_active` as effective positioning)
* Pollock: Egg stage yes; Larva: Depth (DVM optional later)
* Cod: Egg stage yes; Larva: Depth
* Herring: Egg stage no; Larva: Depth
* Rockfish: Egg stage no; Larva: Depth
* Razor clams: Egg stage no; Larva: Depth

---

## Defaults and user-facing simplicity

User-facing inputs per run:

* **Depth mode**: `z_pref`, `w_active`
* **DVM mode**: `z_day`, `z_night`, `w_active`
* **If egg stage is used**: `hatch_time_hours`

Global defaults (not per species):

* band expansion: `dz_min=1 m`, `dz_rel=0.1`, `dz_max=15 m`
* DVM schedule: simple fixed day/night hours (implementation-defined default)

Suggested starting-point presets for UI:

* HAB motile dinoflagellate: DVM; `z_night≈5 m`, `z_day≈25 m`, `w_active≈0.001 m/s`
* HAB PN: Depth; `z_pref≈5 m`, `w_active≈0.0002 m/s`
* Fish larvae: Depth; `z_pref≈10–15 m`, `w_active≈0.003 m/s`
* Invertebrate larvae: Depth; `z_pref≈5 m`, `w_active≈0.001 m/s`

---

## Implementation guidance for the coder (high level)

* Implement shared vertical behavior logic (Depth/DVM) in a reusable location so LarvalFish and HAB can call it.
* Implement internal band expansion from user depths using the global `dz` rule, and clip bands to the **physical** water column only (surface to seabed).
* Apply `w_active` as a deterministic vertical tendency driving particles back into the currently active band.
* In LarvalFish, implement egg stage as a timer using fixed `hatch_time_hours`; eggs are passive until they convert to larvae; larvae then use the shared Depth/DVM logic.
* Keep everything single-species-per-run; no per-particle species registry required for v1.

## Other notes

* the larvalfish.py model exists in the opendrift code base already (harmfulalgalbloom.py does not)
* try to modify only larvalfish.py to implement this behavior, but allow for the original behavior to be used as a user selection where possible. OpenDrift is someone else's codebase and I want to make minimal and only useful changes that dont impact the currently-available behavior if possible. We could make a new model if necessary.
* include references where possible but make sure they are correct and real
* I think both phytoplankton/hab and larvalfish scenarios can be implemented in larvalfish with different options

## changes to particle-tracking-manager (PTM)

* regardless of whether phytoplankton and larvalfish are combined, make there be two separate scenarios in PTM
* each need to have their own json schema that can be exported