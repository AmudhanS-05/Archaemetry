# archaemetry: CPU lane (mesh to drawing)

Turns a 3D mesh of a sherd or vessel into a publication-style profile drawing with a
GO / REVIEW verdict. Pure Python + numpy stack (numpy, scipy, trimesh, matplotlib,
PyYAML): no GPU, runs on Windows, macOS and Linux, Python 3.10-3.13.

```
mesh ─▶ cleanup ─▶ axis ─▶ section cut + thickness ─▶ drawing (SVG/PDF/PNG) + report
        (largest    (pca / rim_arc /   (plane through the axis,   (filled section, axis,
         piece)      normals_refine /   perpendicular thickness)   rim line, exterior view,
                     manual)                                       block scale bar)
```

## Use it

**Field testers:** download the desktop app from the repository's *Actions* or
*Releases* page and read `packaging/FIELD_TEST.md`. No Python needed.

**From source:**
```bash
pip install ./cpu                    # or: pip install -e "./cpu[test]"
archaemetry selftest                 # checks the install against known shapes
archaemetry run sherd.ply -o out --id SAL_1.1            # rim sherd (default: auto axis)
archaemetry run vessel.obj -o out --axis pca             # complete vessel
archaemetry run sherd.ply -c config/sherd.example.yaml   # every option in one file
archaemetry gui                      # the desktop app
archaemetry config                   # print the default config
```
`python run.py --config ... --input ... --outdir ...` still works.

Outputs, per sherd ID: `_drawing.svg/.pdf/.png` (1:1 mm), `_qa.png` (check plot: all
vertices, the chosen cut, thickness), `_report.json` (verdict, flags, axis, rim
diameter, thickness, config used), `_profile.json` (section and exterior in mm),
`_clean.ply`.

## Package layout

| Module | What it does |
|---|---|
| `meshio.py` | Load/save OBJ, PLY, STL, OFF, GLB (trimesh); connected components (scipy). |
| `cleanup.py` | Keep the largest piece (drops turntable and noise) or all pieces above a face count; optional capped hole fill (pymeshlab). |
| `axis.py` | `pca` (tries all three principal directions), `rim_arc` (open meshes only: needs the rim as a mesh border), `normals_refine` (surface-normal fit + profile-consistency refinement on the outer wall, multi-start ambiguity check; experimental), `manual`, `auto`. |
| `profile.py` | Rim orientation; robust band profile; **plane section** through the axis (most complete cut, checked against the robust exterior); split into outer/inner wall; **perpendicular wall thickness**, median over cuts. |
| `drawing.py` | Convention drawing: section filled (or hatched) on one side, exterior view mirrored on the other, axis line, rim line across the diameter, ID, block scale bar, caption. |
| `pipeline.py` | Config (defaults, validation, unknown-key warnings), scale, stages, verdict, outputs. |
| `selftest.py`, `synthetic.py` | Known geometry: closed two-walled sherds with break faces, a rounded lip and varying thickness; a complete two-walled vessel; exterior-only vessels. |
| `gui.py`, `cli.py` | Desktop app and command line. |

`tools/`: research harnesses (`axis_methods_eval.py`, `accuracy_oracle.py`,
`synthetic_validation.py`, `crop_wedge.py`, `diagnose_components.py`).
`packaging/`: PyInstaller spec, app entry point, field-test read-me and feedback form.

## Verdict

REVIEW is raised when: `rim_arc` was requested but could not fix the axis; the
surface-fit axis is ambiguous (near-optimal axes disagree by more than 0.3 mm) or
could not separate the walls; angular coverage is very low; cleanup kept under half
the mesh; or the ware is decorated (GigaMesh station). Thresholds are provisional,
set from synthetic tests, and must be rechecked on real sherds.

## Validation status (synthetic, exact ground truth)

Run `archaemetry selftest --full` to reproduce. Noise 0.05 mm on the surface.

| Case | Verdict | Section error (RMS) | Thickness error (mean) |
|---|---|---|---|
| Complete two-walled vessel | GO | 0.06 mm | 0.02 mm |
| Closed rim sherd, 120° | GO | 0.11 mm | 0.03 mm |
| Closed rim sherd, 90° | GO | 0.05 mm | 0.02 mm |
| Closed rim sherd, 60° | REVIEW | 0.49 mm | 0.03 mm |
| Closed rim sherd, 45° | REVIEW | 0.89 mm | 0.03 mm |
| Exterior-only vessel | GO | 0.06 mm (exterior) | n/a |

**Not yet validated:** real photographed sherds, real breaks, out-of-round hand-made
ware (a known 1-2 mm error floor that no fragment-only method can detect). See
`../docs/AXIS_PROTOTYPES.md`.
