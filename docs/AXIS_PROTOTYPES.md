# Axis prototypes: whole-surface fitting (experimental)

Status: **prototype, synthetic data only.** Default behaviour of the pipeline is unchanged;
opt in with `--axis-method normals_refine` (or `profile.axis_method: normals_refine`).
Code: [`cpu/archaemetry/axis.py`](../cpu/archaemetry/axis.py). Harness:
[`cpu/tools/axis_methods_eval.py`](../cpu/tools/axis_methods_eval.py).

## Why
`rim_arc` needs about 90 degrees of clean rim and refuses below that; body sherds have no
rim at all. Both cases leave a human in the loop. These two methods use the whole sherd
surface instead of the rim edge.

## Methods
1. **`normals`** - on a surface of revolution every surface-normal line meets the axis.
   The axis is the line closest to coplanar with all normal lines (line geometry, a 6x6
   eigenproblem, robust IRLS). Needs no rim.
2. **`refine`** - for the right axis, radius depends on height only, so vertices in a thin
   height band share a radius. Nudge tilt and offset of a seed axis to minimise that
   within-band spread (numpy pattern search, no new dependency).
3. **`fit_axis_multistart`** (what `normals_refine` runs) - refine from several seeds
   (normals, PCA, jittered normals) and keep the best. If several *near-optimal* axes
   disagree, the axis is ambiguous and the verdict becomes REVIEW.

## Results (synthetic vessels with known axis, 0.1 mm vertex noise)
Radial RMS in mm of the drawn profile against the truth.

| Case | pca | rim_arc | normals_refine |
|---|---|---|---|
| Lekythos, rim crop 60 deg | 1.47 | 3.94 (falls back) | **0.03** |
| Lekythos, rim crop 90 deg | 0.36 | 0.07 | 0.07 |
| Lekythos, body crop 150 deg (no rim) | 7.98 | 7.98 | **0.02** |
| Straight jar, body crop 90 deg | 17.2 | 17.2 | **0.03** |
| Cone bowl, rim crop 90 deg | 23.7 | 23.7 | 0.09 |
| Lekythos made 3% out-of-round, rim crop 90 deg | 0.41 | 1.90 | 0.87 |

Through the whole lane (`run.py` to published SVG): a 60 deg rim sherd goes from
1.39 mm RMS (`rim_arc`, REVIEW) to 0.027 mm; 60-120 deg body sherds give 0.01-0.03 mm; a
45 deg rim sherd gives 0.41 mm and is flagged REVIEW.

## Limits found (read before trusting)
- **Synthetic only.** Perfect surfaces of revolution plus Gaussian noise. Not yet run on a
  real mesh or a real break.
- **Out-of-round vessels set an error floor that no fragment-only method can see.** A short
  wedge of a 3% oval vessel is a perfectly consistent surface of revolution about a
  slightly different axis. All methods are off by roughly 1-2.5 mm there, and neither the
  fit spread nor the ambiguity check notices. Only validation against a calliper on real
  sherds can size this.
- **Spheres and cones have no unique axis.** Multi-start disagreement flags most of these
  cases, with some false alarms (a cone rim at 90 deg is flagged but was accurate).
- **Known miss:** a hemisphere-pole crop converged to a wrong axis with no disagreement
  between seeds (13 mm error, not flagged). Globular jars are the risk.
- **Gate thresholds are provisional.** `AMBIGUITY_MM_MAX = 0.3` in `cpu/run.py` and the
  relaxed coverage floor (0.10 instead of 0.25) were set from these synthetic runs.
- A noise-based uncertainty estimate (`axis_uncertainty`) was prototyped and is **not**
  used as a gate: it missed 4 of 6 sphere failures and ignores model error.

## Reproduce
```bash
cd cpu
python tools/axis_methods_eval.py out/ --noise 0.1 --degrees 60,90,150 --methods pca,rim_arc,normals+refine
python tools/synthetic_validation.py out2/      # includes the E_/F_ normals_refine cases
```

## Update: closed two-walled sherds (the real case)

Everything above was measured on **open, exterior-only** meshes (wedges cut from a
complete vessel). A photogrammetry scan of a whole sherd is a **closed** surface: outer
wall, inner wall and break faces. Testing on closed synthetic sherds
(`archaemetry.synthetic.make_sherd`) exposed three problems, now fixed:

1. **`rim_arc` cannot work on a closed mesh.** It finds the rim from the mesh's open
   border, and a closed sherd has none. `auto` therefore uses the surface fit for real
   sherds; `rim_arc` remains for open scans.
2. **The consistency refinement read wall thickness as radius spread.** With both walls
   present, each height band holds two radii. The refinement now uses only the outer
   wall (wall vertices split by the sign of normal . radial; the group at larger radius
   is the outer wall), and the normal-line fit uses only wall vertices (break-face
   normals do not meet the axis).
3. **PCA picked the diameter of squat vessels as the axis.** It now refines from all
   three principal directions and keeps the most rotationally consistent one.

Closed rim sherd (jar profile, 0.1 mm noise), radial RMS against the true axis:

| Arc | Before fix | After fix | Verdict |
|---|---|---|---|
| 120 deg | 3.1 mm spread, ambiguous | 0.10 mm | GO |
| 90 deg | - | 0.20 mm | REVIEW (ambiguity 0.49 mm, conservative) |
| 60 deg | - | 0.67 mm | REVIEW |
| 45 deg | - | 13 mm | REVIEW (wall classification failed: new flag) |

The 45 deg case shows why the new flag exists: with too little surface the outer/inner
split collapses (8% "outer" instead of ~50%) and the axis is wrong without any
disagreement between seeds; the run is now flagged instead of passing.
