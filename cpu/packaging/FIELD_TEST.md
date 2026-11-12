# Archaemetry field test: read me first

Thank you for testing. Archaemetry turns a 3D mesh of a sherd or vessel into a profile
drawing (section filled black, axis, rim line, exterior view, scale bar) with a
GO / REVIEW verdict. It runs on your own computer; nothing is uploaded.

## 1. Install (no Python needed)
- **Windows:** unzip `Archaemetry-windows.zip` anywhere, open the folder, double-click
  `Archaemetry.exe`. If Windows SmartScreen warns about an unknown publisher, choose
  *More info -> Run anyway* (the test build is not code-signed).
- **macOS:** unzip, then **right-click** `Archaemetry.app` -> *Open* -> *Open* (first time
  only; the test build is not notarised). If macOS still refuses, run
  `xattr -dr com.apple.quarantine Archaemetry.app` in Terminal.
- **Linux:** extract `Archaemetry-linux.tar.gz` and run `Archaemetry/Archaemetry`.

## 2. First: run the self-test
Click **Self-test**. It draws known shapes and checks the result against the exact
answer. It should say *Self-test passed*. If it fails, please send us the Log tab.

## 3. Make a drawing
1. **Browse** to a mesh (OBJ, PLY, STL). It should be scaled; if it is in cm or m,
   pick that under *Mesh units*.
2. Type the sherd ID.
3. Say what it is: *rim sherd*, *body sherd*, *complete vessel*, or *already oriented*
   (if you rotated it in MeshLab so the vessel axis is the Y axis).
4. Click **Make drawing**. A 1-2 million triangle mesh takes about a minute.

Outputs (in `archaemetry_out` next to the mesh): `<ID>_drawing.svg / .pdf / .png` at
1:1 scale, `<ID>_qa.png` (check plot), `<ID>_report.json`, `<ID>_profile.json`.

## 4. Read the verdict
- **GO**: the axis was fixed from good evidence. Still look at the drawing.
- **REVIEW**: the program is not confident (too little rim, an ambiguous shape,
  decorated ware...). The reason is listed. Do not use the drawing without checking it.
- If the rim is at the bottom, tick *Rim came out at the bottom: flip* and run again.

## 5. What we need back from you
For each sherd you try, please send:
1. the `<ID>_report.json` file (no images of your data are needed unless you agree),
2. your calliper measurements: rim diameter if known, and wall thickness at 3-5 points
   with their depth below the rim,
3. the feedback form (`FEEDBACK_FORM.md`).

Known limits of this test build: wall thickness needs both the inside and outside of
the sherd in the mesh; the surface-fit axis used for small and body sherds is
experimental; out-of-round hand-made vessels can be off by 1-2 mm without warning.
