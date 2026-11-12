"""Desktop app for field testers: pick a mesh, answer three questions, get the drawing.

Plain tkinter (ships with Python and freezes cleanly with PyInstaller). The pipeline
runs on a worker thread; results come back through a queue polled by the UI thread.
"""
import os
import queue
import subprocess
import sys
import threading
import traceback
import webbrowser

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from . import __version__

GUIDE_URL = "https://amudhans-05.github.io/Archaemetry/field-guide.html"
SHERD_TYPES = [
    ("Rim sherd", "auto"),
    ("Body sherd (no rim)", "normals_refine"),
    ("Complete or near-complete vessel", "pca"),
    ("Already oriented: axis is the Y axis", "manual"),
]
UNITS = [("millimetres", 1.0), ("centimetres", 10.0), ("metres", 1000.0)]
COLORS = {"GO": "#3f7a46", "REVIEW": "#b36b00", "ERROR": "#a33", "": "#555"}


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"Archaemetry {__version__} - pottery profile drawing")
        self.geometry("1100x760")
        self.minsize(900, 600)
        self.q = queue.Queue()
        self.busy = False
        self.last = None
        self._img = None
        self._build()
        self.after(100, self._poll)

    # ------------------------------------------------------------------ layout
    def _build(self):
        pad = dict(padx=8, pady=4)
        left = ttk.Frame(self, padding=10)
        left.pack(side="left", fill="y")
        right = ttk.Frame(self, padding=(0, 10, 10, 10))
        right.pack(side="left", fill="both", expand=True)

        ttk.Label(left, text="1. Mesh file", font=("", 10, "bold")).pack(anchor="w")
        f = ttk.Frame(left); f.pack(fill="x", **pad)
        self.mesh = tk.StringVar()
        ttk.Entry(f, textvariable=self.mesh, width=34).pack(side="left", fill="x", expand=True)
        ttk.Button(f, text="Browse...", command=self._browse).pack(side="left", padx=4)

        ttk.Label(left, text="2. Sherd ID", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.sid = tk.StringVar()
        ttk.Entry(left, textvariable=self.sid, width=40).pack(fill="x", **pad)

        ttk.Label(left, text="3. What is it?", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.kind = tk.StringVar(value=SHERD_TYPES[0][1])
        for label, val in SHERD_TYPES:
            ttk.Radiobutton(left, text=label, value=val, variable=self.kind).pack(anchor="w", padx=12)

        ttk.Label(left, text="4. Mesh units", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.units = tk.StringVar(value=UNITS[0][0])
        ttk.Combobox(left, textvariable=self.units, values=[u for u, _ in UNITS],
                     state="readonly", width=16).pack(anchor="w", **pad)

        opt = ttk.LabelFrame(left, text="Options", padding=6)
        opt.pack(fill="x", pady=(10, 0))
        self.flip = tk.BooleanVar(); self.deco = tk.BooleanVar(); self.right_side = tk.BooleanVar()
        self.hatch = tk.BooleanVar()
        ttk.Checkbutton(opt, text="Rim came out at the bottom: flip", variable=self.flip).pack(anchor="w")
        ttk.Checkbutton(opt, text="Section on the right", variable=self.right_side).pack(anchor="w")
        ttk.Checkbutton(opt, text="Hatch the section instead of solid", variable=self.hatch).pack(anchor="w")
        ttk.Checkbutton(opt, text="Decorated ware", variable=self.deco).pack(anchor="w")

        ttk.Label(left, text="Output folder", font=("", 10, "bold")).pack(anchor="w", pady=(10, 0))
        f = ttk.Frame(left); f.pack(fill="x", **pad)
        self.outdir = tk.StringVar()
        ttk.Entry(f, textvariable=self.outdir, width=34).pack(side="left", fill="x", expand=True)
        ttk.Button(f, text="...", width=3, command=self._browse_out).pack(side="left", padx=4)

        self.run_btn = ttk.Button(left, text="Make drawing", command=self._run)
        self.run_btn.pack(fill="x", pady=(12, 4))
        b = ttk.Frame(left); b.pack(fill="x")
        ttk.Button(b, text="Open output folder", command=self._open_out).pack(side="left", fill="x", expand=True)
        ttk.Button(b, text="Self-test", command=self._selftest).pack(side="left", padx=4)
        ttk.Button(b, text="Field guide", command=lambda: webbrowser.open(GUIDE_URL)).pack(side="left")
        self.pb = ttk.Progressbar(left, mode="determinate")
        self.pb.pack(fill="x", pady=8)

        self.banner = tk.Label(right, text="Choose a mesh to begin.", anchor="w", fg="white",
                               bg=COLORS[""], font=("", 12, "bold"), padx=10, pady=6)
        self.banner.pack(fill="x")
        self.detail = tk.Label(right, text="", anchor="w", justify="left", wraplength=700, padx=4)
        self.detail.pack(fill="x", pady=(4, 6))
        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True)
        self.canvas = tk.Canvas(nb, bg="white", highlightthickness=0)
        self.canvas.bind("<Configure>", lambda e: self._show_image())
        self.log = tk.Text(nb, height=10, wrap="word", font=("Courier", 9))
        nb.add(self.canvas, text="Drawing")
        nb.add(self.log, text="Log")
        self.view = tk.StringVar(value="drawing")
        r = ttk.Frame(right); r.pack(fill="x", pady=(4, 0))
        ttk.Radiobutton(r, text="Drawing", value="drawing", variable=self.view, command=self._show_image).pack(side="left")
        ttk.Radiobutton(r, text="Check plot (all points + cut + thickness)", value="qa", variable=self.view,
                        command=self._show_image).pack(side="left", padx=8)

    # ------------------------------------------------------------------ actions
    def _browse(self):
        p = filedialog.askopenfilename(title="Choose a mesh", filetypes=[
            ("3D meshes", "*.obj *.ply *.stl *.off *.glb"), ("All files", "*.*")])
        if p:
            self.mesh.set(p)
            stem = os.path.splitext(os.path.basename(p))[0]
            if not self.sid.get():
                self.sid.set(stem)
            if not self.outdir.get():
                self.outdir.set(os.path.join(os.path.dirname(p), "archaemetry_out"))

    def _browse_out(self):
        p = filedialog.askdirectory(title="Output folder")
        if p:
            self.outdir.set(p)

    def _config(self):
        k = self.kind.get()
        return {
            "sherd": {"id": self.sid.get().strip() or None},
            "scale": {"mm_per_unit": dict(UNITS)[self.units.get()]},
            "profile": {"axis_method": k, "orientation": "flip" if self.flip.get() else "auto"},
            "drawing": {"section_side": "right" if self.right_side.get() else "left",
                        "fill": "hatch" if self.hatch.get() else "solid"},
            "ware": "decorated" if self.deco.get() else "plain",
        }

    def _run(self):
        if self.busy:
            return
        mesh = self.mesh.get()
        if not os.path.isfile(mesh):
            messagebox.showwarning("Archaemetry", "Choose a mesh file first.")
            return
        out = self.outdir.get() or os.path.join(os.path.dirname(mesh), "archaemetry_out")
        self.outdir.set(out)
        cfg = self._config()
        self._start(f"Working on {os.path.basename(mesh)}...")
        threading.Thread(target=self._work, args=(mesh, out, cfg), daemon=True).start()

    def _work(self, mesh, out, cfg):
        from .pipeline import run
        try:
            rep = run(mesh, out, cfg, log=lambda m: self.q.put(("log", m)))
            self.q.put(("done", (out, rep)))
        except Exception as e:  # report, do not crash the app
            self.q.put(("log", traceback.format_exc()))
            self.q.put(("error", f"{type(e).__name__}: {e}"))

    def _selftest(self):
        if self.busy:
            return
        self._start("Running self-test on known shapes (about a minute)...")

        def work():
            from .selftest import main
            try:
                ok, res = main(log=lambda m: self.q.put(("log", m)))
                self.q.put(("selftest", (ok, res)))
            except Exception as e:
                self.q.put(("log", traceback.format_exc()))
                self.q.put(("error", f"{type(e).__name__}: {e}"))
        threading.Thread(target=work, daemon=True).start()

    def _open_out(self):
        p = self.outdir.get()
        if not p or not os.path.isdir(p):
            return
        if sys.platform.startswith("win"):
            os.startfile(p)  # noqa: S606 - opens the user's own folder
        else:
            subprocess.Popen(["open" if sys.platform == "darwin" else "xdg-open", p])

    # ------------------------------------------------------------------ ui updates
    def _start(self, msg):
        self.busy = True
        self.run_btn.state(["disabled"])
        self.pb.config(mode="indeterminate")
        self.pb.start(12)
        self._banner("", msg)
        self.detail.config(text="")
        self.log.delete("1.0", "end")

    def _finish(self):
        self.busy = False
        self.run_btn.state(["!disabled"])
        self.pb.stop()
        self.pb.config(mode="determinate")
        self.pb["value"] = 0

    def _banner(self, kind, text):
        self.banner.config(text=text, bg=COLORS.get(kind, COLORS[""]))

    def _poll(self):
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind == "log":
                    self.log.insert("end", payload + "\n"); self.log.see("end")
                elif kind == "done":
                    self._finish()
                    out, rep = payload
                    self.last = (out, rep)
                    v = rep["verdict"]
                    self._banner(v, "GO - drawing is ready" if v == "GO"
                                 else "REVIEW - check before using this drawing")
                    lines = [f"! {f}" for f in rep["flags"]] + [f"- {n}" for n in rep["notes"]]
                    p = rep["profile"]
                    t = rep["wall_thickness_mm"]
                    lines.append(f"Rim diameter {p['rim_diameter_mm']} mm, height {p['height_mm']} mm"
                                 + (f", wall {t['min']}-{t['max']} mm" if t else ", no wall thickness (one surface)"))
                    lines.append(f"Saved to {out}")
                    self.detail.config(text="\n".join(lines))
                    self._show_image()
                elif kind == "selftest":
                    self._finish()
                    ok, res = payload
                    self._banner("GO" if ok else "ERROR", "Self-test passed: this install measures correctly"
                                 if ok else "Self-test FAILED: please send the log to the project")
                    self.detail.config(text="\n".join(
                        f"{r['case']}: {r['verdict']}  section error {r.get('section_rms', '-')} mm, "
                        f"thickness error {r.get('thickness_mean_err', '-')} mm" for r in res))
                elif kind == "error":
                    self._finish()
                    self._banner("ERROR", "Something went wrong")
                    self.detail.config(text=payload + "\nSee the Log tab.")
        except queue.Empty:
            pass
        self.after(100, self._poll)

    def _show_image(self):
        self.canvas.delete("all")
        if not self.last:
            return
        out, rep = self.last
        key = "png" if self.view.get() == "drawing" else "qa"
        name = rep["outputs"].get(key)
        if not name:
            return
        path = os.path.join(out, name)
        try:
            img = tk.PhotoImage(file=path)
        except tk.TclError:
            return
        cw, ch = max(self.canvas.winfo_width(), 50), max(self.canvas.winfo_height(), 50)
        k = max(1, -(-img.width() // cw), -(-img.height() // ch))
        self._img = img.subsample(k, k)
        self.canvas.create_image(cw // 2, ch // 2, image=self._img)


def main():
    App().mainloop()


if __name__ == "__main__":
    main()
