# IMFPY User Guide

Welcome to the IMFPY Dust Trajectory Toolkit!  This user guide brings together everything you need to install the application, configure your environment, learn the workflow, and debug issues when things go wrong.

Use the table of contents below to jump straight to the topic you need:

- [Installation](installation.md)
- [Quick start tour](#quick-start-tour)
- [Tutorial library](tutorials.md)
- [Troubleshooting matrix](troubleshooting.md)
- [Frequently asked questions](faq.md)
- [Reference cheatsheet](#reference-cheatsheet)

---

## Quick start tour

### 1. Launch the application

```bash
python -m imfpy.gui.main
```

The first launch compiles the Fortran backend when a compiler is available.  A status dialog reports progress; allow several minutes on slower laptops.  When the window opens you will see three panes:

1. **Controls sidebar** – configure particle counts, integration parameters, and electromagnetic fields.
2. **Particle table** – edit initial state vectors, masses, charges, and β coefficients.
3. **3D viewer** – orbit plot that updates when a simulation finishes.

Open **Help → Quick Start Guide** (or press `F1`) from inside the app at any time to revisit this tour.

### 2. Run the default scenario

1. Leave the default four-particle circular orbit ensemble intact.
2. Click **Run Simulation**.
3. Watch the status bar: it changes from *Ready* → *Running simulation…* → *Simulation finished* when done.
4. Inspect the 3D viewer.  Drag with the left mouse button to orbit, use the scroll wheel to zoom, and right-click for the Matplotlib context menu.

### 3. Customise the scenario

1. Increase **Particles** to 16 (the table expands automatically).
2. Press **Populate Circular Orbit** to refresh the defaults, then edit a few rows to introduce eccentric or retrograde orbits.
3. Adjust `Δt` to 10.0 seconds and `Steps` to 10000 for higher resolution.
4. Modify the electric or magnetic field vectors to model different solar wind conditions.
5. Press **Run Simulation** again and compare the new curves.

### 4. Save results

- Use **File → Export Result…** (available after your first run) to persist the current simulation to an `.npz` archive.
- The CLI runner supports `--output path` for scripted pipelines (see [tutorials](tutorials.md#automating-with-the-cli)).

---

## Reference cheatsheet

| Task | Where to go |
| --- | --- |
| Install dependencies | [Installation guide](installation.md) |
| GPU setup | [Installation → Optional GPU support](installation.md#optional-gpu-support) |
| Understand every GUI control | [GUI deep dive](tutorials.md#gui-deep-dive) |
| Automate experiments | [Automating with the CLI](tutorials.md#automating-with-the-cli) |
| Resolve build/runtime errors | [Troubleshooting matrix](troubleshooting.md) |
| FAQs about physics and numerics | [Frequently asked questions](faq.md) |

Need something else?  File an issue on GitHub or contact the maintainers – we welcome feedback!
