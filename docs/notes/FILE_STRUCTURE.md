# Repository organization

The project assets have been grouped by purpose to make navigation easier:

- `src/imfpy/` – Python source files that implement simulations, visualisations and utilities.
- `src/fortran/` – FORTRAN source files used by the collisional integrator workflow.
- `scripts/` – Helper scripts for database conversion and other command line utilities.
- `data/raw/` – Simulation inputs and representative outputs (e.g. `subtest.8`, `in.dat`).
- `data/sql/` – SQLite databases and SQL exports generated from charge balance studies.
- `docs/papers/` – Reference papers that inform the modelling work.
- `docs/presentations/` – Presentation slide decks and related documents.
- `docs/notes/` – Miscellaneous documentation, including this overview.
- `media/figures/` – Static images and figure exports from the simulations.
- `media/videos/` – Movies and animations of particle trajectories and magnetic fields.
- `notebooks/` – Mathematica notebooks capturing exploratory field calculations.

This layout mirrors the original numbering scheme while removing duplicated placeholder files and
gathering related resources into dedicated folders.
