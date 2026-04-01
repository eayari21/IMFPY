#!/usr/bin/env python3
"""
Publication-quality field-line plots for the heliospheric electric and magnetic fields.

This version fixes the geometric inconsistency:

- The heliopause nose / parabolic cavity appears only in planes that contain the z-axis.
- Therefore, to match the flow-line cavity visually, use x-z and y-z cuts,
  not x-y at z = 0.

Fields used:
- Electric field from Eq. (13)
- Magnetic field from Eq. (27)

Reference:
ISD_Filtering_GRL-2.pdf, Sections 2.2 and 2.3
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from matplotlib.patches import FancyArrowPatch


# =============================================================================
# Style
# =============================================================================
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 15,
    "axes.labelsize": 19,
    "axes.titlesize": 21,
    "axes.linewidth": 1.4,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "xtick.major.size": 7,
    "ytick.major.size": 7,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.width": 1.0,
    "ytick.minor.width": 1.0,
    "legend.frameon": False,
})


# =============================================================================
# Parameters
# =============================================================================
@dataclass
class Params:
    z0: float = 100.0
    B0: float = 1.0
    v0: float = 1.0
    c: float = 1.0

    eps: float = 1e-10
    ds: float = 1.0
    max_steps: int = 7000
    min_speed: float = 1e-10

    # Plot windows
    x_lim: float = 200.0
    y_lim: float = 200.0
    z_min: float = -220.0
    z_max: float = 140.0

    # Seed density
    n_edge: int = 34
    n_halo: int = 34

    # Arrow placement
    arrows_per_line: int = 2
    arrow_scale: int = 11


P = Params()


# =============================================================================
# Geometry
# =============================================================================
def R_val(x, y, z):
    return np.sqrt(x*x + y*y + z*z)


def r_cyl(x, y):
    return np.sqrt(x*x + y*y)


def S2_val(x, y, z, z0):
    """
    Eq. (9): S^2 = r^2 + 2 z0^2 (z/R - 1)
    """
    R = np.maximum(R_val(x, y, z), P.eps)
    r = r_cyl(x, y)
    return r*r + 2.0 * z0*z0 * (z / R - 1.0)


def outside_heliopause(x, y, z, z0):
    return S2_val(x, y, z, z0) > 0.0


# =============================================================================
# Fields from the manuscript
# =============================================================================
def electric_field(x, y, z, p: Params):
    """
    Electric field from Eq. (13), using the compact helper from Eq. (14).
    """
    z0, B0, v0, c, eps = p.z0, p.B0, p.v0, p.c, p.eps

    r = np.maximum(r_cyl(x, y), eps)
    R = np.maximum(R_val(x, y, z), eps)
    S2raw = S2_val(x, y, z, z0)

    valid = S2raw > 0.0
    S = np.sqrt(np.maximum(S2raw, eps))

    helper = (-S / (r*r)) + (1.0 / S) - (z0*z0 * z / (S * R**3))
    pref = B0 * v0 / (c * r)

    Ex = pref * (helper * x * y)
    Ey = pref * (S + helper * y*y)
    Ez = pref * ((z0*z0 / (r*r * R**3 * S)) * y)

    Ex = np.where(valid, Ex, np.nan)
    Ey = np.where(valid, Ey, np.nan)
    Ez = np.where(valid, Ez, np.nan)

    return Ex, Ey, Ez


def magnetic_field_eq27(x, y, z, p: Params):
    """
    Magnetic field from Eq. (27), which is the Bz = 0 case.
    """
    z0, B0, eps = p.z0, p.B0, p.eps

    r = np.maximum(r_cyl(x, y), eps)
    R = np.maximum(R_val(x, y, z), eps)
    S2raw = S2_val(x, y, z, z0)

    valid = S2raw > 0.0
    S = np.sqrt(np.maximum(S2raw, eps))

    A = -1.0 + (z0*z0 * z / R**3)
    T = (S / (r*r)) - (1.0 / S) + (z0*z0 * z / (S * R**3))
    U = (-S / (r*r)) + (1.0 / S) - (z0*z0 * z / (S * R**3))

    pref = (B0 / r) * A

    Bx = pref * (-S + T * y*y)
    By = pref * (U * x * y)
    Bz = np.zeros_like(Bx)

    Bx = np.where(valid, Bx, np.nan)
    By = np.where(valid, By, np.nan)
    Bz = np.where(valid, Bz, np.nan)

    return Bx, By, Bz


# =============================================================================
# Plane-restricted projected fields
# =============================================================================
def projected_field_xz(x, z, which="E", p: Params = P):
    """
    x-z cut at y = 0
    """
    y = 0.0
    if which == "E":
        Ex, _, Ez = electric_field(x, y, z, p)
        return np.array([Ex, Ez], dtype=float)
    elif which == "B":
        Bx, _, Bz = magnetic_field_eq27(x, y, z, p)
        return np.array([Bx, Bz], dtype=float)
    raise ValueError("which must be 'E' or 'B'")


def projected_field_yz(y, z, which="E", p: Params = P):
    """
    y-z cut at x = 0
    """
    x = 0.0
    if which == "E":
        _, Ey, Ez = electric_field(x, y, z, p)
        return np.array([Ey, Ez], dtype=float)
    elif which == "B":
        _, By, Bz = magnetic_field_eq27(x, y, z, p)
        return np.array([By, Bz], dtype=float)
    raise ValueError("which must be 'E' or 'B'")


# =============================================================================
# RK4 field-line tracer
# =============================================================================
def rk4_step_2d(pos, h, field_func, normalize=True):
    def tangent(q):
        F = np.asarray(field_func(q[0], q[1]), dtype=float)
        if np.any(~np.isfinite(F)):
            return np.array([np.nan, np.nan], dtype=float)

        mag = np.hypot(F[0], F[1])
        if mag < P.min_speed:
            return np.array([np.nan, np.nan], dtype=float)

        return F / mag if normalize else F

    k1 = tangent(pos)
    if np.any(~np.isfinite(k1)):
        return np.array([np.nan, np.nan])

    k2 = tangent(pos + 0.5 * h * k1)
    if np.any(~np.isfinite(k2)):
        return np.array([np.nan, np.nan])

    k3 = tangent(pos + 0.5 * h * k2)
    if np.any(~np.isfinite(k3)):
        return np.array([np.nan, np.nan])

    k4 = tangent(pos + h * k3)
    if np.any(~np.isfinite(k4)):
        return np.array([np.nan, np.nan])

    return pos + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def trace_field_line(seed, field_func, inside_domain, ds, max_steps):
    def march(seed_pos, h):
        pts = [np.asarray(seed_pos, dtype=float)]
        pos = np.asarray(seed_pos, dtype=float)

        for _ in range(max_steps):
            new_pos = rk4_step_2d(pos, h, field_func, normalize=True)

            if np.any(~np.isfinite(new_pos)):
                break
            if not inside_domain(new_pos[0], new_pos[1]):
                break

            pts.append(new_pos.copy())
            pos = new_pos

        return np.array(pts)

    forward = march(seed, +ds)
    backward = march(seed, -ds)

    if len(backward) > 1:
        backward = backward[::-1][:-1]

    return np.vstack([backward, forward])


# =============================================================================
# Domain tests
# =============================================================================
def inside_xz(x, z, p: Params):
    if x < -p.x_lim or x > p.x_lim or z < p.z_min or z > p.z_max:
        return False
    return outside_heliopause(x, 0.0, z, p.z0)


def inside_yz(y, z, p: Params):
    if y < -p.y_lim or y > p.y_lim or z < p.z_min or z > p.z_max:
        return False
    return outside_heliopause(0.0, y, z, p.z0)


# =============================================================================
# Seed generation
# =============================================================================
def make_edge_seeds(xmin, xmax, zmin, zmax, n):
    xs = np.linspace(xmin, xmax, n)
    zs = np.linspace(zmin, zmax, n)

    seeds = []
    for z in zs:
        seeds.append((xmin, z))
        seeds.append((xmax, z))
    for x in xs:
        seeds.append((x, zmin))
        seeds.append((x, zmax))

    return np.array(seeds, dtype=float)


def make_halo_seeds_meridional(z0, n):
    """
    Seeds just outside the heliopause in a meridional plane.
    Heliospause boundary in nondimensional form:
        s = (2 - rho^2) / sqrt(4 - rho^2)
    """
    rho = np.linspace(0.0, 1.96, n // 2)
    s = (2.0 - rho**2) / np.sqrt(np.maximum(4.0 - rho**2, 1e-12))

    a = 1.03
    x = a * z0 * rho
    z = a * z0 * s

    upper = np.column_stack([+x, z])
    lower = np.column_stack([-x, z])

    return np.vstack([upper, lower])


# =============================================================================
# Plot helpers
# =============================================================================
def add_heliopause_meridional(ax, horiz_label="x", z0=100.0, color="white"):
    q = np.linspace(-1.999, 1.999, 1400)
    z = z0 * (2.0 - q*q) / np.sqrt(np.maximum(4.0 - q*q, 1e-12))
    h = z0 * q

    z_floor = -260.0
    h_fill = np.concatenate([h, h[::-1]])
    z_fill = np.concatenate([z, np.full_like(z, z_floor)])

    ax.fill(h_fill, z_fill, color=color, zorder=6)
    ax.plot(h, z, color="black", lw=1.8, zorder=7)


def format_axes(ax, xlabel, ylabel, xlim, ylim, equal=False):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.minorticks_on()
    ax.tick_params(direction="in", which="both", top=True, right=True)
    if equal:
        ax.set_aspect("equal", adjustable="box")


def add_arrows_to_line(ax, line, color, n_arrows=2, mutation_scale=10, zorder=5):
    """
    Add small arrowheads along an already-integrated polyline.
    """
    if len(line) < 8:
        return

    seg = np.diff(line, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    total = np.sum(seglen)
    if total <= 0:
        return

    cum = np.concatenate([[0.0], np.cumsum(seglen)])
    targets = np.linspace(0.28 * total, 0.72 * total, n_arrows)

    for t in targets:
        idx = np.searchsorted(cum, t) - 1
        idx = np.clip(idx, 0, len(line) - 2)

        p0 = line[idx]
        p1 = line[idx + 1]

        if np.allclose(p0, p1):
            continue

        arrow = FancyArrowPatch(
            posA=(p0[0], p0[1]),
            posB=(p1[0], p1[1]),
            arrowstyle='->',
            mutation_scale=mutation_scale,
            lw=0.0,
            color=color,
            zorder=zorder
        )
        ax.add_patch(arrow)


def plot_lines(ax, lines, color, lw=1.35, alpha=0.9, arrows_per_line=2, arrow_scale=10):
    for line in lines:
        if len(line) < 2:
            continue
        ax.plot(line[:, 0], line[:, 1], color=color, lw=lw, alpha=alpha, zorder=2)
        add_arrows_to_line(
            ax,
            line,
            color=color,
            n_arrows=arrows_per_line,
            mutation_scale=arrow_scale,
            zorder=4
        )


# =============================================================================
# Build one panel
# =============================================================================
def build_meridional_panel(ax, plane="xz", p: Params = P):
    if plane == "xz":
        horiz_min, horiz_max = -p.x_lim, p.x_lim
        seeds_edge = make_edge_seeds(horiz_min, horiz_max, p.z_min, p.z_max, p.n_edge)
        seeds_halo = make_halo_seeds_meridional(p.z0, p.n_halo)

        def dom(a, z):
            return inside_xz(a, z, p)

        Efunc = lambda a, z: projected_field_xz(a, z, which="E", p=p)
        Bfunc = lambda a, z: projected_field_xz(a, z, which="B", p=p)
        xlabel = r"$x\ \mathrm{[AU]}$"
        title = r"$x$-$z$ cut at $y=0$"

    elif plane == "yz":
        horiz_min, horiz_max = -p.y_lim, p.y_lim
        seeds_edge = make_edge_seeds(horiz_min, horiz_max, p.z_min, p.z_max, p.n_edge)
        seeds_halo = make_halo_seeds_meridional(p.z0, p.n_halo)

        def dom(a, z):
            return inside_yz(a, z, p)

        Efunc = lambda a, z: projected_field_yz(a, z, which="E", p=p)
        Bfunc = lambda a, z: projected_field_yz(a, z, which="B", p=p)
        xlabel = r"$y\ \mathrm{[AU]}$"
        title = r"$y$-$z$ cut at $x=0$"

    else:
        raise ValueError("plane must be 'xz' or 'yz'")

    E_lines = []
    B_lines = []

    for seed in np.vstack([seeds_edge, seeds_halo]):
        if dom(seed[0], seed[1]):
            line = trace_field_line(seed, Efunc, dom, p.ds, p.max_steps)
            if len(line) > 12:
                E_lines.append(line)

    for seed in np.vstack([seeds_edge, seeds_halo]):
        if dom(seed[0], seed[1]):
            line = trace_field_line(seed, Bfunc, dom, p.ds, p.max_steps)
            if len(line) > 12:
                B_lines.append(line)

    plot_lines(
        ax, E_lines, color="red",
        lw=1.25, alpha=0.85,
        arrows_per_line=p.arrows_per_line,
        arrow_scale=p.arrow_scale
    )
    plot_lines(
        ax, B_lines, color="blue",
        lw=1.25, alpha=0.85,
        arrows_per_line=p.arrows_per_line,
        arrow_scale=p.arrow_scale
    )

    add_heliopause_meridional(ax, z0=p.z0, color="white")

    format_axes(
        ax,
        xlabel=xlabel,
        ylabel=r"$z\ \mathrm{[AU]}$",
        xlim=(horiz_min, horiz_max),
        ylim=(p.z_min, p.z_max),
        equal=False
    )

    ax.set_title(title)
    ax.text(
        0.02, 0.98,
        "red: electric field lines\nblue: magnetic field lines",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=13,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.75", alpha=0.95),
    )


# =============================================================================
# Main
# =============================================================================
def main():
    # Two publication-quality panels
    fig1, ax1 = plt.subplots(figsize=(9.0, 7.2))
    build_meridional_panel(ax1, plane="xz", p=P)
    fig1.tight_layout()
    fig1.savefig("field_lines_xz.pdf")
    fig1.savefig("field_lines_xz.png")

    fig2, ax2 = plt.subplots(figsize=(9.0, 7.2))
    build_meridional_panel(ax2, plane="yz", p=P)
    fig2.tight_layout()
    fig2.savefig("field_lines_yz.pdf")
    fig2.savefig("field_lines_yz.png")

    plt.show()


if __name__ == "__main__":
    main()
