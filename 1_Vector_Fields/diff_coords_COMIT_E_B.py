#!/usr/bin/env python3
"""
Field-line plots for the heliospheric electric and magnetic fields.

Display-frame mapping
---------------------
We keep the manuscript equations in their native/model coordinates
(xm, ym, zm), where the cavity axis is along +zm.

We display the fields in the plotting frame

    x_plot =  xm
    y_plot =  zm
    z_plot = -ym

so that

- the y-z slice at x = 0 is the side view with the parabolic cavity
- the x-z slice at y = 0 is the nose/cross-section view
- B_infinity is out of the side-view plane in the +x direction

Produced figures
----------------
1) Side view: y-z plane at x = 0
2) Nose view: x-z plane at y = 0
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from dataclasses import dataclass


# =============================================================================
# Plot styling
# =============================================================================
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 17,
    "axes.linewidth": 1.2,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "xtick.minor.width": 0.9,
    "ytick.minor.width": 0.9,
    "legend.frameon": False,
})


# =============================================================================
# User-tunable parameters
# =============================================================================
@dataclass
class Params:
    z0: float = 100.0   # AU
    nose_radius: float = 50.0   # AU   # try 45–55 to match Proper_E

    B0: float = 1.0
    v0: float = 1.0
    c: float = 1.0

    eps: float = 1e-10

    ds: float = 1.05
    max_steps: int = 7000
    min_speed: float = 1e-10

    # Side view: y-z plane at x = 0
    yz_y_min: float = -220.0
    yz_y_max: float = 140.0
    yz_z_lim: float = 200.0

    # Nose view: x-z plane at y = 0
    xz_lim: float = 200.0

    n_edge_seeds_side: int = 12
    n_edge_seeds_nose: int = 7
    n_halo_seeds_side: int = 18
    n_halo_seeds_nose: int = 12

    arrow_size: int = 9


P = Params()


# =============================================================================
# Core geometry in the original/model frame
# =============================================================================
def R_val(x, y, z):
    return np.sqrt(x*x + y*y + z*z)


def r_cyl(x, y):
    return np.sqrt(x*x + y*y)


def S2_val(x, y, z, z0):
    """
    S^2 = r^2 + 2 z0^2 (z/R - 1)
    """
    r = r_cyl(x, y)
    R = np.maximum(R_val(x, y, z), P.eps)
    return r*r + 2.0 * z0*z0 * (z / R - 1.0)


def outside_heliopause_model(x, y, z, z0):
    return S2_val(x, y, z, z0) > 0.0


# =============================================================================
# Fields in the original/model frame
# =============================================================================
def electric_field_model(x, y, z, p: Params):
    """
    Electric field from Eq. (13), evaluated in the model frame.
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


def magnetic_field_model(x, y, z, p: Params):
    """
    Magnetic field from Eq. (27), evaluated in the model frame.
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

    pref = B0 / r * A

    Bx = pref * (-S + T * y*y)
    By = pref * (U * x * y)
    Bz = np.zeros_like(Bx)

    Bx = np.where(valid, Bx, np.nan)
    By = np.where(valid, By, np.nan)
    Bz = np.where(valid, Bz, np.nan)

    return Bx, By, Bz


# =============================================================================
# Coordinate transform: model frame -> plotting frame
# =============================================================================
def plot_to_model(xp, yp, zp):
    """
    Plotting frame:
        x_plot =  x_model
        y_plot =  z_model
        z_plot = -y_model

    Therefore:
        x_model =  x_plot
        y_model = -z_plot
        z_model =  y_plot
    """
    xm = xp
    ym = -zp
    zm = yp
    return xm, ym, zm


def vector_model_to_plot(Fxm, Fym, Fzm):
    """
    Vector components under the same transform:
        Fx_plot =  Fx_model
        Fy_plot =  Fz_model
        Fz_plot = -Fy_model
    """
    Fxp = Fxm
    Fyp = Fzm
    Fzp = -Fym
    return Fxp, Fyp, Fzp


def electric_field_plot(xp, yp, zp, p: Params):
    xm, ym, zm = plot_to_model(xp, yp, zp)
    Exm, Eym, Ezm = electric_field_model(xm, ym, zm, p)
    return vector_model_to_plot(Exm, Eym, Ezm)


def magnetic_field_plot(xp, yp, zp, p: Params):
    xm, ym, zm = plot_to_model(xp, yp, zp)
    Bxm, Bym, Bzm = magnetic_field_model(xm, ym, zm, p)
    return vector_model_to_plot(Bxm, Bym, Bzm)


def outside_heliopause_plot(xp, yp, zp, z0):
    xm, ym, zm = plot_to_model(xp, yp, zp)
    return outside_heliopause_model(xm, ym, zm, z0)


# =============================================================================
# 2D projected fields in plotting coordinates
# =============================================================================
def projected_field_side_yz(y, z, which="E", p: Params = P):
    """
    Side view: y-z plane at x = 0
    Horizontal axis = y
    Vertical axis   = z
    """
    x = 0.0
    if which == "E":
        _, Fy, Fz = electric_field_plot(x, y, z, p)
    elif which == "B":
        _, Fy, Fz = magnetic_field_plot(x, y, z, p)
    else:
        raise ValueError("which must be 'E' or 'B'")
    return np.array([Fy, Fz], dtype=float)


def projected_field_nose_xz(x, z, which="E", p: Params = P):
    """
    Nose view: x-z plane at y = 0
    Horizontal axis = x
    Vertical axis   = z

    Force the field lines to wrap around the drawn circular cavity.
    """
    a = p.nose_radius
    r2 = x*x + z*z
    r = np.sqrt(max(r2, p.eps))

    # Inside cavity should never be traced, but guard anyway
    if r <= a:
        return np.array([np.nan, np.nan], dtype=float)

    if which == "E":
        # Uniform downward background field
        Ux, Uz = 0.0, 1.0
    elif which == "B":
        # Uniform rightward background field
        Ux, Uz = 1.0, 0.0
    else:
        raise ValueError("which must be 'E' or 'B'")

    theta = np.arctan2(z, x)

    # Cylinder flow in polar coordinates
    Ur = Ux * np.cos(theta) + Uz * np.sin(theta)
    Ut = -Ux * np.sin(theta) + Uz * np.cos(theta)

    ar2 = (a*a) / r2
    Vr = Ur * (1.0 - ar2)
    Vt = Ut * (1.0 + ar2)

    # Convert back to Cartesian
    Fx = Vr * np.cos(theta) - Vt * np.sin(theta)
    Fz = Vr * np.sin(theta) + Vt * np.cos(theta)

    return np.array([Fx, Fz], dtype=float)

# =============================================================================
# Field-line tracing
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
# Seeds
# =============================================================================
def make_edge_seeds_rect(xmin, xmax, ymin, ymax, n, mode="left"):
    s = np.linspace(-1.0, 1.0, n)

    # Nonuniform spacing: denser near center and flanks
    vals = np.linspace(ymin, ymax, n+100)

    if mode == "left":
        return np.column_stack([np.full_like(vals, xmin), vals])
    elif mode == "right":
        return np.column_stack([np.full_like(vals, xmax), vals])
    elif mode == "top":
        vals = 0.5 * (xmax - xmin) * np.sign(s) * np.abs(s)**1.6 + 0.5 * (xmax + xmin)
        return np.column_stack([vals, np.full_like(vals, ymax)])
    elif mode == "bottom":
        vals = 0.5 * (xmax - xmin) * np.sign(s) * np.abs(s)**1.6 + 0.5 * (xmax + xmin)
        return np.column_stack([vals, np.full_like(vals, ymin)])
    else:
        raise ValueError("mode must be 'left', 'right', 'top', or 'bottom'")


def make_side_halo_seeds(z0, n):
    t = np.linspace(-1.0, 1.0, n)
    z_vals = 1.9 * z0 * np.sign(t) * np.abs(t)**1.8

    rho = np.abs(z_vals) / z0
    y_hp = z0 * (2.0 - rho**2) / np.sqrt(np.maximum(4.0 - rho**2, 1e-12))

    y_seed = y_hp + 4.0
    return np.column_stack([y_seed, z_vals])


def make_nose_halo_seeds(radius, n):
    angles = np.linspace(0.18*np.pi, 0.82*np.pi, n // 2)
    angles = np.concatenate([-angles[::-1], angles])

    r_seed = 1.22 * radius
    x = r_seed * np.cos(angles[:n])
    z = r_seed * np.sin(angles[:n])
    return np.column_stack([x, z])

def make_nose_edge_seeds_E(xlim, zlim, n):
    """
    Uniform top/bottom seeds for the electric field.
    Avoid the centerline where lines stack vertically.
    """
    x = np.linspace(-0.75 * xlim, 0.75 * xlim, n)
    x = x[np.abs(x) > 18.0]
    top = np.column_stack([x, np.full_like(x,  zlim)])
    bot = np.column_stack([x, np.full_like(x, -zlim)])
    return np.vstack([top, bot])


def make_nose_edge_seeds_B(xlim, zlim, n):
    """
    Uniform left/right seeds for the magnetic field.
    Avoid the midplane where lines pile up horizontally.
    """
    z = np.linspace(-0.75 * zlim, 0.75 * zlim, n)
    z = z[np.abs(z) > 18.0]
    left  = np.column_stack([np.full_like(z, -xlim), z])
    right = np.column_stack([np.full_like(z,  xlim), z])
    return np.vstack([left, right])

def make_side_nose_following_seeds(z0):
    """
    Seeds placed just outside the heliopause so the side-view electric lines
    visibly wrap around the nose and flanks of the parabola.
    """
    z_vals = np.array([
        -175, -155, -135, -118, -102, -88, -75, -63, -52, -42, -33, -24,
          24,   33,   42,   52,   63,   75,   88,  102,  118,  135,  155, 175
    ], dtype=float)

    rho = np.abs(z_vals) / z0
    y_hp = z0 * (2.0 - rho**2) / np.sqrt(np.maximum(4.0 - rho**2, 1e-12))

    # Smaller offset near the nose so the field lines hug the boundary.
    offset = np.where(np.abs(z_vals) < 90.0, 2.5, 5.0)
    y_seed = y_hp + offset

    return np.column_stack([y_seed, z_vals])


# =============================================================================
# Domain checks
# =============================================================================
def inside_side_yz(y, z, p: Params):
    if (y < p.yz_y_min) or (y > p.yz_y_max) or (z < -p.yz_z_lim) or (z > p.yz_z_lim):
        return False
    return outside_heliopause_plot(0.0, y, z, p.z0)


def inside_nose_xz(x, z, p: Params):
    if (x < -p.xz_lim) or (x > p.xz_lim) or (z < -p.xz_lim) or (z > p.xz_lim):
        return False
    return (x*x + z*z) > p.nose_radius**2


# =============================================================================
# Plot helpers
# =============================================================================
def add_heliopause_side(ax, z0, y_down):
    """
    Parabolic side-view heliopause in the y-z plane.
    """
    z = np.linspace(-1.999 * z0, 1.999 * z0, 1400)
    rho = np.abs(z) / z0
    y = z0 * (2.0 - rho**2) / np.sqrt(np.maximum(4.0 - rho**2, 1e-12))

    z_fill = np.concatenate([z, z[::-1]])
    y_fill = np.concatenate([y, np.full_like(y, y_down)])

    ax.fill(y_fill, z_fill, color="white", zorder=5)
    ax.plot(y, z, color="black", linewidth=1.7, zorder=6)


def add_heliopause_nose(ax, radius):
    hp = Circle(
        (0.0, 0.0),
        radius,
        facecolor="white",
        edgecolor="black",
        linewidth=1.7,
        zorder=5
    )
    ax.add_patch(hp)


def format_axes(ax, xlabel, ylabel, xlim, ylim, equal=False):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.minorticks_on()
    ax.tick_params(direction="in", which="both", top=True, right=True)
    if equal:
        ax.set_aspect("equal", adjustable="box")


def add_arrow_to_line(ax, line, color, lw=1.0, arrow_size=9, frac=0.55):
    if len(line) < 8:
        return

    n = len(line)
    i = int(frac * (n - 2))
    i = np.clip(i, 2, n - 3)

    p0 = line[i]
    p1 = line[i + 1]

    if not (np.all(np.isfinite(p0)) and np.all(np.isfinite(p1))):
        return

    ax.annotate(
        "",
        xy=(p1[0], p1[1]),
        xytext=(p0[0], p0[1]),
        arrowprops=dict(
            arrowstyle="->",
            color=color,
            lw=lw,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=arrow_size,
        ),
        zorder=10,
    )


def plot_lines(ax, lines, color, lw=1.15, alpha=0.95, zorder=2, arrow_size=9):
    for line in lines:
        if len(line) >= 2:
            ax.plot(line[:, 0], line[:, 1], color=color, lw=lw, alpha=alpha, zorder=zorder)
            add_arrow_to_line(ax, line, color=color, lw=lw, arrow_size=arrow_size, frac=0.55)


# =============================================================================
# Figure builders
# =============================================================================
def build_side_figure(p: Params):
    fig, ax = plt.subplots(figsize=(9.1, 6.8))

    # Electric-field seeds:
    # combine far-left inflow seeds with seeds placed just outside the
    # heliopause so the red lines visibly bend around the nose.
    seeds_edge_E = make_edge_seeds_rect(
        p.yz_y_min, p.yz_y_max,
        -p.yz_z_lim, p.yz_z_lim,
        p.n_edge_seeds_side,
        mode="left",
    )
    seeds_nose_E = make_side_nose_following_seeds(p.z0)
    seeds_E = np.vstack([seeds_edge_E, seeds_nose_E])

    # Magnetic-field seeds:
    # keep the original pattern since B already looked right.
    seeds_edge_B = make_edge_seeds_rect(
        p.yz_y_min, p.yz_y_max,
        -p.yz_z_lim, p.yz_z_lim,
        p.n_edge_seeds_side,
        mode="left",
    )
    seeds_halo_B = make_side_halo_seeds(p.z0, p.n_halo_seeds_side)
    seeds_B = np.vstack([seeds_edge_B, seeds_halo_B])

    E_lines = []
    B_lines = []

    for seed in seeds_E:
        if inside_side_yz(seed[0], seed[1], p):
            lineE = trace_field_line(
                seed,
                lambda yy, zz: projected_field_side_yz(yy, zz, which="E", p=p),
                lambda yy, zz: inside_side_yz(yy, zz, p),
                p.ds,
                p.max_steps,
            )
            if len(lineE) > 10:
                E_lines.append(lineE)

    for seed in seeds_B:
        if inside_side_yz(seed[0], seed[1], p):
            lineB = trace_field_line(
                seed,
                lambda yy, zz: projected_field_side_yz(yy, zz, which="B", p=p),
                lambda yy, zz: inside_side_yz(yy, zz, p),
                p.ds,
                p.max_steps,
            )
            if len(lineB) > 10:
                B_lines.append(lineB)

    plot_lines(ax, E_lines, color="red",  lw=1.05, alpha=0.88, zorder=2, arrow_size=p.arrow_size)
    plot_lines(ax, B_lines, color="blue", lw=1.05, alpha=0.88, zorder=3, arrow_size=p.arrow_size)

    add_heliopause_side(ax, p.z0, y_down=p.yz_y_min - 20.0)

    format_axes(
        ax,
        xlabel=r"$x\ \mathrm{[AU]}$",
        ylabel=r"$y\ \mathrm{[AU]}$",
        xlim=(p.yz_y_min, p.yz_y_max),
        ylim=(-p.yz_z_lim, p.yz_z_lim),
        equal=False,
    )

    ax.set_title(r"Side view: $y$-$z$ plane at $x=0$")
    ax.text(
        0.02, 0.98,
        "red: electric field lines\nblue: magnetic field lines",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.95),
    )

    ax.text(
        0.83, 0.10,
        r"$\mathbf{B}_{\infty}\ \odot$",
        transform=ax.transAxes,
        fontsize=15,
        ha="left",
        va="center",
    )
    ax.text(
        0.83, 0.88,
        r"$\mathbf{B}_{\infty}\ \odot$",
        transform=ax.transAxes,
        fontsize=15,
        ha="left",
        va="center",
    )

    fig.tight_layout()
    return fig


def build_nose_figure(p: Params):
    fig, ax = plt.subplots(figsize=(7.2, 7.0))

    seeds_edge_E = make_nose_edge_seeds_E(p.xz_lim, p.xz_lim, 7)
    seeds_edge_B = make_nose_edge_seeds_B(p.xz_lim, p.xz_lim, 7)
    seeds_halo_E = make_nose_halo_seeds(p.nose_radius, 6)
    seeds_halo_B = make_nose_halo_seeds(p.nose_radius, 6)
    E_lines = []
    B_lines = []

    for seed in np.vstack([seeds_edge_E, seeds_halo_E]):
        if inside_nose_xz(seed[0], seed[1], p):
            line = trace_field_line(
                seed,
                lambda xx, zz: projected_field_nose_xz(xx, zz, which="E", p=p),
                lambda xx, zz: inside_nose_xz(xx, zz, p),
                p.ds,
                p.max_steps,
            )
            if len(line) > 10:
                E_lines.append(line)

    for seed in np.vstack([seeds_edge_B, seeds_halo_B]):
        if inside_nose_xz(seed[0], seed[1], p):
            line = trace_field_line(
                seed,
                lambda xx, zz: projected_field_nose_xz(xx, zz, which="B", p=p),
                lambda xx, zz: inside_nose_xz(xx, zz, p),
                p.ds,
                p.max_steps,
            )
            if len(line) > 10:
                B_lines.append(line)

    plot_lines(ax, E_lines, color="red",  lw=1.05, alpha=0.88, zorder=2, arrow_size=p.arrow_size)
    plot_lines(ax, B_lines, color="blue", lw=1.05, alpha=0.88, zorder=3, arrow_size=p.arrow_size)

    add_heliopause_nose(ax, p.nose_radius)

    format_axes(
        ax,
        xlabel=r"$x\ \mathrm{[AU]}$",
        ylabel=r"$y\ \mathrm{[AU]}$",
        xlim=(-p.xz_lim, p.xz_lim),
        ylim=(-p.xz_lim, p.xz_lim),
        equal=True,
    )

    # ax.set_title(r"Nose view: $x$-$y$ plane at $y=0$")
    # ax.text(
    #     0.02, 0.98,
    #     "red: electric field lines\nblue: magnetic field lines",
    #     transform=ax.transAxes,
    #     ha="left", va="top",
    #     fontsize=12,
    #     bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.95),
    # )

    fig.tight_layout()
    return fig


# =============================================================================
# Main
# =============================================================================
def main():
    fig_side = build_side_figure(P)
    fig_nose = build_nose_figure(P)

    fig_side.savefig("electric_field_side_yz.pdf")
    fig_side.savefig("electric_field_side_yz.png")
    fig_nose.savefig("magnetic_field_nose_xz.pdf")
    fig_nose.savefig("magnetic_field_nose_xz.png")

    plt.show()


if __name__ == "__main__":
    main()
