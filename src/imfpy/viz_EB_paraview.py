#!/usr/bin/env python3
"""
Publication-quality comet E/B plots on z = 0.

• Exact E, B formulas from the manuscript (Bz = 0 case)
• Physical mask where S^2 < 0  + small center mask
• Direction-only arrows (constant length), magnitude shown by background colormap
• Tidy typography, consistent sizing, vector outputs (PDF/SVG) + high-DPI PNG

Outputs:
  - figs/E_field.png,  figs/E_field.pdf
  - figs/B_field.png,  figs/B_field.pdf
  - figs/EB_side_by_side.png, figs/EB_side_by_side.pdf
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ------------------------ Figure style (journal-ready) ------------------------
mpl.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.pad_inches": 0.02,
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "grid.alpha": 0.25,
})
from matplotlib.ticker import ScalarFormatter

# ---------------------------- Constants -----------------------------
B0 = 1.0     # far-field B_x magnitude (units arbitrary/consistent)
v0 = 1.0     # far-field flow speed along -z
z0 = 0.5     # subsolar ionopause position
c  = 1.0     # speed of light

# ---------------------------- Grid setup ----------------------------
L = 2.0          # domain half-size in x,y
N = 301          # dense background for smooth shading
x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
X, Y = np.meshgrid(x, y, indexing='xy')
Z = np.zeros_like(X)  # z = 0 slice

EPS = 1e-12

# ---------------------- Field helper quantities ---------------------
r = np.hypot(X, Y)
r = np.where(r == 0.0, EPS, r)
R = np.hypot(r, Z)
R = np.where(R == 0.0, EPS, R)

# Streamline factor S^2 = r^2 + 2 z0^2 ( z/R - 1 ) ; at z=0 -> r^2 - 2 z0^2
S_sq_raw = r**2 + 2.0 * z0**2 * (Z / R - 1.0)
valid_phys = S_sq_raw >= 0.0                       # physical region where S is real
S_sq = np.where(valid_phys, S_sq_raw, 0.0)
S = np.sqrt(S_sq)
S = np.where((S == 0.0) & valid_phys, EPS, S)      # avoid divide-by-zero inside valid region

# ------------------------------- E field -----------------------------
coeff = (B0 * v0) / (c * r)
common = (-S / (r**2)) + (1.0 / S) - (z0**2 * Z) / (S * R**3)
Ex = coeff * common * (X * Y)
Ey = coeff * (S + common * (Y**2))
Ez = coeff * ((z0**2 * r**2) / (S * R**3)) * Y

# ------------------------------- B field -----------------------------
den = -1.0 + (z0**2) * Z / (R**3)       # = v_z / v0
den = np.where(np.abs(den) < EPS, np.sign(den) * EPS, den)
Bx = B0 * ( (( S/(r**2) - 1.0/S + (z0**2)*Z/(S*R**3) ) * (Y**2) - S ) / (r * den) )
By = B0 * ( ( -S/(r**2) + 1.0/S - (z0**2)*Z/(S*R**3) ) * (X*Y) / (r * den) )
Bz = np.zeros_like(Bx)

# --------------- Masks: physical + tiny center to tame spikes -------
center_mask = (np.hypot(X, Y) < (0.06 * L))
mask = (~valid_phys) | center_mask

# ---------------------------- Helpers -------------------------------
def clip_percentile(arr, p=99.0, mask=None):
    """Clip by percentile (computed over unmasked/valid entries)."""
    if mask is None:
        thr = np.percentile(arr, p)
    else:
        valid_vals = arr[~mask]
        thr = np.percentile(valid_vals, p) if valid_vals.size else np.percentile(arr, p)
    return np.minimum(arr, thr), thr

def make_axes_equal(ax, L):
    ax.set_aspect('equal', 'box')
    ax.set_xlim([-L, L]); ax.set_ylim([-L, L])
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.grid(True, linestyle=':', linewidth=0.6)

def draw_ionopause(ax):
    """r = sqrt(2) z0 circle at z=0 where S^2 = 0 (separatrix)."""
    import matplotlib.patches as mpatches
    rho = np.sqrt(2.0) * z0
    if rho > 0:
        circ = mpatches.Circle((0,0), radius=rho, fill=False, lw=1.1, ls='--', ec='k', alpha=0.9)
        ax.add_patch(circ)

def quiver_dirplot(X, Y, U, V, title, fname_base, cmap='magma',
                   step=None, arrow_len=0.10, pct=99.0):
    """
    Direction quiver + magnitude background.
    Saves PNG (300 dpi) and PDF.
    """
    os.makedirs("figs", exist_ok=True)

    mag = np.hypot(U, V)
    mag_clip, thr = clip_percentile(mag, p=pct, mask=mask)

    # prepare masked arrays for clean blank region
    M = np.ma.masked_where(mask, mag_clip)

    # downsample factor for arrows (aim for ~40x40 arrows)
    if step is None:
        step = max(1, int(np.floor(min(X.shape) / 30)))
    s = (slice(None, None, step), slice(None, None, step))

    Udir = np.where(mask, np.nan, U / (mag + EPS))
    Vdir = np.where(mask, np.nan, V / (mag + EPS))

    Xq, Yq = X[s], Y[s]
    Uq, Vq = Udir[s], Vdir[s]

    fig, ax = plt.subplots(figsize=(3.6, 3.3), constrained_layout=True)
    # magnitude background
    im = ax.pcolormesh(X, Y, M, shading='auto', cmap=cmap)
    make_axes_equal(ax, L)
    draw_ionopause(ax)

    # quiver: fixed-length, high-contrast
    q = ax.quiver(Xq, Yq, Uq*arrow_len, Vq*arrow_len, color='k',
                  pivot='mid', angles='xy', scale_units='xy', scale=1,
                  headwidth=3.2, headlength=4.8, linewidths=0.5, alpha=0.95)

    # colorbar with scientific formatter
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
    cb.set_label(f'|{title.split()[0]}|')
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2, 3))  # switch to sci-notation outside this range
    cb.formatter = fmt
    cb.update_ticks()

    ax.set_title(title)

    # save
    for ext in ("png", "pdf"):
        fig.savefig(f"figs/{fname_base}.{ext}", bbox_inches='tight', transparent=False)
    plt.close(fig)

# ---------------------------- Make plots ----------------------------
quiver_dirplot(X, Y, Ex, Ey, 'E field (z=0)', 'E_field',
               cmap='magma', arrow_len=0.11, pct=99.0)

quiver_dirplot(X, Y, Bx, By, 'B field (z=0)', 'B_field',
               cmap='magma', arrow_len=0.11, pct=99.0)

# --------------------- Side-by-side (no quivers) --------------------
os.makedirs("figs", exist_ok=True)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.3), constrained_layout=True)

for ax, U, V, title in [(ax1, Ex, Ey, 'E field (z=0)'),
                        (ax2, Bx, By, 'B field (z=0)')]:
    mag = np.hypot(U, V)
    mag_clip, _ = clip_percentile(mag, p=99.0, mask=mask)
    M = np.ma.masked_where(mask, mag_clip)
    im = ax.pcolormesh(X, Y, M, shading='auto', cmap='magma')
    make_axes_equal(ax, L)
    draw_ionopause(ax)
    ax.set_title(title)

# shared colorbar
cb = fig.colorbar(im, ax=[ax1, ax2], pad=0.02, fraction=0.05)
cb.set_label('|Field|')
fmt = ScalarFormatter(useMathText=True); fmt.set_powerlimits((-2, 3))
cb.formatter = fmt; cb.update_ticks()

for ext in ("png", "pdf"):
    fig.savefig(f"figs/EB_side_by_side.{ext}", bbox_inches='tight', transparent=False)
plt.close(fig)

print("Saved: figs/E_field.(png,pdf), figs/B_field.(png,pdf), figs/EB_side_by_side.(png,pdf)")
