import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import matplotlib as mpl

# ============================================================
# Style
# ============================================================
mpl.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.linewidth": 1.0,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

# ============================================================
# Load data
# ============================================================
data = np.loadtxt("coord.dat")

t  = data[:, 0]   # years
x  = data[:, 1]   # AU
y  = data[:, 2]   # AU
vx = data[:, 4]
vy = data[:, 5]
vz = data[:, 6]

speed = np.sqrt(vx**2 + vy**2 + vz**2)

# ============================================================
# Split trajectories
# ============================================================
dt = np.diff(t)
breaks = np.where(dt < 0)[0] + 1
tracks = np.split(np.arange(len(t)), breaks)

# decimate particle population
tracks = tracks[::30]

# ============================================================
# Select an "interesting" trajectory
# ============================================================
def trajectory_score(idx):
    xm = x[idx]
    ym = y[idx]
    vm = speed[idx]

    if len(idx) < 50:
        return -np.inf

    r = np.sqrt(xm**2 + ym**2)
    r_min = np.min(r)
    dv    = np.max(vm) - np.min(vm)
    dphi  = np.ptp(np.unwrap(np.arctan2(ym, xm)))

    return 2.0 / (r_min + 0.5) + 1.5 * dv + dphi


scores = np.array([trajectory_score(idx) for idx in tracks])
best_idx = tracks[np.argmax(scores)]

print(f"Selected trajectory score = {scores.max():.3f}")

# ============================================================
# Extract ONLY that trajectory
# ============================================================
xm = x[best_idx]
ym = y[best_idx]
tm = t[best_idx] - t[best_idx][0]   # relative time

# ============================================================
# HARD decimation (THIS MATTERS)
# ============================================================
STEP = 10
xm = xm[::STEP]
ym = ym[::STEP]
tm = tm[::STEP]

# ============================================================
# Perihelion
# ============================================================
r = np.sqrt(xm**2 + ym**2)
i_peri = np.argmin(r)

x_peri = xm[i_peri]
y_peri = ym[i_peri]
t_peri = tm[i_peri]

# ============================================================
# Build segments (ONLY ONE TRAJECTORY)
# ============================================================
pts = np.column_stack((xm, ym))
segs = np.stack([pts[:-1], pts[1:]], axis=1)

# time coloring (clipped)
TMAX = 500.0
t_color = np.clip(tm[:-1], 0, TMAX)
norm = Normalize(vmin=0, vmax=TMAX)

# ============================================================
# STATIC PLOT
# ============================================================
fig, ax = plt.subplots(figsize=(4.4, 5.4))

lc = LineCollection(
    segs,
    cmap="plasma",
    norm=norm,
    linewidth=1.6,
    alpha=0.95
)
lc.set_array(t_color)
ax.add_collection(lc)

# Sun
ax.plot(0, 0, "o", ms=6, mfc="#d4af37", mec="black", zorder=10)
ax.text(0.4, 0.3, "Sun", fontsize=9)

# Perihelion
ax.plot(x_peri, y_peri, "x", color="black", ms=7, zorder=11)
ax.annotate(
    f"Perihelion\n(t = {t_peri:.1f} yr)",
    xy=(x_peri, y_peri),
    xytext=(x_peri + 1.5, y_peri - 1.5),
    arrowprops=dict(arrowstyle="->", lw=0.8),
    fontsize=8
)

# Zoom
RZOOM = 6.0
ax.set_xlim(x_peri - RZOOM, x_peri + RZOOM)
ax.set_ylim(y_peri - RZOOM, y_peri + RZOOM)

ax.set_aspect("equal")
ax.set_xlabel("x [AU]")
ax.set_ylabel("y [AU]")

cbar = fig.colorbar(lc, ax=ax, pad=0.02, fraction=0.05)
cbar.set_label("Time since entry [yr] (clipped)")

ax.set_title("Interstellar Dust Trajectory (Selected Grain)")

plt.tight_layout()
plt.savefig("isd_single_trajectory_clean.png")
plt.savefig("isd_single_trajectory_clean.pdf")
plt.show()
