import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# ============================================================
# Physical constants (Juhász & Horányi 2013)
# ============================================================
AU = 1.496e11
Omega = 2.86533e-6       # rad/s (solar rotation)
w = 400e3               # m/s (solar wind)
rs = 2.5 * 6.96e8       # source surface radius
alpha = np.deg2rad(25)  # HCS tilt angle
phi0 = 0.0

# ============================================================
# Grid
# ============================================================
r = np.linspace(0.3 * AU, 1.2 * AU, 140)
phi = np.linspace(0, 2*np.pi, 360)
R, PHI = np.meshgrid(r, phi)

# ============================================================
# Current sheet equation (Eq. 3)
# ============================================================
def current_sheet_theta(R, PHI, t):
    arg = PHI - phi0 - Omega*t + Omega*(R - rs)/w
    theta = np.pi/2 - np.arctan(np.tan(alpha) * np.sin(arg))
    return theta

# ============================================================
# Figure
# ============================================================
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
fig.patch.set_facecolor("black")
ax.set_facecolor("black")

# ============================================================
# Precompute color normalization (fixed, not flickering)
# ============================================================
lat_max = np.rad2deg(alpha)
norm = mpl.colors.Normalize(vmin=-lat_max, vmax=lat_max)
cmap = mpl.cm.coolwarm

mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([])

cbar = fig.colorbar(
    mappable,
    ax=ax,
    shrink=0.6,
    pad=0.05
)
cbar.set_label(
    r"HCS latitude offset  [$^\circ$]",
    color="white",
    fontsize=11
)
cbar.ax.yaxis.set_tick_params(color="white")
plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="white")

# ============================================================
# Time annotation (figure-level, survives ax.cla())
# ============================================================
time_text = fig.text(
    0.05, 0.93, "",
    color="white",
    fontsize=12
)


# ============================================================
# Animation update
# ============================================================
def update(frame):
    ax.cla()
    ax.set_facecolor("black")

    # Time (seconds → days)
    t = frame * 3.0e5
    t_days = t / 86400.0

    theta = current_sheet_theta(R, PHI, t)
    lat_offset = np.rad2deg(theta - np.pi/2)

    X = (R/AU) * np.sin(theta) * np.cos(PHI)
    Y = (R/AU) * np.sin(theta) * np.sin(PHI)
    Z = (R/AU) * np.cos(theta)

    ax.plot_surface(
        X, Y, Z,
        facecolors=cmap(norm(lat_offset)),
        linewidth=0,
        antialiased=True,
        alpha=0.9
    )

    # Sun
    ax.scatter(0, 0, 0, color="gold", s=140)

    # Axes limits & view
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_zlim(-1.3, 1.3)
    ax.set_axis_off()
    ax.view_init(elev=25, azim=frame * 0.6)

    # Time label
    time_text.set_text(f"Time = {t_days:5.1f} days")

    return []

ani = FuncAnimation(
    fig,
    update,
    frames=180,
    interval=60,
    blit=False
)

plt.show()

ani.save("current_sheet.mp4", dpi=200)
