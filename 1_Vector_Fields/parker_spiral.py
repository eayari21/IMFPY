import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ============================================================
# Physical constants
# ============================================================
AU = 1.496e11
B0 = 3.5e-9          # Tesla at 1 AU
Omega = 2.86533e-6  # rad/s (solar rotation)
w = 400e3           # m/s (solar wind)
r0 = AU
rs = 2.5 * 6.96e8

# ============================================================
# Animation tuning (THIS is the secret sauce)
# ============================================================
dt = 3.0e4          # seconds per frame  (~8.3 hours)
FPS = 30            # perceived smoothness
FRAMES = 300

# ============================================================
# Grid (Cartesian, moderate resolution for speed)
# ============================================================
n = 240
x = np.linspace(-1.3, 1.3, n)
y = np.linspace(-1.3, 1.3, n)
X, Y = np.meshgrid(x, y)

R = np.sqrt(X**2 + Y**2) * AU
PHI = np.arctan2(Y, X)

# Mask inside the Sun
R[R < 0.1 * AU] = np.nan

# ============================================================
# Parker spiral magnetic field
# ============================================================
def parker_field(R, PHI, t):
    Br = B0 * (r0 / R)**2
    Bphi = -Br * Omega * (R - rs) / w

    PHI_t = PHI - Omega * t

    Bx = Br * np.cos(PHI_t) - Bphi * np.sin(PHI_t)
    By = Br * np.sin(PHI_t) + Bphi * np.cos(PHI_t)

    return Bx, By

# ============================================================
# Figure
# ============================================================
fig, ax = plt.subplots(figsize=(7, 7))
fig.patch.set_facecolor("black")

# ============================================================
# Animation update
# ============================================================
def update(frame):
    ax.cla()
    ax.set_facecolor("black")

    t = frame * dt
    t_days = t / 86400.0

    Bx, By = parker_field(R, PHI, t)
    speed = np.sqrt(Bx**2 + By**2)

    ax.streamplot(
        X, Y,
        Bx, By,
        color=np.log10(speed),
        cmap="plasma",
        density=1.25,       # lighter = faster
        linewidth=1.1,
        arrowsize=1.6       # BIGGER arrowheads
    )

    # Sun
    ax.scatter(0, 0, s=240, color="gold", zorder=10)

    # Formatting
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.axis("off")

    # Time annotation
    ax.text(
        0.03, 0.95,
        f"Time = {t_days:5.1f} days",
        transform=ax.transAxes,
        color="white",
        fontsize=12
    )

    return []

# ============================================================
# Animate
# ============================================================
ani = FuncAnimation(
    fig,
    update,
    frames=FRAMES,
    interval=1000 / FPS,
    blit=False
)

plt.show()
ani.save("parker_spiral.mp4", dpi=200)
