#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ============================================================
# Compact A&A journal style
# ============================================================

rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 1.0,
    "figure.dpi": 200
})

# ============================================================
# Constants
# ============================================================

G      = 6.67430e-11
M_sun  = 1.98847e30
GM     = G*M_sun
AU     = 1.496e11
c      = 299792458.0

v_inf  = 26e3
w_vec  = np.array([400e3, 0.0, 0.0])  # solar wind velocity for Lorentz term

B0     = 3.5e-9
v0     = 400e3

# ------------------------------------------------------------
# In your LaTeX model, ionopause radius at z=0 is sqrt(2)*z0.
# If you want ~40 AU cavity radius at z=0,
# choose z0 = 40 AU / sqrt(2) ≈ 28.3 AU.
# ------------------------------------------------------------
z0     = (40.0/np.sqrt(2.0)) * AU

# ============================================================
# Asymmetric beta curves (matching screenshot structure)
# ============================================================

def asym_log(m, peak, m0, sL, sR):
    logm = np.log10(m)
    logm0 = np.log10(m0)
    sigma = np.where(logm < logm0, sL, sR)
    return peak*np.exp(-(logm-logm0)**2/(2*sigma**2))

def set_stream_alpha(sp, alpha):
    """
    Matplotlib streamplot() returns a StreamplotSet with:
      - sp.lines  : LineCollection
      - sp.arrows : PatchCollection or list of patches (version-dependent)
    This sets transparency in a version-robust way.
    """
    if hasattr(sp, "lines") and sp.lines is not None:
        sp.lines.set_alpha(alpha)

    if hasattr(sp, "arrows") and sp.arrows is not None:
        arrows = sp.arrows
        # Newer mpl: PatchCollection-like
        if hasattr(arrows, "set_alpha"):
            arrows.set_alpha(alpha)
        else:
            # Older mpl: list of patches
            try:
                for a in arrows:
                    a.set_alpha(alpha)
            except TypeError:
                pass


# Carbon
def beta_carbon_0(m):  return asym_log(m,3.2,2e-17,0.6,1.0)
def beta_carbon_45(m): return asym_log(m,2.7,2e-17,0.7,1.1)
def beta_carbon_70(m): return asym_log(m,2.0,2e-17,0.8,1.2)

# Silicate
def beta_silicate_0(m):  return asym_log(m,0.85,5e-17,0.7,1.1)
def beta_silicate_45(m): return asym_log(m,0.55,5e-17,0.8,1.2)
def beta_silicate_70(m): return asym_log(m,0.30,5e-17,0.9,1.3)

# Astrosilicates
def beta_astrosil(m):     return asym_log(m,1.25,3e-17,0.65,1.05)
def beta_adapted(m):      return asym_log(m,1.6,3e-17,0.65,1.05)

# ============================================================
# FULL E and B fields (Horányi & Mendis-style; matches your LaTeX)
# ============================================================

def compute_fields(x, y, z=0.0):
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    # Promote z to same shape as x/y if scalar
    if z.shape == ():
        z = z + np.zeros_like(x)

    r2 = x**2 + y**2
    r  = np.sqrt(r2)
    R  = np.sqrt(r2 + z**2)

    # Avoid exact zeros (only for algebraic stability; cavity is handled by S-mask)
    r_safe = np.where(r == 0.0, 1e-30, r)
    R_safe = np.where(R == 0.0, 1e-30, R)

    # Streamlines term
    S2 = r2 + 2.0*(z0**2)*(z/R_safe - 1.0)

    # Cavity where S is not real (S2 <= 0)
    cavity = S2 <= 0.0
    S = np.sqrt(np.where(cavity, np.nan, S2))

    # Common prefactor for E
    prefE = (B0*v0)/(c*r_safe)

    # Helper factor appearing repeatedly
    with np.errstate(divide="ignore", invalid="ignore"):
        A = (-S/(r_safe**2) + 1.0/S - (z0**2 * z)/(R_safe**3 * S))

        Ex = prefE * (A * x * y)
        Ey = prefE * (S + A * (y**2))
        Ez = prefE * ((z0**2 * r2)/(R_safe**3 * S) * y)

    # Mask cavity hard to NaN (streamplot behaves well with NaNs)
    Ex = np.where(cavity, np.nan, Ex)
    Ey = np.where(cavity, np.nan, Ey)
    Ez = np.where(cavity, np.nan, Ez)

    # Magnetic field (Bz=0, Horányi 1985 form from your LaTeX)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = r_safe * (-1.0 + (z0**2 * z)/(R_safe**3))

        # bug-proof: if denom hits 0 exactly, keep a signed tiny number
        tiny = 1e-30
        denom = np.where(np.abs(denom) < tiny, np.where(denom >= 0, tiny, -tiny), denom)

        prefB = B0 / denom

        Bx = prefB * (-S + (S/(r_safe**2) - 1.0/S + (z0**2 * z)/(S * R_safe**3)) * (y**2))
        By = prefB * ((-S/(r_safe**2) + 1.0/S - (z0**2 * z)/(S * R_safe**3)) * (x * y))
        Bz = np.zeros_like(Bx)

    Bx = np.where(cavity, np.nan, Bx)
    By = np.where(cavity, np.nan, By)
    Bz = np.where(cavity, np.nan, Bz)

    E = np.stack((Ex, Ey, Ez), axis=0)
    B = np.stack((Bx, By, Bz), axis=0)
    return E, B

# ============================================================
# Full acceleration (gravity + rad pressure + Lorentz)
# ============================================================

def acceleration(state, m, Q, beta_func):
    x, y, vx, vy = state
    z = 0.0

    r_vec = np.array([x, y, z], dtype=float)
    v_vec = np.array([vx, vy, 0.0], dtype=float)

    r = np.linalg.norm(r_vec)
    if not np.isfinite(r) or r < 1e-30:
        r = 1e-30

    beta = float(beta_func(m))
    a_grav = -GM*(1.0-beta)*r_vec/(r**3)

    E, B = compute_fields(x, y, z)
    E = np.asarray(E, dtype=float).reshape(3,)
    B = np.asarray(B, dtype=float).reshape(3,)

    v_rel = v_vec - w_vec
    a_L = (Q/m) * (E + np.cross(v_rel, B))

    a_total = a_grav + a_L
    return np.array([vx, vy, a_total[0], a_total[1]], dtype=float)

# ============================================================
# RK4
# ============================================================

def rk4_step(state, dt, m, Q, beta_func):
    k1 = acceleration(state,               m, Q, beta_func)
    k2 = acceleration(state + 0.5*dt*k1,   m, Q, beta_func)
    k3 = acceleration(state + 0.5*dt*k2,   m, Q, beta_func)
    k4 = acceleration(state + dt*k3,       m, Q, beta_func)
    return state + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0

# ============================================================
# Integrate trajectory
# ============================================================

def integrate(m, Q, beta_func):
    state = np.array([-100*AU, 30*AU, v_inf, 0.0], dtype=float)
    dt   = 5e6
    tmax = 5e10
    traj = []

    t = 0.0
    while t < tmax:
        traj.append(state.copy())

        rr = np.sqrt(state[0]**2 + state[1]**2)
        if rr < 5*AU or rr > 150*AU:
            break

        state = rk4_step(state, dt, m, Q, beta_func)
        t += dt

        if not np.all(np.isfinite(state)):
            break

    return np.array(traj)

# ============================================================
# Field line plot
# ============================================================

def plot_fields():
    x = np.linspace(-100*AU, 100*AU, 320)
    y = np.linspace(-100*AU, 100*AU, 320)
    X, Y = np.meshgrid(x, y)

    E, B = compute_fields(X, Y, 0.0)
    Ex, Ey = E[0], E[1]
    Bx, By, Bz = B[0], B[1], B[2]

    fig, ax = plt.subplots(figsize=(4.6, 4.2))

    spE = ax.streamplot(X/AU, Y/AU, Ex, Ey, color="red",  density=1.0, linewidth=0.7, arrowsize=1.4)
    spB = ax.streamplot(X/AU, Y/AU, Bx, By, color="blue", density=1.0, linewidth=0.7, arrowsize=1.4)

    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
    ax.contour(X/AU, Y/AU, Bmag, levels=18, colors="blue", linewidths=0.35, alpha=0.55)

    cavity_r = np.sqrt(2.0)*z0/AU
    cavity = plt.Circle((0, 0), cavity_r, color="lightgrey", alpha=0.35, ec="grey", lw=0.8)
    ax.add_artist(cavity)

    sun = plt.Circle((0, 0), 2, color="gold", zorder=6)
    ax.add_artist(sun)

    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_aspect("equal")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.set_title("Electric (red) & Magnetic (blue) Fields (z=0 slice)")

    plt.tight_layout()
    plt.show()

# ============================================================
# Trajectory plot
# ============================================================

def plot_trajectories():
    masses = [5e-18, 2e-17, 1e-16]
    Q = 1e-16

    fig, ax = plt.subplots(figsize=(4.6, 4.2))

    for m in masses:
        traj = integrate(m, Q, beta_adapted)
        if traj.size == 0:
            continue
        ax.plot(traj[:, 0]/AU, traj[:, 1]/AU, lw=1.2, label=f"m={m:.1e}")

    sun = plt.Circle((0, 0), 2, color="gold", zorder=5)
    ax.add_artist(sun)

    cavity_r = np.sqrt(2.0)*z0/AU
    cavity = plt.Circle((0, 0), cavity_r, fill=False, ec="grey", lw=0.8, alpha=0.8)
    ax.add_artist(cavity)

    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_aspect("equal")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("Dust Trajectories (Full Equation; RK4)")

    plt.tight_layout()
    plt.show()

# ============================================================
# Combined plot (fields + trajectories superimposed)
#   - streamplot alpha handled via set_stream_alpha (no 'alpha' kwarg)
#   - trajectories drawn with black halo + order chosen so none are hidden
# ============================================================

def plot_fields_with_trajectories():
    # Grid for fields
    x = np.linspace(-100*AU, 100*AU, 320)
    y = np.linspace(-100*AU, 100*AU, 320)
    X, Y = np.meshgrid(x, y)

    E, B = compute_fields(X, Y, 0.0)
    Ex, Ey = E[0], E[1]
    Bx, By, Bz = B[0], B[1], B[2]

    fig, ax = plt.subplots(figsize=(4.9, 4.5))

    # Fields in the background
    spE = ax.streamplot(X/AU, Y/AU, Ex, Ey, color="red",  density=1.0, linewidth=0.60, arrowsize=1.15, zorder=1)
    spB = ax.streamplot(X/AU, Y/AU, Bx, By, color="blue", density=1.0, linewidth=0.60, arrowsize=1.15, zorder=2)

    # Make fields lighter so trajectories pop
    set_stream_alpha(spE, 0.45)
    set_stream_alpha(spB, 0.45)

    # |B| contours (optional visual aid)
    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
    ax.contour(X/AU, Y/AU, Bmag, levels=18, colors="blue", linewidths=0.30, alpha=0.35, zorder=0)

    # Cavity
    cavity_r = np.sqrt(2.0)*z0/AU
    cavity = plt.Circle((0, 0), cavity_r, color="lightgrey", alpha=0.35, ec="grey", lw=0.8, zorder=3)
    ax.add_artist(cavity)

    # Trajectories on top:
    # draw largest mass first, smallest last, so none get hidden
    masses = [1e-16, 2e-17, 5e-18]   # green, orange, blue (blue ends up on top)
    Q = 1e-16

    for m in masses:
        traj = integrate(m, Q, beta_adapted)
        if traj.size == 0:
            continue

        xtraj = traj[:, 0] / AU
        ytraj = traj[:, 1] / AU

        # black halo for visibility
        ax.plot(xtraj, ytraj, lw=3.2, color="k", zorder=7)
        # colored trajectory
        ax.plot(xtraj, ytraj, lw=1.9, zorder=8, label=f"m={m:.1e}")

    # Sun marker on very top
    sun = plt.Circle((0, 0), 2, color="gold", zorder=9)
    ax.add_artist(sun)

    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_aspect("equal")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    ax.set_title("Fields (streamlines) with Dust Trajectories (z=0 slice)")

    plt.tight_layout()
    plt.show()

# ============================================================
# Beta curve figure (all materials)
# ============================================================

def plot_beta_curves():
    m = np.logspace(-19, -12, 600)

    fig, ax = plt.subplots(figsize=(4.9, 3.8))

    ax.loglog(m, beta_carbon_0(m),  'k-',  label='Carbon p=0%')
    ax.loglog(m, beta_carbon_45(m), 'k--', label='Carbon p=45%')
    ax.loglog(m, beta_carbon_70(m), 'k:',  label='Carbon p=70%')

    ax.loglog(m, beta_adapted(m),   'r-',  label='Adapted Astron. sil.')
    ax.loglog(m, beta_astrosil(m),  'b:',  label='Astron. sil.')

    ax.loglog(m, beta_silicate_0(m),  color='limegreen',           label='Silicate p=0%')
    ax.loglog(m, beta_silicate_45(m), color='limegreen', linestyle='--', label='Silicate p=45%')
    ax.loglog(m, beta_silicate_70(m), color='limegreen', linestyle=':',  label='Silicate p=70%')

    ax.set_xlabel("Mass [kg]")
    ax.set_ylabel(r"$\beta$")
    ax.set_xlim(1e-19, 1e-12)
    ax.set_ylim(1e-2, 4)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.set_title(r"$\beta$ curves for multiple materials")

    plt.tight_layout()
    plt.show()

# ============================================================
# Run everything
# ============================================================

if __name__ == "__main__":
    plot_fields()
    plot_trajectories()
    plot_fields_with_trajectories()
    plot_beta_curves()
