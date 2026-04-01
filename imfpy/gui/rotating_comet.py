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
w_vec  = np.array([400e3, 0.0, 0.0])

B0     = 3.5e-9
v0     = 400e3

E0x = 0.0
E0y = 1.0   # or -1.0 if you want downward
B0x = 1.0
B0y = 0.0

# ============================================================
# Magnetic field rotation angle (radians)
# ============================================================

B_ROTATION = 0.0

# ============================================================
# Parameter sweep configuration
# ============================================================

B_ROTATIONS = [0.0, np.pi/2, np.pi, 3*np.pi/2]
B_STRENGTHS = [2.0e-9, 3.5e-9, 5.0e-9]

# ============================================================
# Ionopause geometry
# ============================================================

z0 = (40.0/np.sqrt(2.0)) * AU

# ============================================================
# Utility functions
# ============================================================

def rotate_about_z(vec, theta):
    cth = np.cos(theta)
    sth = np.sin(theta)
    R = np.array([[ cth, -sth, 0.0],
                  [ sth,  cth, 0.0],
                  [ 0.0,  0.0, 1.0]])
    return np.tensordot(R, vec, axes=(1, 0))


def find_extreme_mass(beta_func):
    mgrid = np.logspace(-19, -12, 5000)
    beta_vals = beta_func(mgrid)
    idx = np.argmax(beta_vals)
    return mgrid[idx], beta_vals[idx]


def asym_log(m, peak, m0, sL, sR):
    logm = np.log10(m)
    logm0 = np.log10(m0)
    sigma = np.where(logm < logm0, sL, sR)
    return peak*np.exp(-(logm-logm0)**2/(2*sigma**2))


def charge_from_mass(m):
    rho = 3000.0
    a = (3*m/(4*np.pi*rho))**(1/3)
    phi = 5.0
    eps0 = 8.854e-12
    return 4*np.pi*eps0*a*phi


def set_stream_alpha(sp, alpha):
    if hasattr(sp, "lines") and sp.lines is not None:
        sp.lines.set_alpha(alpha)
    if hasattr(sp, "arrows") and sp.arrows is not None:
        arrows = sp.arrows
        if hasattr(arrows, "set_alpha"):
            arrows.set_alpha(alpha)
        else:
            try:
                for a in arrows:
                    a.set_alpha(alpha)
            except TypeError:
                pass

# ============================================================
# Beta curves
# ============================================================

def beta_carbon_0(m):  return asym_log(m,3.2,2e-17,0.6,1.0)
def beta_carbon_45(m): return asym_log(m,2.7,2e-17,0.7,1.1)
def beta_carbon_70(m): return asym_log(m,2.0,2e-17,0.8,1.2)

def beta_silicate_0(m):  return asym_log(m,0.85,5e-17,0.7,1.1)
def beta_silicate_45(m): return asym_log(m,0.55,5e-17,0.8,1.2)
def beta_silicate_70(m): return asym_log(m,0.30,5e-17,0.9,1.3)

def beta_astrosil(m): return asym_log(m,1.25,3e-17,0.65,1.05)
def beta_adapted(m):  return asym_log(m,1.6,3e-17,0.65,1.05)

# ============================================================
# Field model
# ============================================================

def compute_fields(x, y, z=0.0, B_rotation=0.0):
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    if z.shape == ():
        z = z + np.zeros_like(x)

    x = x - 40*AU
    r2 = x**2 + y**2
    r  = np.sqrt(r2)
    R  = np.sqrt(r2 + z**2)

    r_safe = np.where(r == 0.0, 1e-30, r)
    R_safe = np.where(R == 0.0, 1e-30, R)

    r = np.sqrt(r2)
    r_safe = np.where(r == 0.0, 1e-30, r)

    # Use x as heliospheric axis so the z=0 slice shows the paraboloid
    # Shift so nose is near origin (tune if needed)
    # Move the nose left of the Sun
    # Nose location
    x_nose = -60 * AU

    # Shift coordinates so the cavity nose sits at x_nose
    x_ax = x - x_nose

    # Cylindrical radius about the x-axis
    r2 = y**2 + z**2
    R  = np.sqrt(x_ax**2 + r2)
    R_safe = np.where(R == 0.0, 1e-30, R)

    # Rounded comet-like cavity opening to +x
    z0_eff = 0.75 * z0
    S2 = r2 + 2.0 * z0_eff**2 * (-x_ax / R_safe - 1.0)

    cavity = S2 <= 0.0
    S = np.sqrt(np.where(cavity, np.nan, S2))

    # Distance-like coordinate outside the cavity
    S_pos = np.sqrt(np.maximum(S2, 0.0))

    # Relaxation length: larger = fields stay distorted farther away
    L_relax = 35.0 * AU

    # f = 1 near boundary, f -> 0 far away
    f = np.exp(-(S_pos / L_relax)**2)
    f = np.where(cavity, np.nan, f)

    prefE = (B0*v0)/(c*r_safe)

    with np.errstate(divide="ignore", invalid="ignore"):
        A = (-S/(r_safe**2) + 1.0/S - (z0**2 * z)/(R_safe**3 * S))
        Ex = prefE * (A * x * y)
        Ey = prefE * (S + A * (y**2))
        Ez = prefE * ((z0**2 * r2)/(R_safe**3 * S) * y)

    Ex = np.where(cavity, np.nan, Ex)
    Ey = np.where(cavity, np.nan, Ey)
    Ez = np.where(cavity, np.nan, Ez)

    with np.errstate(divide="ignore", invalid="ignore"):
        denom = r_safe * (-1.0 + (z0**2 * z)/(R_safe**3))
        tiny = 1e-30
        denom = np.where(np.abs(denom) < tiny,
                         np.where(denom >= 0, tiny, -tiny),
                         denom)
        prefB = B0 / denom
        Bx = prefB * (-S + (S/(r_safe**2) - 1.0/S +
             (z0**2 * z)/(S * R_safe**3)) * (y**2))
        By = prefB * ((-S/(r_safe**2) + 1.0/S -
             (z0**2 * z)/(S * R_safe**3)) * (x * y))
        Bz = np.zeros_like(Bx)

    Bx = np.where(cavity, np.nan, Bx)
    By = np.where(cavity, np.nan, By)
    Bz = np.where(cavity, np.nan, Bz)

    # Normalize distorted fields before blending
    Emag = np.sqrt(Ex**2 + Ey**2)
    Bmag = np.sqrt(Bx**2 + By**2)

    Emag_safe = np.where((~np.isfinite(Emag)) | (Emag == 0.0), 1.0, Emag)
    Bmag_safe = np.where((~np.isfinite(Bmag)) | (Bmag == 0.0), 1.0, Bmag)

    Ex_d = Ex / Emag_safe
    Ey_d = Ey / Emag_safe
    Bx_d = Bx / Bmag_safe
    By_d = By / Bmag_safe

    # Blend distorted fields near cavity with uniform far-field components
    Ex = f * Ex_d + (1.0 - f) * E0x
    Ey = f * Ey_d + (1.0 - f) * E0y

    Bx = f * Bx_d + (1.0 - f) * B0x
    By = f * By_d + (1.0 - f) * B0y

    E = np.stack((Ex, Ey, Ez), axis=0)
    B = np.stack((Bx, By, Bz), axis=0)


    if B_rotation != 0.0:
        B = rotate_about_z(B, B_rotation)

    return E, B

# ============================================================
# Dynamics
# ============================================================

def acceleration(state, m, Q, beta_func):
    x, y, z, vx, vy, vz = state

    r_vec = np.array([x, y, z], dtype=float)
    v_vec = np.array([vx, vy, vz], dtype=float)

    r = np.linalg.norm(r_vec)
    if not np.isfinite(r) or r < 1e-30:
        r = 1e-30

    beta = float(beta_func(m))
    a_grav = -GM*(1.0-beta)*r_vec/(r**3)

    E, B = compute_fields(x, y, z, B_rotation=B_ROTATION)
    E = np.asarray(E, dtype=float).reshape(3,)
    B = np.asarray(B, dtype=float).reshape(3,)

    Q_eff = charge_from_mass(m)
    a_L = (Q_eff/m) * (E + np.cross(v_vec, B))

    a_total = a_grav + a_L

    return np.array([vx, vy, vz,
                     a_total[0], a_total[1], a_total[2]],
                     dtype=float)




def rk4_step(state, dt, m, Q, beta_func):
    k1 = acceleration(state,               m, Q, beta_func)
    k2 = acceleration(state + 0.5*dt*k1,   m, Q, beta_func)
    k3 = acceleration(state + 0.5*dt*k2,   m, Q, beta_func)
    k4 = acceleration(state + dt*k3,       m, Q, beta_func)
    return state + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0


def integrate(m, Q, beta_func):
    state = np.array([-100*AU, 15*AU, 0.0,
                       v_inf, 0.0, 0.0], dtype=float)

    dt   = 5e5
    tmax = 5e10
    traj = []

    t = 0.0
    while t < tmax:
        traj.append(state.copy())
        rr = np.linalg.norm(state[:3])
        if rr < 5*AU or rr > 150*AU:
            break
        state = rk4_step(state, dt, m, Q, beta_func)
        t += dt
        if not np.all(np.isfinite(state)):
            break

    return np.array(traj)

# ============================================================
# Plotting functions (UNCHANGED behavior)
# ============================================================

def plot_fields():
    x = np.linspace(-100*AU, 100*AU, 320)
    y = np.linspace(-100*AU, 100*AU, 320)
    X, Y = np.meshgrid(x, y)

    E, B = compute_fields(X, Y, 0.0, B_rotation=B_ROTATION)
    Ex, Ey = E[0], E[1]
    Bx, By, Bz = B[0], B[1], B[2]

    fig, ax = plt.subplots(figsize=(4.6, 4.2))

    spE = ax.streamplot(X/AU, Y/AU, Ex, Ey,
                        color="red", density=1.0,
                        linewidth=0.7, arrowsize=1.4)
    
        # --- cavity boundary (same geometry as in compute_fields) ---
    z0_eff = 0.75 * z0

    # Must match the coordinate shifts used in compute_fields
    Xs = X - 40*AU
    x_nose = -60 * AU
    Xax = Xs - x_nose

    # z = 0 slice, so cylindrical radius about x-axis is just |y|
    r2_cav = Y**2
    R_cav = np.sqrt(Xax**2 + r2_cav)
    R_cav = np.where(R_cav == 0.0, 1e-30, R_cav)

    S2_cav = r2_cav + 2.0 * z0_eff**2 * (-Xax / R_cav - 1.0)

    ax.contour(X/AU, Y/AU, S2_cav,
               levels=[0.0], colors="black",
               linewidths=1.6, zorder=5)
    # spB = ax.streamplot(X/AU, Y/AU, Bx, By,
                        # color="blue", density=1.0,
                        # linewidth=0.7, arrowsize=1.4)

    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
    # ax.contour(X/AU, Y/AU, Bmag,
    #            levels=18, colors="blue",
    #            linewidths=0.35, alpha=0.55)

    

    sun = plt.Circle((0,0), 2,
                     color="gold", zorder=6)
    ax.add_artist(sun)

    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    ax.set_aspect("equal")
    ax.set_xlabel("z [AU]")
    ax.set_ylabel("y [AU]")
    # ax.set_title("Electric (red) & Magnetic (blue) Fields (z=0 slice)")

    plt.tight_layout()
    plt.show()


def plot_trajectories():
    m_extreme, beta_extreme = find_extreme_mass(beta_adapted)
    print(f"Extreme beta = {beta_extreme:.3f} at m = {m_extreme:.3e} kg")

    masses = [m_extreme]
    Q = 1e-16

    fig, ax = plt.subplots(figsize=(4.6, 4.2))

    for m in masses:
        traj = integrate(m, Q, beta_adapted)
        if traj.size == 0:
            continue
        ax.plot(traj[:,0]/AU, traj[:,1]/AU,
                lw=1.2, label=f"m={m:.1e}")

    sun = plt.Circle((0,0), 2,
                     color="gold", zorder=5)
    ax.add_artist(sun)

    

    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    ax.set_aspect("equal")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("Dust Trajectories (Full Equation; RK4)")

    plt.tight_layout()
    plt.show()


def plot_fields_with_trajectories():
    x = np.linspace(-100*AU, 100*AU, 320)
    y = np.linspace(-100*AU, 100*AU, 320)
    X, Y = np.meshgrid(x, y)

    E, B = compute_fields(X, Y, 0.0,
                          B_rotation=B_ROTATION)
    Ex, Ey = E[0], E[1]
    Bx, By, Bz = B[0], B[1], B[2]

    fig, ax = plt.subplots(figsize=(4.9, 4.5))

    spE = ax.streamplot(X/AU, Y/AU, Ex, Ey,
                        color="red", density=1.0,
                        linewidth=0.60, arrowsize=1.15, zorder=1)
    spB = ax.streamplot(X/AU, Y/AU, Bx, By,
                        color="blue", density=1.0,
                        linewidth=0.60, arrowsize=1.15, zorder=2)

    set_stream_alpha(spE, 0.45)
    set_stream_alpha(spB, 0.45)

    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
    ax.contour(X/AU, Y/AU, Bmag,
               levels=18, colors="blue",
               linewidths=0.30, alpha=0.35, zorder=0)

    

    masses = [1e-16, 2e-17, 5e-18]
    Q = 1e-16

    for m in masses:
        traj = integrate(m, Q, beta_adapted)
        if traj.size == 0:
            continue

        xtraj = traj[:,0]/AU
        ytraj = traj[:,1]/AU

        ax.plot(xtraj, ytraj,
                lw=3.2, color="k", zorder=7)
        ax.plot(xtraj, ytraj,
                lw=1.9, zorder=8,
                label=f"m={m:.1e}")

    sun = plt.Circle((0,0), 2,
                     color="gold", zorder=9)
    ax.add_artist(sun)

    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    ax.set_aspect("equal")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.legend(fontsize=7,
              frameon=False, loc="lower left")
    ax.set_title("Fields (streamlines) with Dust Trajectories (z=0 slice)")

    plt.tight_layout()
    plt.show()


def plot_beta_curves():
    m = np.logspace(-19, -12, 600)

    fig, ax = plt.subplots(figsize=(4.9, 3.8))

    ax.loglog(m, beta_carbon_0(m),'k-',label='Carbon p=0%')
    ax.loglog(m, beta_carbon_45(m),'k--',label='Carbon p=45%')
    ax.loglog(m, beta_carbon_70(m),'k:',label='Carbon p=70%')

    ax.loglog(m, beta_adapted(m),'r-',label='Adapted Astron. sil.')
    ax.loglog(m, beta_astrosil(m),'b:',label='Astron. sil.')

    ax.loglog(m, beta_silicate_0(m),
              color='limegreen',label='Silicate p=0%')
    ax.loglog(m, beta_silicate_45(m),
              color='limegreen',linestyle='--',
              label='Silicate p=45%')
    ax.loglog(m, beta_silicate_70(m),
              color='limegreen',linestyle=':',
              label='Silicate p=70%')

    ax.set_xlabel("Mass [kg]")
    ax.set_ylabel(r"$\beta$")
    ax.set_xlim(1e-19,1e-12)
    ax.set_ylim(1e-2,4)
    ax.legend(fontsize=7,frameon=False,
              loc="upper right")
    ax.set_title(r"$\beta$ curves for multiple materials")

    plt.tight_layout()
    plt.show()

# ============================================================
# Parameter sweep
# ============================================================

def run_parameter_sweep():
    global B_ROTATION, B0

    original_rotation = B_ROTATION
    original_B0 = B0

    for B_rot in B_ROTATIONS:
        for B_strength in B_STRENGTHS:

            print("\n====================================")
            print(f"B rotation = {B_rot:.2f} rad")
            print(f"B strength = {B_strength:.2e} T")
            print("====================================")

            B_ROTATION = B_rot
            B0 = B_strength

            plot_fields_with_trajectories()

    B_ROTATION = original_rotation
    B0 = original_B0

# ============================================================
# Run everything
# ============================================================

if __name__ == "__main__":
    # B_ROTATION = np.pi/4
    # plot_fields()

    plot_fields()
    plot_trajectories()
    plot_fields_with_trajectories()
    plot_beta_curves()
    run_parameter_sweep()
