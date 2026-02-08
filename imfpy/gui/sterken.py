#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Physical constants
# ===============================
AU = 1.495978707e11      # m
G = 6.67430e-11
M_sun = 1.98847e30
GM = G * M_sun

# Interstellar dust speed at infinity
v_inf = 26e3             # m/s

# ===============================
# Analytic beta-cone function
# ===============================
def beta_cone(beta, x):
    """
    Analytic beta-cone envelope from Sterken et al.
    x : distance along flow [m]
    returns z(x) [m]
    """
    mu_eff = GM * (1 - beta)
    a = -mu_eff / v_inf**2          # hyperbolic semi-major axis
    z = np.sqrt(np.maximum(0, -2 * a * x))
    return z

# ===============================
# Plot setup
# ===============================
plt.figure(figsize=(7, 7))

# x-range (Sterken scale)
x = np.linspace(-6, 6, 1000) * AU

# Beta values and colors (Sterken-like)
beta_vals = np.arange(1.1, 2.05, 0.1)
colors = plt.cm.rainbow(np.linspace(0, 1, len(beta_vals)))

# ===============================
# Plot beta-cones
# ===============================
for beta, c in zip(beta_vals, colors):
    z = beta_cone(beta, x)
    plt.plot(x / AU,  z / AU, ls=":", color=c, lw=2)
    plt.plot(x / AU, -z / AU, ls=":", color=c, lw=2)
    plt.text(-5.8, (beta - 1.0) * 4.5,
             f"β = {beta:.1f}", color=c, fontsize=9)

# ===============================
# Sun
# ===============================
plt.scatter(0, 0, s=120, color="gold", zorder=10)
plt.text(0.15, 0.1, "Sun", fontsize=10)

# ===============================
# Earth orbit
# ===============================
theta = np.linspace(0, 2*np.pi, 400)
plt.plot(np.cos(theta), np.sin(theta), 'k-', lw=1)
plt.text(1.05, -0.15, "Earth orbit", fontsize=9)

# ===============================
# Axes & labels
# ===============================
plt.xlabel("X [AU]  (along dust flow)")
plt.ylabel("Z [AU]")
plt.title("Interstellar Dust β-Cones (Sterken et al. style)")

plt.xlim(-6, 6)
plt.ylim(-6, 6)
plt.gca().set_aspect("equal", adjustable="box")
plt.grid(False)

plt.tight_layout()
plt.show()
