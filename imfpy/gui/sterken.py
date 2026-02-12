#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# Constants
# ==========================================================
AU = 1.495978707e11
G  = 6.67430e-11
M_sun = 1.98847e30
GM = G * M_sun
v_inf = 26e3

# flow tilt from paper (~7 degrees in Z)
tilt = np.deg2rad(7.0)

# ==========================================================
# Beta paraboloid
# ==========================================================
def beta_paraboloid(beta, z):
    a = 2.0 * GM * (beta - 1.0) / v_inf**2
    return a - (z**2) / (2.0 * a)

# ==========================================================
# Setup
# ==========================================================
beta_vals = np.arange(1.1, 2.05, 0.1)
colors = plt.cm.rainbow(np.linspace(0, 1, len(beta_vals)))

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# ==========================================================
# LEFT PANEL — XZ PLANE (THIS IS THE KEY FIX)
# ==========================================================
ax = axs[0]

z = np.linspace(-6, 6, 2000) * AU

for beta, c in zip(beta_vals, colors):
    x = beta_paraboloid(beta, z)

    # apply flow tilt
    z_tilt = z * np.cos(tilt)
    x_tilt = x + z * np.sin(tilt)

    ax.plot(x_tilt/AU, z_tilt/AU, ':', color=c, lw=2)

    ax.text(3.75, (beta-1.0)*5.0, f'β = {beta:.1f}', color=c, fontsize=15)

# Earth orbit (projection into XZ)
theta = np.linspace(0, 2*np.pi, 600)
ax.plot(np.cos(theta), np.sin(theta), color='blue', lw=1)
ax.text(0.2, -1.35, 'Earth orbit', color='blue', fontsize = 15)

# Sun
ax.scatter(0, 0, s=120, color='gold', zorder=5)
ax.text(-.9, 0.1, 'Sun', fontsize=14)

ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_aspect('equal')
ax.set_xlabel('X [AU]')
ax.set_ylabel('Z [AU]')
ax.set_title('XZ plane (along ISD flow)')

# ==========================================================
# RIGHT PANEL — YZ CROSS SECTION
# ==========================================================
ax = axs[1]

x_cut = 5.0 * AU   # upstream slice as in paper
phi = np.linspace(0, 2*np.pi, 800)

for beta, c in zip(beta_vals, colors):
    a = 2.0 * GM * (beta - 1.0) / v_inf**2
    r = np.sqrt(2 * a * x_cut) / AU

    ax.plot(r*np.cos(phi), r*np.sin(phi), ':', color=c, lw=2)

    ax.text(3.75, 5.0 * (beta - 1.0),
            f'β = {beta:.1f}',
            color=c, fontsize=15)    

# Earth orbit in ecliptic frame
x = np.cos(phi)
y = np.sin(phi)
z = np.zeros_like(phi)

# rotate orbit about Y-axis by ISD tilt
z_rot = -x * np.sin(.1*tilt)

# YZ projection
# rotate in Y–Z plane to match Sterken orientation
theta = np.deg2rad(-7)   # ~10–20° works; Sterken is ~15°

y_final =  y * np.cos(theta) - z_rot * np.sin(theta)
z_final =  y * np.sin(theta) + z_rot * np.cos(theta)

ax.plot(y_final, z_final, color='blue', lw=1)
ax.text(-1.2, -.6, 'Earth orbit', color='blue', fontsize=15)





# Sun
ax.scatter(0, 0, s=120, color='gold', zorder=5)
ax.text(-0.9, 0.15, 'Sun', fontsize=14)

ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_aspect('equal')
ax.set_xlabel('Y [AU]')
ax.set_ylabel('Z [AU]')
ax.set_title('YZ plane (cross section upstream)')

plt.tight_layout()
plt.show()
