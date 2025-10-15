#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

"""
A Python script to describe the beta parameter for ISD grains.
__author__      = Ethan Ayari, 
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""

from scipy.stats import norm, lognorm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ||
# ||
# ||
plt.style.use("seaborn-v0_8-pastel")

# ||
# ||
# ||
def calculate_beta(L, Q_pr, G, M, c, rho, s, ideal=False):
    """
    Calculate beta using the given parameters.

    Parameters:
    - L: Constant parameter
    - Q_pr: Constant parameter
    - G: Gravitational constant
    - M: Mass of the body
    - c: Speed of light
    - rho: mass density
    - s: radius

    Returns:
    - beta: Calculated value using the provided equation
    """
    # Log-normal distribution for the Gaussian term
    stddev = 1.5
    # print(f"std = {stddev}")
    mean = .1
    dist=lognorm([stddev],loc=mean)

    # Modified formula using the log-normal term
    if ideal:
        return .5*rho/s**2
        # prop = G * M / s**2
        # beta = prop * (5.7e-5 * s**4/s / (rho * s))
    # rho * norm.pdf(s, loc=mean, scale=stddev)
    beta = ((5/rho)*dist.pdf(s))
    return beta


# ||
# ||
# ||
def turning_radius(G, B, M, v_i, phi):
    r = (4 * G * (B - 1) * M) / (v_i ** 2) * (1 + np.cos(phi))
    return r


# ||
# ||
# ||
def calculate_trajectory(G, B, M, v_i, phi, theta, initial_position):
    r_values = turning_radius(G, B, M, v_i, phi)
    x_values = initial_position[0] + r_values * np.sin(phi + theta) * np.cos(phi)
    y_values = initial_position[1] + r_values * np.sin(phi + theta) * np.sin(phi)
    z_values = initial_position[2] + r_values * np.cos(phi + theta)
    return x_values, y_values, z_values


# ||
# ||
# ||
def plot_deflected_trajectories(x_values, y_values, z_values, initial_position, heliosphere_radius):
    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

     # Get a color map from plt.cm.inferno
    color_map = plt.cm.inferno(np.linspace(0.3, 0.7, len(beta_values)))

    # Plot each trajectory separately with a different color
    for i, color in zip(range(len(beta_values)), color_map):
        ax.plot(x_values[i], y_values[i], z_values[i], color=color, label=r'$\beta$ = {:.2f}'.format(beta_values[i]))

    # # Plot Heliosphere sphere with radius 100
    # u = np.linspace(0, 2 * np.pi, 100)
    # v = np.linspace(0, np.pi, 100)
    # x_heliosphere = heliosphere_radius * np.outer(np.cos(u), np.sin(v))
    # y_heliosphere = heliosphere_radius * np.outer(np.sin(u), np.sin(v))
    # z_heliosphere = heliosphere_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    # ax.plot_surface(x_heliosphere, y_heliosphere, z_heliosphere, color='red', alpha=0.5)

    # Plot ISD sphere at the initial position
    ax.scatter(*initial_position, color='blue', label='ISD $r_{I}$')
    ax.text(*initial_position, ' ISD $r_{I}$', color='blue', fontsize=10, ha='left')



    # # Plot each trajectory separately
    # for i in range(len(beta_values)):
    #     r_values = np.sqrt(x_values[i]**2 + y_values[i]**2 + z_values[i]**2)

    #     # Plot points outside Heliosphere
    #     outside_heliosphere = r_values > heliosphere_radius
    #     ax.plot(x_values[i][outside_heliosphere], y_values[i][outside_heliosphere], z_values[i][outside_heliosphere], label=f'Beta = {beta_values[i]:.2f}')

    ax.set_xlabel('X [AU]', font="Times New Roman", fontsize=15, fontweight='bold')
    ax.set_ylabel('Y [AU]', font="Times New Roman", fontsize=15, fontweight='bold')
    ax.set_zlabel('Z [AU]', font="Times New Roman", fontsize=15, fontweight='bold')

    plt.title('Effect of Beta on the Trajectories of ISD Grains', fontsize=20, fontweight='bold')
    plt.legend()
    plt.show()


    # ||
    # ||
    # ||
if __name__ == "__main__":
    # Define range of beta values
    beta_values = np.linspace(0.1, 2.0, 8)  # Adjust the number of points as needed

    # Define other parameters
    G = 6.67430e-11  # Gravitational constant in m^3/kg/s^2
    B = beta_values
    M = 1e15  # Assume a typical mass for the sun in kg
    v_i = 26  # Initial velocity in m/s
    phi = np.linspace(0, np.pi, 100)  # Angle in radians
    # phi=0

    # Common initial position phase shift
    theta = np.radians(253)  # Convert degrees to radians
    phi_initial = np.radians(8)  # Convert degrees to radians
    initial_r = 150

    # Calculate initial position in Cartesian coordinates
    initial_position = (initial_r * np.sin(phi_initial) * np.cos(theta), initial_r * np.sin(phi_initial) * np.sin(theta), initial_r * np.cos(phi_initial))
    print(f"Initial position = {initial_position}")

    # Specify the Heliosphere radius
    heliosphere_radius = 100.0

    # Calculate trajectory for each beta value using the common initial position
    x_values, y_values, z_values = zip(*(calculate_trajectory(G, B_val, M, v_i, phi, theta, initial_position) for B_val in beta_values))

    # Visualize turning radii, plotting points outside the Heliosphere
    plot_deflected_trajectories(x_values, y_values, z_values, initial_position, heliosphere_radius)

    # ||
    # ||
    # ||
    # Generate values for s in the specified range
    s_values = np.logspace(-2, 1, 100)


    # Calculate beta values for different rho values
    beta_values_1 = [calculate_beta(1.0, 2, 6.67430e-11, 1.9e30, 3.0e8, 3.0, s_val) for s_val in s_values]
    beta_values_2 = [calculate_beta(1.0, 2, 6.67430e-11, 1.9e30, 3.0e8, 1.0, s_val, ideal=True) for s_val in s_values]
    beta_values_3 = [calculate_beta(1.0, 2, 6.67430e-11, 1.9e30, 3.0e8, 5.0, s_val) for s_val in s_values]

    print(f"1 = {beta_values_1}, 2 = {beta_values_2}, 3 = {beta_values_3}")


    # Create the plot
    plt.plot(s_values, beta_values_1, label=r'Graphite: $\rho$ = 2.26 g/$cm^{3}$', color='blue')
    plt.plot(s_values, beta_values_2, label=r'Ideal material: $\rho$ = 3 g/$cm^{3}$', color='red')
    plt.plot(s_values, beta_values_3, label=r'Olivine: $\rho$ = 3.22 g/$cm^{3}$', color='green')

    # Set logarithmic scale for both axes
    plt.xscale('log')
    plt.yscale('log')

    # Set x-axis limits
    plt.xlim(0.11, 10)

    # Set y-axis limits
    plt.ylim(0.01, 10)

    # Set labels and title with the specified parameters
    plt.xlabel(r's [$\mu$ m]', fontfamily='Times New Roman', fontsize=15, fontweight='bold')
    plt.ylabel(r'$\beta$', fontfamily='Times New Roman', fontsize=15, fontweight='bold')
    plt.title(r'$\beta$ vs. Particle Radius $s$', fontfamily='Times New Roman', fontsize=15, fontweight='bold')

    # Show the plot
    plt.legend(loc='upper right', fontsize=12, frameon=False, fancybox=True, edgecolor='black')
    plt.show()
