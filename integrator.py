#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

title = """

//  ================================================================================
//  ||                                                                            ||
//  ||              integrator                                                    ||
//  ||              ------------------------------------------------------        ||
//  ||                                 I N T E G R A T O R                        ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||                __author__      = Ethan Ayari                               ||
//  ||                IMPACT/LASP, CU Boulder                                     ||
//  ||                                                                            ||
//  ||                For: PHYS 7810 Plasma final project                         ||
//  ||                                                                            ||
//  ||                2024                                                        ||
//  ||                                                                            ||
//  ||                                                                            ||
//  ||                Works with Python 3.10.4                                    ||
//  ||                                                                            ||
//  ================================================================================


Integrate the equation of motion of charged particles in an electromagnetic field

"""
print(title)

"""
A Python script to visualize a simple flow.
__author__      = Ethan Ayari, 
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""

# ||
# ||
# ||
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
from scipy.integrate import odeint
from numba import njit

plt.style.use("seaborn-v0_8-pastel")

# || Lorentz force acting on a state vector
@njit
def lorentz_force(t, y, q, E, B):
    """
    Function to compute the derivative of the state vector y for the Lorentz force equation.
    
    Parameters:
        t (float): Time.
        y (array): State vector [x, y, z, vx, vy, vz].
        q (float): Charge of the particle.
        E (array): Electric field vector [Ex, Ey, Ez].
        B (array): Magnetic field vector [Bx, By, Bz].
        
    Returns:
        dydt (array): Derivative of the state vector.
    """
    # Unpack state vector
    x_, y_, z_, vx_, vy_, vz_ = y
    
    # Compute Lorentz force components
    Fx = q * (E[0] + vy_ * B[2] - vz_ * B[1])
    Fy = q * (E[1] + vz_ * B[0] - vx_ * B[2])
    Fz = q * (E[2] + vx_ * B[1] - vy_ * B[0])
    
    # Return derivative of the state vector
    return [vx_, vy_, vz_, Fx, Fy, Fz]

# || Executable Code
def integrate_lorentz_force(q, E, B, initial_conditions, t):
    """
    Function to integrate Newton's law for the Lorentz force equation.
    
    Parameters:
        q (float): Charge of the particle.
        E (array): Electric field vector [Ex, Ey, Ez].
        B (array): Magnetic field vector [Bx, By, Bz].
        initial_conditions (array): Initial conditions [x0, y0, z0, vx0, vy0, vz0].
        t (array): Time array.
        
    Returns:
        y (array): State vector at each time step.
    """
    # Integrate using odeint
    y = odeint(lorentz_force, initial_conditions, t, args=(q, E, B))
    return y


# || IMF
# %% Function to calculate the magnetic field components
def interplanetary_magnetic_field(x, y, U, R, B0, r0, phi0, theta, t):
    rsq = x**2 + y**2
    B_R = B0 * (r0 / rsq)**2 * np.cos(np.pi * t / 11 + phi0)
    B_T = B0 * (r0 / rsq) * np.cos(theta) * np.cos(np.pi * t / 11 + phi0)

    B_x = U * y + B_R * x / rsq - B_T * y / rsq
    B_y = -U * x + B_R * y / rsq + B_T * x / rsq

    return B_x, B_y


# || IMF Display
def display_interplanetary_magnetic_field(U, R, B0, r0, phi0, theta):

    # Create a grid for the streamlines within the circular region
    a = R  # Radius of the circle
    x = np.linspace(-a, a, 100)
    y = np.linspace(-a, a, 100)
    X, Y = np.meshgrid(x, y)

    # Create the figure and axes
    fig, ax = plt.subplots()
    color = np.log(1 + np.sqrt(interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0]**2 + interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1]**2))
    stream = ax.streamplot(X, Y, interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0], interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1],
                        color=color, linewidth=1, cmap=plt.cm.inferno, density=2, arrowstyle='->', arrowsize=1.5)

    # Create color normalization for the colormap
    norm = LogNorm(vmin=color.min(), vmax=color.max())

    # Create ScalarMappable for color mapping
    sm = plt.cm.ScalarMappable(cmap=plt.cm.inferno, norm=norm)
    sm.set_array([])

    # Add colorbar
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('IMF Magnitude')

    # Add time display using plt.figtext
    time_text = plt.figtext(0.925, 0.9, 'Time: 0 years', ha='center', va='center', color='black', fontsize=12)

    # Animation function
    def update(frame):
        ax.clear()
        # Recalculate magnetic field at each timestep
        B_x, B_y = interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)
        color = np.log(1+np.sqrt(B_x**2 + B_y**2))
        stream = ax.streamplot(X, Y, B_x, B_y,
                        color=color, linewidth=1, cmap=plt.cm.inferno, density=2, arrowstyle='->', arrowsize=1.5)
        time_text.set_text(f'Time: {frame} years')
        plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')

        return stream.lines, time_text

    # Create animation
    animation = FuncAnimation(fig, update, frames=np.arange(1, 100, 1), interval=25, blit=False)

    # Save the animation as an MP4 file
    animation.save('updated_parker_spiral.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

    # Show the plot
    plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')
    # plt.show()

# %%
# || Executable Code
if __name__ == "__main__":

    # || IMF Parameters
    U = 1.0
    R = 1.0
    B0 = 1.0
    r0 = 0.1
    phi0 = 0.0
    theta = np.pi / 4  # Example value for theta
    # display_interplanetary_magnetic_field(U, R, B0, r0, phi0, theta)

    # || Integrator Parameters

    # Define parameters
    q = 1  # Charge of the particle (e.g., electron)
    E = np.array([0, 0, 1])  # Electric field vector [Ex, Ey, Ez]
    B = np.array([0, 1, 0])  # Magnetic field vector [Bx, By, Bz]
    initial_conditions = np.array([0, 0, 0, 1, 1, 1])  # Initial conditions [x0, y0, z0, vx0, vy0, vz0]
    t = np.linspace(0, 10, 100)  # Time array
    
    # Integrate the Lorentz force equation
    y = integrate_lorentz_force(q, E, B, initial_conditions, t)
    
    # Print the final state vector
    print("Final state vector:")
    print(y[-1])

