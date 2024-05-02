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
A Python script to integrate particle trajectories along a given path.
__author__      = Ethan Ayari, 
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""

# ||
# ||
# ||
import subprocess
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
# @njit
def lorentz_force(y, t, q, m, E, U, R, B0, r0, phi0, theta):
    x, y, z, vx, vy, vz = y
    
    # Compute magnetic field components at the current location
    B_x, B_y, B_z = interplanetary_magnetic_field(x, y, U, R, B0, r0, phi0, theta, t)
    
    # Compute Lorentz force components
    Fx = q/m * (E[0] + vy * B_z - vz * B_y)
    Fy = q/m * (E[1] + vz * B_x - vx * B_z)
    Fz = q/m * (E[2] + vx * B_y - vy * B_x)
    
    # Return derivative of the state vector
    return [vx, vy, vz, Fx, Fy, Fz]


# || Integrator
def integrate_lorentz_force(q, m, E, initial_conditions, t, U, R, B0, r0, phi0, theta):
    y = np.zeros((len(t), 6))
    y[0] = initial_conditions  # Set initial conditions
    
    for i in range(1, len(t)):
        # Compute magnetic field components at the current location
        B_x, B_y, B_z = interplanetary_magnetic_field(y[i-1, 0], y[i-1, 1], U, R, B0, r0, phi0, theta, t[i])
        print(f"B = {B_x}, {B_y}, {B_z}")
        
        # Convert y[i-1] to a list
        y_im1_list = y[i-1].tolist()
        
        # Integrate using Euler's method
        dt = t[i] - t[i-1]

        formatted_string = f"y[i] = {y[i]}\n" \
                   f"y[i-1] = {y[i-1]}\n" \
                   f"dt = {dt}\n" \
                   f"t[i-1] = {t[i-1]}\n" \
                   f"q = {q}\n" \
                   f"E = {E}\n" \
                   f"U = {U}\n" \
                   f"R = {R}\n" \
                   f"B0 = {B0}\n" \
                   f"r0 = {r0}\n" \
                   f"phi0 = {phi0}\n" \
                   f"theta = {theta}\n"

        # print(formatted_string)


        y[i] = np.array(y_im1_list) + dt * np.array(lorentz_force(y_im1_list, t[i-1], q, m, E, U, R, B0, r0, phi0, theta))
    
        
    
    return y


# || IMF
# %% Function to calculate the magnetic field components
def interplanetary_magnetic_field(x, y, U, R, B0, r0, phi0, theta, t):
    rsq = x**2 + y**2

    B_R = B0 * (r0 / rsq)**2 * np.cos(np.pi * t / 11 + phi0)
    B_T = B0 * (r0 / rsq) * np.cos(theta) * np.cos(np.pi * t / 11 + phi0)

    B_x = U * y + B_R * x / rsq - B_T * y / rsq
    B_y = -U * x + B_R * y / rsq + B_T * x / rsq

    return B_x, B_y, 0

    


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
        B_x, B_y = interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)[0], interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)[1]
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
    # animation.save('updated_parker_spiral.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

    # Show the plot
    plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.show()
    return animation

# || Animate state vector in 3d
def animate_state_vector_3d(state_vector, timestamps):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set initial plot data
    line, = ax.plot([], [], [], lw=2)
    ax.set_xlim3d(np.min(state_vector[:, 0]), np.max(state_vector[:, 0]))
    ax.set_ylim3d(np.min(state_vector[:, 1]), np.max(state_vector[:, 1]))
    ax.set_zlim3d(np.min(state_vector[:, 2]), np.max(state_vector[:, 2]))

    # Function to update the plot for each frame
    def update(frame):
        line.set_data(state_vector[:frame, 0], state_vector[:frame, 1])
        line.set_3d_properties(state_vector[:frame, 2])
        # ax.set_title(f'Time: {timestamps[frame]:.2f}', fontsize=20, fontweight='bold', fontname='Times New Roman')
        ax.set_xlabel('X', fontsize=20, fontweight='bold', fontname='Times New Roman')
        ax.set_ylabel('Y', fontsize=20, fontweight='bold', fontname='Times New Roman')
        ax.set_zlabel('Z', fontsize=20, fontweight='bold', fontname='Times New Roman')
        return line,

    # Create animation
    ani = FuncAnimation(fig, update, frames=len(timestamps), interval=10, blit=True)
    ani.save('spiral_w_trajectory.mp4', fps=10, extra_args=['-vcodec', 'libx264'])
    
    # plt.show()

def animate_state_vector_xy(state_vector, timestamps):
    fig, ax = plt.subplots()

    # Set initial plot data
    line, = ax.plot([], [], lw=2)
    # ax.set_xlim(np.min(state_vector[:, 0]), np.max(state_vector[:, 0]))
    # ax.set_ylim(np.min(state_vector[:, 1]), np.max(state_vector[:, 1]))
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    # Function to update the plot for each frame
    def update(frame):
        line.set_data(state_vector[:frame, 0], state_vector[:frame, 1])
        # ax.set_title(f'Time: {timestamps[frame]:.2f}', fontsize=20, fontweight='bold', fontname='Times New Roman')
        ax.set_xlabel('X', fontsize=20, fontweight='bold', fontname='Times New Roman')
        ax.set_ylabel('Y', fontsize=20, fontweight='bold', fontname='Times New Roman')
        return line,

    # Create animation
    ani = FuncAnimation(fig, update, frames=len(timestamps), interval=50, blit=True)
    
    plt.show()
    return ani

# || Combine the animations using ffmpeg
def combine_animations(particle_animation, magnetic_field_animation):
    # Save the animations as temporary files
    particle_animation.save('particle_animation.mp4', fps=10, extra_args=['-vcodec', 'libx264'])
    magnetic_field_animation.save('magnetic_field_animation.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

    # Combine the animations using ffmpeg
    command = 'ffmpeg -i magnetic_field_animation.mp4 -i particle_animation.mp4 -filter_complex overlay=0:0 all_animations.mp4'
    subprocess.run(command, shell=True)

# || Superimpose particle trajectory on magnetic field
def animate_all(state_vector):
    U = 1.0
    R = 1.0
    B0 = 1.0
    r0 = 0.1
    phi0 = 0.0
    theta = np.pi / 4  # Example value for theta

    # Create a grid for the streamlines within the circular region
    a = R  # Radius of the circle
    x = np.linspace(-a, a, 100)
    y = np.linspace(-a, a, 100)
    X, Y = np.meshgrid(x, y)

    # Create the figure and axes
    fig, ax = plt.subplots()
    color = np.log(1 + np.sqrt(interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0] ** 2 +
                                interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1] ** 2))
    stream = ax.streamplot(X, Y, interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0],
                           interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1],
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
        B_x, B_y = interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)[0], \
                   interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)[1]
        color = np.log(1 + np.sqrt(B_x ** 2 + B_y ** 2))
        stream = ax.streamplot(X, Y, B_x, B_y,
                               color=color, linewidth=1, cmap=plt.cm.inferno, density=2, arrowstyle='->',
                               arrowsize=1.5)
        time_text.set_text(f'Time: {frame} years')
        plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')

        # Plot state vector
        plt.plot(state_vector[:frame, 0], state_vector[:frame, 1], 'r-', alpha=0.5)

        return stream.lines, time_text

    # Create animation
    animation = FuncAnimation(fig, update, frames=np.arange(1, state_vector.shape[0], 1), interval=25, blit=False)

    # Show the plot
    plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')
    # plt.show()
    animation.save('spiral_w_trajectory.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

    return animation




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
    # magnetic_field_animation = display_interplanetary_magnetic_field(U, R, B0, r0, phi0, theta)

    # || Integrator Parameters

    # Define parameters
    q = 1  # Charge of the particle (e.g., electron)
    m = 1 
    E = np.array([0, 0, 0])  # Electric field vector [Ex, Ey, Ez]

    initial_conditions = np.array([.5, .5, -.5, -1, 1, 0])  # Initial conditions [x0, y0, z0, vx0, vy0, vz0]
    t = np.linspace(0, 1, 100)  # Time array

    # Integrate the Lorentz force equation
    y = integrate_lorentz_force(q, m, E, initial_conditions, t, U, R, B0, r0, phi0, theta)

    # Print the final state vector
    print("Final state vector:")
    print(y)

    timestamps = np.linspace(0, 1, len(y))

    # Animate the state vector in both 3d and 2d
    animate_state_vector_3d(y, timestamps)
    # particle_animation = animate_state_vector_xy(y[:, :2], timestamps)
    animate_all(y)

    # Overlay particle trajectory onto spiral
    # magnetic_field_animation = display_interplanetary_magnetic_field(U, R, B0, r0, phi0, theta)
    # combine_animations(particle_animation, magnetic_field_animation)
