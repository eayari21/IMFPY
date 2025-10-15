#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

title = """

//  ================================================================================
//  ||                                                                            ||
//  ||              vector_field                                                  ||
//  ||              ------------------------------------------------------        ||
//  ||                           V E C T O R  F I E L D                           ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||                __author__      = Ethan Ayari                               ||
//  ||                IMPACT/LASP, CU Boulder                                     ||
//  ||                                                                            ||
//  ||                For: NASA FINESST 2024                                      ||
//  ||                                                                            ||
//  ||                2024                                                        ||
//  ||                                                                            ||
//  ||                                                                            ||
//  ||                Works with Python 3.10.4                                    ||
//  ||                                                                            ||
//  ================================================================================


Visualize vector fields relevant for the IMF and interstellar magnetic field.

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

# Global variable to store streamplot lines
stream = None
cbar = None


# ||
# ||
# ||
# %%Function to calculate the vector field V(R)
def vector_field(R, a, v0):
    r = np.linalg.norm(R)
    V = (a * R) / (r**3) - v0
    return V


# ||
# ||
# ||
# %%
def create_vector_field(a, v0):
    # Generate 3D grid of points
    x, y, z = np.meshgrid(np.linspace(-200, 200, 100), np.linspace(-200, 200, 100), np.linspace(-200, 200, 100))

    

    # Calculate the vector field at each point in the grid
    V_x, V_y, V_z = np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            for k in range(x.shape[2]):
                R = np.array([x[i, j, k], y[i, j, k], z[i, j, k]])
                V = vector_field(R, a, v0)
                V_x[i, j, k], V_y[i, j, k], V_z[i, j, k] = V

    return x, y, z, V_x, V_y, V_z


# ||
# ||
# ||
# %%
def visualize_vector_field(x, y, z, V_x, V_y, V_z):
    # Create a 3D quiver plot
    fig = plt.figure(figsize=(12, 5))

    # Plot 1: 3D quiver plot
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.quiver(x, y, z, V_x, V_y, V_z, length=.5, normalize=True, color='red', arrow_length_ratio=0.5)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title('Vector Field in 3D')

    # Plot 2: Cross-section in the x-y plane (z=0)
    ax2 = fig.add_subplot(122)
    ax2.quiver(x[:, :, -1], y[:, -1, :], V_x[:, :, -1], V_y[:, :, -1], scale=20, color='red', scale_units='xy', angles='xy')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_title('Cross-section in X-Y Plane (Z=0)')

    # Show the plots
    plt.tight_layout()
    plt.show()


# ||
# ||
# ||
# %%
def vector_field_2d(X, Y, a, v0):
    R = np.array([X, Y, 0])  # Assume z=0 for a 2D visualization
    r = np.linalg.norm(R)
    V = (a * R) / (r**3) - v0
    return V[:2]  # Take only the x and y components


# ||
# ||
# ||
# %%
def visualize_vector_field_2d(a, v0, xlim=(-2, 2), ylim=(-2, 2), quiver_scale=20):
    # Generate a grid of points
    x = np.linspace(xlim[0], xlim[1], 20)
    y = np.linspace(ylim[0], ylim[1], 20)
    X, Y = np.meshgrid(x, y)

    # Calculate the vector field at each point in the grid
    V_x, V_y = np.zeros_like(X), np.zeros_like(Y)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            V = vector_field_2d(X[i, j], Y[i, j], a, v0)
            V_x[i, j], V_y[i, j] = V

    # Visualize the 2D vector field
    plt.quiver(X, Y, V_x, V_y, scale=quiver_scale, color='blue', scale_units='xy', angles='xy')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('2D Vector Field Visualization')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.grid(True)
    plt.show()

# || IMF
# %% Function to calculate the magnetic field components
def magnetic_field(x, y, U, R, B0, r0, phi0, theta, t):
    rsq = x**2 + y**2
    B_R = B0 * (r0 / rsq)**2 * np.cos(np.pi * t / 11 + phi0)
    B_T = B0 * (r0 / rsq) * np.cos(theta) * np.cos(np.pi * t / 11 + phi0)
    
    B_x = U * y + B_R * x / rsq - B_T * y / rsq
    B_y = -U * x + B_R * y / rsq + B_T * x / rsq
    
    return B_x, B_y

# ||
# ||
# || INTERSTELLAR B-FIELD
# %%Function to calculate the stream function for flow over a cylinder
def stream_function(x, y, U, R):
    rsq = x ** 2 + y ** 2
    psi = U * x * ( 1.0 - R ** 2 / rsq )
    return np.where( rsq >= R * R, psi, 0.0 )

# velocity compoents
def velocity_components(x, y, U, R):
    rsq = x ** 2 + y ** 2
    Rsq = R ** 2
    u = U * ( 1 - Rsq / rsq + 2 * Rsq * y ** 2 / ( rsq * rsq ) )
    v = -2 * U * Rsq * x * y / ( rsq * rsq )
    u = np.where( rsq >= Rsq, u, 0.0 )
    v = np.where( rsq >= Rsq, v, 0.0 )
    return u, v

def plot_streamlines(ax, X, Y, stream_function, U, R):
    psi = stream_function(X, Y, U, R)
    levels = np.linspace(np.min(psi), np.max(psi), 50)
    # ax.contour(X, Y, psi, levels=levels, colors='b')
    
    V_x = partial_derivative(stream_function, X, Y, U, R, direction=0)
    V_y = partial_derivative(stream_function, X, Y, U, R, direction=1)

    color = np.sqrt(V_x**2 + V_y**2)
    
    stream = ax.streamplot(X, Y, -V_x, V_y, color=color, linewidth=1, cmap=plt.cm.inferno,
                            density=2, arrowstyle='->', arrowsize=1.5)
    
    return stream

def partial_derivative(f, X, Y, U, R, direction):
    h = 1e-5  # small perturbation
    if direction == 0:
        return (f(X, Y + h, U, R) - f(X, Y - h, U, R)) / (2 * h)
    elif direction == 1:
        return (f(X + h, Y, U, R) - f(X - h, Y, U, R)) / (2 * h)
    else:
        raise ValueError("Invalid direction. Use 0 for X-direction or 1 for Y-direction.")

    

# ||
# %% Function to plot velocity vectors
def plot_velocity_vectors(ax, X, Y, velocity_components, U, R ):
    spacing = 5
    x_points = X[::spacing, ::spacing]
    y_points = Y[::spacing, ::spacing]

    psi_x, psi_y = velocity_components(x_points, y_points, U, R)
    ax.quiver(x_points, y_points, psi_x, psi_y, color='r', scale=0.3 )


# ||
# ||
# ||
 # %%
def field_line_equations(y, t, a, v0):
    # Define the system of ordinary differential equations for field lines
    x, y = y
    V = vector_field_2d(x, y, a, v0)
    dxdt, dydt = V[0], V[1]
    return [dxdt, dydt]


# ||
# ||
# ||
# %%
def visualize_field_lines(a, v0, fig, ax, xlim=(-200, 200), ylim=(-200, 200), num_lines=10, show=False):
    global cbar
    # Add a labeled colorbar
    if cbar:
        ax.clear
        # plt.delaxes(cbar.ax)

    ax.clear()
    # Generate a grid of points
    x = np.linspace(xlim[0], xlim[1], 200)
    y = np.linspace(ylim[0], ylim[1], 200)
    X, Y = np.meshgrid(x, y)

    # Calculate the vector field at each point in the grid
    V_x, V_y = np.zeros_like(X), np.zeros_like(Y)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            V = vector_field_2d(X[i, j], Y[i, j], a, v0)
            V_x[i, j], V_y[i, j] = V

    color = 2 * np.log(np.hypot(V_x, V_y))
    # Visualize the 2D vector field
    # Clear previous streamplot lines
    
    stream = ax.streamplot(X, Y, V_x, V_y,  color=color, linewidth=1, cmap=plt.cm.inferno,
              density=2, arrowstyle='->', arrowsize=1.5)
    

    
    # Add a circle to the specific 'ax' axes
    ionocircle = plt.Circle((0, 0), a, color='red', fill=False, linewidth=3, label="Ionopause")
    nucleocircle = plt.Circle((0, 0), .1*a, color='black', fill=True, label="Nucleopause")
    # ax.add_patch(ionocircle)
    # ax.add_patch(nucleocircle)


    
    if cbar:
        pass
    else:
        cbar = plt.colorbar(stream.lines, ax=ax)
        cbar.set_label('Vector Field Magnitude')

    ax.set_xlabel('X [AU]', font="Times New Roman", fontsize=15, fontweight='bold')
    ax.set_ylabel('Y [AU]', font="Times New Roman", fontsize=15, fontweight='bold')
    plt.suptitle(r'Heliospheric Flow Lines $V = \frac{(\alpha * r)}{(r^{3})} - v0$', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.xlim(xlim)
    plt.ylim(ylim)
    ax.grid(True)
    if show:
        plt.show()
    return stream

# ||
# %%
def update(val):
        a = slider_a.val
        v0[0] = slider_v0_x.val
        stream = visualize_field_lines(a, v0, fig, ax)

        plt.show()


# %%
# || Executable Code
if __name__ == "__main__":
    # || Heliospheric/ISM dust flowlines
    # Set values for a and v0
    # a = 20000
    # v0 = np.array([-1.5, 0.0, 0.0])  # v0 is a 3D vector
    # x, y, z, V_x, V_y, V_z = create_vector_field(a, v0)
    # # visualize_vector_field(x, y, z, V_x, V_y, V_z)

    # # Visualize the 2D vector field
    # # visualize_vector_field_2d(a, v0)

    # # Visualize the 2D vector field and field lines
    # # visualize_field_lines(a, v0, show=True)

    # # Create sliders using matplotlib.widgets
    # fig, ax = plt.subplots(figsize=(5, 2))
    # plt.subplots_adjust(bottom=0.25)

    # # ax_a = plt.axes([0.1, 0.1, 0.65, 0.03])
    # # slider_a = Slider(ax_a, r'$\alpha$', 0, 2.0, valinit=a)

    # # ax_v0 = plt.axes([0.1, 0.05, 0.65, 0.03])
    # # slider_v0_x = Slider(ax_v0, 'v0', -1.0, 1.0, valinit=-1.0)


    # visualize_field_lines(a, v0, fig, ax)
    # # # Update the plot when sliders are changed

    # # slider_a.on_changed(update)
    # # slider_v0_x.on_changed(update)

    # plt.show()

    # # || VLISM B-FIELD CALCS:

    # # Create a figure and axis
    # fig, ax = plt.subplots()

    # # Cylinder properties
    # D = 150  # Diameter of the cylinder
    # R = D / 2  # Radius of the cylinder
    # U = 0.01  # Fluid velocity

    # # Set plot limits
    # ax.set_xlim(-D, D)
    # ax.set_ylim(-D, D)

    # # Set aspect ratio to be equal
    # ax.set_aspect('equal')

    # # Create a grid for the streamlines
    # x = np.linspace( -2*D, 2*D, 200 )
    # y = np.linspace( -2*D, 2*D, 200 )
    # X, Y = np.meshgrid(x, y)

    # stream = plot_streamlines( ax, X, Y, stream_function, U, R )
    # # plot_velocity_vectors( ax, X, Y, velocity_components, U, R )
        
    # ax.set_xlim(-2*D, 2*D)
    # ax.set_ylim(-2*D, 2*D)
    # ax.set_aspect('equal')
    # cbar = plt.colorbar(stream.lines, ax=ax)
    # cbar.set_label('VLISM B-Field Magnitude')

    # ax.set_xlabel('X [AU]', font="Times New Roman", fontsize=15, fontweight='bold')
    # ax.set_ylabel('Y [AU]', font="Times New Roman", fontsize=15, fontweight='bold')
    # plt.suptitle("VLISM Magnetic Field Approximation", font="Times New Roman", fontsize=20, fontweight='bold')

    # plt.show()

    # || IMF Display

   # Parameters
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
    color = np.log(1 + np.sqrt(magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0]**2 + magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1]**2))
    stream = ax.streamplot(X, Y, magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0], magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1],
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
        B_x, B_y = magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)
        color = np.log(1+np.sqrt(B_x**2 + B_y**2))
        stream = ax.streamplot(X, Y, B_x, B_y,
                        color=color, linewidth=1, cmap=plt.cm.inferno, density=2, arrowstyle='->', arrowsize=1.5)
        time_text.set_text(f'Time: {frame} years')
        plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')

        return stream.lines, time_text

    # Create animation
    animation = FuncAnimation(fig, update, frames=np.arange(1, 100, 1), interval=100, blit=False)

    # Save the animation as an MP4 file
    # animation.save('updated_parker_spiral.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

    # Show the plot
    plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')
    plt.show()