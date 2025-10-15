#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

title = """

//  ================================================================================
//  ||                                                                            ||
//  ||              particle_flyer                                                ||
//  ||              ------------------------------------------------------        ||
//  ||                           P A R T I C L E   F L Y E R                      ||
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


GUI Interface for the particle flyer

"""
print(title)

"""
A Python script to visualize a simple flow.
__author__      = Ethan Ayari, 
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""


import sys
from PyQt6.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QVBoxLayout, QPushButton
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from .integrator import (
    display_interplanetary_magnetic_field,
    interplanetary_magnetic_field,
    integrate_lorentz_force,
    animate_state_vector_3d,
    animate_all,
)
from matplotlib.animation import FuncAnimation

animation_layout = QVBoxLayout()


class ParticleFlyer(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        vbox = QVBoxLayout()

        # Create labels and text boxes for IMF parameters
        parameters = ["U", "R", "B0", "r0", "phi0", "theta", "q", "m", "E", "Initial Conditions"]
        default_values = [1.0, 1.0, 1.0, 0.1, 0.0, 0.7853981634, 1, 1, "0,0,0", "0.5,0.5,-0.5,-1,1,0"]
        self.text_boxes = {}

        for param, default_value in zip(parameters, default_values):
            label = QLabel(param)
            textbox = QLineEdit()
            textbox.setText(str(default_value))
            vbox.addWidget(label)
            vbox.addWidget(textbox)
            self.text_boxes[param] = textbox

        # Create a button to start the calculation
        button = QPushButton("Start Calculation")
        button.clicked.connect(self.start_calculation)
        vbox.addWidget(button)

        self.setLayout(vbox)
        self.setWindowTitle('IMF Dust Trajectory Calculator')
        # Create the figure and axes
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)


    def clearLayout(self):
        for i in reversed(range(self.layout().count())):
            self.layout().itemAt(i).widget().setParent(None)

    def start_calculation(self):
        
        # Retrieve values from text boxes
        U = float(self.text_boxes["U"].text())
        R = float(self.text_boxes["R"].text())
        B0 = float(self.text_boxes["B0"].text())
        r0 = float(self.text_boxes["r0"].text())
        phi0 = float(self.text_boxes["phi0"].text())
        theta = float(self.text_boxes["theta"].text())
        q = float(self.text_boxes["q"].text())
        m = float(self.text_boxes["m"].text())

        E = np.array(list(map(float, self.text_boxes["E"].text().split(','))))
        initial_conditions = np.array(list(map(float, self.text_boxes["Initial Conditions"].text().split(','))))
        # Create a grid for the streamlines within the circular region
        a = R  # Radius of the circle
        x = np.linspace(-a, a, 100)
        y = np.linspace(-a, a, 100)
        X, Y = np.meshgrid(x, y)

        # Define time array
        t = np.linspace(0, 1, 100)

        # Integrate the Lorentz force equation
        y = integrate_lorentz_force(q, m, E, initial_conditions, t, U, R, B0, r0, phi0, theta)

        # Print the final state vector
        print("Final state vector:")
        print(y)


        color = np.log(1 + np.sqrt(interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0] ** 2 +
                                    interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1] ** 2))
        stream = self.ax.streamplot(X, Y, interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[0],
                            interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, 0)[1],
                            color=color, linewidth=1, cmap=plt.cm.inferno, density=2, arrowstyle='->', arrowsize=1.5)

        # Create color normalization for the colormap
        norm = LogNorm(vmin=color.min(), vmax=color.max())

        # Create ScalarMappable for color mapping
        sm = plt.cm.ScalarMappable(cmap=plt.cm.inferno, norm=norm)
        sm.set_array([])

        # Add colorbar
        cbar = plt.colorbar(sm, ax=self.ax)
        cbar.set_label('IMF Magnitude')

        # Add time display using plt.figtext
        time_text = plt.figtext(0.925, 0.9, 'Time: 0 years', ha='center', va='center', color='black', fontsize=12)

        # Animation function
        def update(frame):
            self.ax.clear()
            # Recalculate magnetic field at each timestep
            B_x, B_y = interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)[0], \
                    interplanetary_magnetic_field(X, Y, U, R, B0, r0, phi0, theta, frame)[1]
            color = np.log(1 + np.sqrt(B_x ** 2 + B_y ** 2))
            stream = self.ax.streamplot(X, Y, B_x, B_y,
                                color=color, linewidth=1, cmap=plt.cm.inferno, density=2, arrowstyle='->',
                                arrowsize=1.5)
            time_text.set_text(f'Time: {frame} years')
            plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
            plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
            plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')

            # Plot state vector
            plt.plot(y[:frame, 0], y[:frame, 1], 'r-', alpha=0.5)

            return stream.lines, time_text

        # Create animation
        animation = FuncAnimation(self.fig, update, frames=np.arange(1, y.shape[0], 1), interval=25, blit=False)

        # Show the plot
        plt.xlabel('X [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.ylabel('Y [100*AU]', font="Times New Roman", fontsize=20, fontweight='bold')
        plt.title('Parker Spiral B-field Model', font="Times New Roman", fontsize=20, fontweight='bold')
        # Clear existing layout

        # Display the animation
        # Create a new window for the animation
        

        self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    ex = ParticleFlyer()
    ex.show()
    animation_window = QWidget()
    animation_layout.addWidget(ex.canvas)
    animation_layout.addWidget(ex.toolbar)
    animation_window.setLayout(animation_layout)
    animation_window.setWindowTitle('Particle Animation')
    animation_window.show()
    sys.exit(app.exec())
