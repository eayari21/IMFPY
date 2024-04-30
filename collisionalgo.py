import os
import sys

import pandas as pd
import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

"""
Python script to process FORTRAN integrator data simulating charged particles
in the Interplanetary Magnetic Field (IMF).
__author__      = Ethan Ayari,
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust
"""


# %%CALL FORTRAN 95 SCRIPT
def fly_particle(q, m, pitch, x, y, z, KE):
    # Remove unwanted files from previous runs

    # Apply desired initial conditions
    script = open("earthcopy.f95", "r")
    scriptlines = script.readlines()
    scriptlines[10] = "       l_p=" + pitch + "\n"
    scriptlines[53] = "      l_p=" + pitch + "\n"
    scriptlines[11] = "       pm=" + m + "\n"
    scriptlines[9] = "       q=" + q + "\n"
    scriptlines[40] = "      y=" + y + "\n"
    scriptlines[41] = "      z=" + z + "\n"
    scriptlines[39] = "      x=" + x + "\n"
    scriptlines[43] = "      xe=" + KE + "\n"
    script = open("earthcopy.f95", "w")
    script.writelines(scriptlines)
    script.close()

    if(sys.platform == 'darwin'):
        # Run the sim
        gfortran = "/opt/local/bin/gfortran"
        os.system(gfortran+' earthcopy.f95')
        os.system("./a.out")
        data = pd.read_csv("subtest.8",
                           header=None, delim_whitespace=True)

    # if(sys.platform == 'win32'):
    #     gfortran = "C:\\MinGW\\bin\\gfortran.exe"
    #     os.system(gfortran+' earthcopy.f95')
    #     os.system("F:\\Terella\\fortranintegrator\\a.exe")
    #     data = pd.read_csv("F:\\Terella\\fortranintegrator\\subtest.8",
    #                        header=None, delim_whitespace=True)
    return data.transpose()


# %%PROCESS OUTPUT FROM FORTRAN 95 SCRIPT
def process_data(data, col, lab):
    data = data  # .transpose()
    # t = data.iloc[0]
    x = data.iloc[1]
    y = data.iloc[2]
    z = data.iloc[3]
    # vx = data.iloc[4]
    # vy = data.iloc[5]
    # vz = data.iloc[6]

    u = np.linspace(0, np.pi, 30)
    v = np.linspace(0, 2 * np.pi, 30)
    xsphere = 2.5*np.outer(np.sin(u), np.sin(v))
    ysphere = 2.5*np.outer(np.sin(u), np.cos(v))
    zsphere = 2.5*np.outer(np.cos(u), np.ones_like(v))

    # magv = np.sqrt(vx**2 + vy**2 + vz**2)
    """
    plt.plot((x+y)**2,z)
    plt.xlabel("Radius")
    plt.ylabel("Z")
    plt.show()
    """
    # mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = Axes3D(fig)
    # ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(xsphere, ysphere, zsphere, color='gray',
                      label='Aluminum (Al) Sphere')
    ax.plot(x, y, z, c='r')
    ax.set_xlabel('X', fontsize=20, fontweight='bold')
    ax.set_ylabel('Y', fontsize=20, fontweight='bold')
    ax.set_zlabel('Z', fontsize=20, fontweight='bold')
    ax.set_title(lab, fontsize=20, fontweight='bold')
    ax.legend()
    plt.show()


"""
    if os.path.exists("F:\\Terella\\fortranintegrator\\a.exe"):
        os.remove("F:\\Terella\\fortranintegrator\\a.exe")
    if os.path.exists("F:\\Terella\\fortranintegrator\\subtest.8"):
        os.remove("F:\\Terella\\fortranintegrator\\subtest.8")

"""


# %%CALL SIMULATION
# if __name__ == "__main__":
df1 = fly_particle("1", "1", "pi/2.", "6", "6", "1",
                   "1")

process_data(df1, 'm', 'Proton w/ angle of pi/2')

t = df1.iloc[0]
x = df1.iloc[1]
y = df1.iloc[2]
z = df1.iloc[3]
dataSet = np.array([x, y, z])
numDataPoints = len(t)

fig = plt.figure()
ax = plt.axes(projection='3d')


# %%ANIMATION FUNCTION FOR 3D DATA
def animate_func(num):
    ax.clear()  # Clears the figure to update the line, point, title, and axes
    # # Plot the Aluminum sphere
    # u = np.linspace(0, np.pi, 30)
    # v = np.linspace(0, 2 * np.pi, 30)
    # xsphere = 2.5*np.outer(np.sin(u), np.sin(v))
    # ysphere = 2.5*np.outer(np.sin(u), np.cos(v))
    # zsphere = 2.5*np.outer(np.cos(u), np.ones_like(v))
    # ax.plot_wireframe(xsphere, ysphere, zsphere, color='gray',
    #                   label='Aluminum (Al) Sphere')

    # Updating Trajectory Line (num+1 due to Python indexing)
    ax.plot3D(x[:num+1], y[:num+1], z[:num+1], c='red')

    # Updating Point Location
    ax.scatter(x[num], y[num], z[num],
               c='purple', marker='o', label="Electron")

    # Adding Constant Origin
    ax.plot3D(x[0], y[0], z[0], c='black', marker='o')

    # Setting Axes Limits
    ax.set_xlim3d([-8, 8])
    ax.set_ylim3d([-8, 8])
    ax.set_zlim3d([-3, 3])

    # Adding Figure Labels
    ax.set_title('Trajectory \nTime = ' + str(np.round(t[num],
                 decimals=7)) + ' sec', fontsize=20, fontweight='bold')
    ax.set_xlabel('X', fontsize=20, fontweight='bold')
    ax.set_ylabel('Y', fontsize=20, fontweight='bold')
    ax.set_zlabel('Z', fontsize=20, fontweight='bold')

    ax.legend()


line_ani = animation.FuncAnimation(fig, animate_func, interval=10,
                                   frames=numDataPoints)
plt.show()


# Saving the Animation
f = r"1_Videos/imfanimation.mp4"
writergif = animation.FFMpegWriter(fps=numDataPoints/60)
line_ani.save(f, writer=writergif)


"""
    df2 = fly_particle("4.8e-10","1.67e-24","pi/4.")
    process_data(df2,'y','Electron w/ pitch angle of pi/4')
    df3 = fly_particle("4.8e-10","1.67e-24","pi.")
    process_data(df3,'c','Electron w/ pitch angle of pi')
    df4 = fly_particle("-4.8e-10","9.1e-28","pi/2.")
    process_data(df4,'r','Electron w/ pitch angle of pi/2')
    df5 = fly_particle("-4.8e-10","9.1e-28","pi/6.")
    process_data(df5,'g','Electron w/ pitch angle of pi/6')
    """
