import numpy as np
import matplotlib.pyplot as plt

AU = 1.495978707e11
G = 6.67430e-11
M_sun = 1.98847e30
GM = G * M_sun


def accel(r, beta):
    norm = np.linalg.norm(r)
    return -(1 - beta) * GM * r / norm**3


# Initial conditions
r = np.array([3 * AU, 0.0])
v = np.array([0.0, -26e3])

beta = 0.3
dt = 1800.0  # 30 min
steps = 20000

xs, ys = [], []

for _ in range(steps):
    a = accel(r, beta)
    v += a * dt
    r += v * dt

    xs.append(r[0] / AU)
    ys.append(r[1] / AU)

plt.figure(figsize=(6,6))
plt.plot(xs, ys)
plt.scatter(0, 0, color="orange", label="Sun")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.axis("equal")
plt.legend()
plt.title(f"Radiation + Gravity Trajectory (β = {beta})")
plt.show()
