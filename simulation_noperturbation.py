from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

# LEO orbit parameters
a = (6371 + 180) * u.km  # Earth radius + 180 km altitude
ecc = 0 * u.one  # circular orbit
inc = 0 * u.deg  # inclination
raan = 0 * u.deg  # right ascension of the ascending node
argp = 0 * u.deg  # argument of perigee
nu = 0 * u.deg  # true anomaly

# Creating the orbit
leo_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)

# Generate orbit data points
num_points = 100  # number of points to plot
period = leo_orbit.period
times = period * np.linspace(0, 1, num_points)

rr = [leo_orbit.propagate(time).r.to(u.km).value for time in times]  # Convert to km

# Extract x, y, z components and plot
x, y, z = zip(*rr)

# Plot the result
plt.figure(figsize=(8, 8))
plt.plot(x, y)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("LEO Orbit")
plt.grid(True)
plt.axis('equal')
plt.show()

