from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

# Define the LEO orbit parameters
a = (6371 + 180) * u.km  # semi-major axis
ecc = 0 * u.one  # eccentricity
inc = 0 * u.deg  # inclination
raan = 0 * u.deg  # right ascension of the ascending node
argp = 0 * u.deg  # argument of perigee
nu = 0 * u.deg  # true anomaly

# Create the orbit
leo_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)

# Define atmospheric drag perturbation function
def atmospheric_drag_acceleration(t0, state, k, R, C_D, A, m):
    from poliastro.constants import rho0_earth, H0_earth
    from poliastro.core.perturbations import atmospheric_drag_exponential
    from astropy.constants import R_earth
    from astropy import units as u
    
    # Convert state from km to m
    r = (state[:3] * u.km).to(u.m).value
    v = (state[3:] * u.km / u.s).to(u.m / u.s).value
    rho = rho0_earth.to(u.kg / u.m**3).value
    H0 = H0_earth.to(u.m).value
    R = R.to(u.m).value
    
    accel = atmospheric_drag_exponential(t0, r, v, k, R, C_D, A, m, rho, H0)
    
    # Convert acceleration back to km/s^2
    return np.array(accel) * u.m / u.s**2 to u.km / u.s**2

# Propagation parameters
C_D = 2.2  # Drag coefficient
A = 1 * u.m**2  # Cross-sectional area (m^2)
m = 100 * u.kg  # Mass (kg)

# Propagate the orbit including atmospheric drag
times = np.linspace(0, leo_orbit.period.to(u.s).value, num=1000) * u.s
rr = cowell(
    Earth.k,
    leo_orbit.r,
    leo_orbit.v,
    times,
    ad=atmospheric_drag_acceleration,
    R=Earth.R,
    C_D=C_D,
    A=A,
    m=m
)

# Extracting x, y positions for plotting
x, y, _ = rr[:, :3].T  # Transposing and ignoring z component

# Plotting
plt.figure(figsize=(10, 8))
plt.plot(x, y)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('LEO Orbit with Atmospheric Drag')
plt.grid(True)
plt.axis('equal')
plt.show()

