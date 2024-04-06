from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.plotting import StaticOrbitPlotter
from astropy import units as u

# LEO orbit parameters
a = (6371 + 180) * u.km  # Earth radius + 180 km altitude
ecc = 0 * u.one  # circular orbit
inc = 0 * u.deg  # inclination, for example purposes
raan = 0 * u.deg  # arbitrary for this example
argp = 0 * u.deg  # not significant for circular orbits
nu = 0 * u.deg  # arbitrary for circular orbits

# Creating the orbit
leo_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)

# Plotting the orbit
plotter = StaticOrbitPlotter()
plotter.plot(leo_orbit)

