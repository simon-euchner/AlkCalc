################################################################################
### Python script for quickly plotting states                                ###
###                                                                          ###
### Author of this file: Simon Euchner                                       ###
################################################################################

# ---------------------------------------------------------------------------- #
### Packages
import numpy as np
from matplotlib import pyplot as plt
from pyalkcalc import plot_state
# ---------------------------------------------------------------------------- #

print("--- Input ------------")
species = input("Enter species: ")
print("Enter quantum numbers :")
n = int(input("n = ")); l = int(input("l = ")); j = float(input("j = "))
tmax = float(input("Enter maximal t: "))
N = int(input("Enter number of points: "))
ts = np.linspace(0, np.sqrt(tmax) - 1e-7, N)**2
plot_state(species, n, l, j, ts)
plt.show()
print("----------------------")
