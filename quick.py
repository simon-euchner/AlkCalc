################################################################################
### Python script for quickly calculating electronic properties              ###
###                                                                          ###
### Author of this file: Simon Euchner                                       ###
################################################################################

# ---------------------------------------------------------------------------- #
### Packages
from pyalkcalc import (
        energy,
        radial_matrix_element,
        oscillator_strength,
        lifetime,
)
# ---------------------------------------------------------------------------- #

print("--- Input ------------")
species = input("Enter species: ")
p = float(input("Enter power p of radial matrix element <a|r**p|b>: "))
print("Enter quantum numbers for |a> = |na,la,ja> and |b> = |nb,lb,jb>:")
na = int(input("na = ")); la = int(input("la = ")); ja = float(input("ja = "))
nb = int(input("nb = ")); lb = int(input("lb = ")); jb = float(input("jb = "))
print("--- Results ----------")
Ea = energy(species, na, la, ja)
Eb = energy(species, nb, lb, jb)
rp = radial_matrix_element(species, na, la, ja, p, nb, lb, jb)
f = oscillator_strength(species, na, la, ja, nb, lb, jb)
print(
        "Energies Ea and Eb of the states |a> and |b>:\n"
        "Ea = {:1.6E} Hartree\n"
        "Eb = {:1.6E} Hartree\n"
        "Eb - Ea = {:1.6E} Hartree\n"
        "Radial matrix element:\n"
        "<a|r**p|b> = {:f} aB**p\n"
        "Oscillator:\n"
        "f = {:f}\n".format(Ea, Eb, Eb - Ea, rp, f)
)
lfts = str(input("Compute lifetime of |a>? (y/N)?: "))
if lfts == "y" or lfts == "Y":
    tauRAD = lifetime(0, species, na, 0, la, ja)
    tauBBR = lifetime(300, species, na, 30, la, ja)
    print(
            "Radiative Lifetime of |a>:\n"
            "tau = {:f} ns\n"
            "Lifetime of |a> at room temperature:\n"
            "tau = {:f} ns".format(tauRAD, tauBBR)
    )
    print("----------------------")
else:
    print("----------------------")
