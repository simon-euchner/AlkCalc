################################################################################
### Setup                                                                    ###
###                                                                          ###
### Author of this file: Simon Euchner                                       ###
################################################################################

# ---------------------------------------------------------------------------- #
### Packages
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
# ---------------------------------------------------------------------------- #

extension = Extension(
        name = "PyAlkCalc._alkcalc",
        sources = ["src/PyAlkCalc/_alkcalc.pyx"],
        include_dirs = [numpy.get_include()],
        libraries = ["alkcalc"],
        library_dirs = ["/home/simon/Files/GitHub/Alkcalc/lib"],

    )

setup(
    ext_modules=cythonize([extension], language_level=3),
)
