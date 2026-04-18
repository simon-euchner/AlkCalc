################################################################################
### Setup                                                                    ###
###                                                                          ###
### Author of this file: Simon Euchner                                       ###
################################################################################

# ---------------------------------------------------------------------------- #
### Packages
from setuptools import setup, Extension
from Cython.Build import cythonize
# ---------------------------------------------------------------------------- #

extension = Extension(
        name = "PyAlkCalc._alkcalc",
        sources = ["src/PyAlkCalc/_alkcalc.pyx"],
        include_dirs = ["/home/simon/Files/GitHub/AlkCalc/interface"],
        library_dirs = ["/home/simon/Files/GitHub/AlkCalc/lib"],
        libraries = ["alkcalc"],
        extra_link_args=[
            "-Wl,-rpath,/home/simon/Files/GitHub/AlkCalc/lib"
        ]
    )

setup(
    ext_modules=cythonize([extension], language_level=3),
)
