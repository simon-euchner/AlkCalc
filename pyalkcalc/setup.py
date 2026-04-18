from setuptools import Extension, find_packages, setup
from Cython.Build import cythonize

extensions = [
    Extension(
        name="pyalkcalc._core",
        sources=["src/pyalkcalc/_core.pyx"],
        libraries=["alkcalc"],
    )
]

setup(
    name="pyalkcalc",
    version="0.1.0",
    package_dir={"": "src"},
    packages=find_packages("src"),
    ext_modules=cythonize(extensions, language_level=3),
)
