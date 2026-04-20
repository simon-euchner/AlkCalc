# pyalkcalc

Python interface to [*AlkCalc*](https://github.com/simon-euchner/AlkCalc) — an
alkali-metal atom and alkaline-earth-metal ion calculator written in C.

**Author:** Simon Euchner [euchner.se[at]gmail.com]

**License:** GPLv3


## Overview

The Python package *pyalkcalc* exposes the *AlkCalc* C library as a clean Python
API. It further provides some additional tools, not offered by the C library.

For the underlying physics, the finite-element method, and a full reference
manual of all quantities, see `theory/theory.pdf` in the *AlkCalc* project
directory.


## Prerequisites

> **The Python package** ***pyalkcalc*** **requires** ***AlkCalc*** **to be
>   installed and set up on your system.**
> Without it, *pyalkcalc* will not work. In particular:
>
> 1. *AlkCalc* must be compiled and the shared library must be available
>    system-wide (or on your `LD_LIBRARY_PATH`).
> 2. The eigenenergy and radial eigenstate data for the atom or ion species of
>    interest must have been generated beforehand using the data generation
>    pipeline of *AlkCalc*. See the *AlkCalc* README and `theory/theory.pdf` for
>    instructions.

*pyalkcalc* itself requires:

- Python >= 3.10
- NumPy
- Cython (build only)


## Installation

```bash
pip install pyalkcalc
```

Or from source:

```bash
git clone https://github.com/simon-euchner/AlkCalc
cd AlkCalc/pyalkcalc
pip install .
```

Plotting functionality requires matplotlib, which is not installed
automatically. If you want to use it:

```bash
pip install matplotlib
```


## License

The software *pyalkcalc* is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version. See [LICENSE](LICENSE) for the full license text.


## Contact

Simon Euchner [euchner.se[at]gmail.com]
