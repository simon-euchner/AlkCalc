"""
Tools for AlkCalc: additional functionality not offered by the C library
"""


# ==============================================================================
# Imports
# ==============================================================================

from .api import (
        state,
)


# ==============================================================================
# Data visualisation
# ==============================================================================

def plot_state(
        species: str,
        n: int,
        l: int,
        j: float,
        ax: "matplotlib.axes.Axes | None" = None
) -> None:
    """
    Plot radial eigenstate.

    The total wave-function is defined as the product of the angular spinor and
    the radial eigenstate Rnlsj. This function plots the state
    fnlsj(t) = sqrt(aB) * r * Rnlsj(r) , where t = r / aB with Bohr's radius aB.
    Note that fnlsj is dimensionless, because the dimension of Rnlsj is
    1/aB**(3/2).

    Parameters
    ----------
    species : str
        String to specify atom/ion species, e.g., 1H for Hydrogen, or 88SR+ for
        the 88Sr+ ion.
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2|.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If `None`, a new figure is created.
    """

    ### Import pyplot or throw an error
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError(
                "Matplotlib is required for plotting, but it is not installed."
        )

    ### Get state
    s = state(species, n, l, j)

    ### Create figure
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(s.t, [0, *s.fnlsj, 0], lw=1)
    ax.set_xlabel(r"$t=r\slash{a_\mathrm{B}}$")
    ax.set_ylabel(r"$\sqrt{a_\mathrm{B}}\times{r}\times{R}_{n,l,s,j}(r)$")
