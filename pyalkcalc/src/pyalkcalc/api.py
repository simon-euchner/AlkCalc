"""
Public Python API for AlkCalc.
"""

# ==============================================================================
# Imports
# ==============================================================================

from ._core import (
        _enlsj_core,
        _StateCore,
        _rp_core,
        _CGCore,
        _SpinorUncoupledBasis,
        _SpinorCoupledBasis,
        _tau_core,
        _fitof_core,
)
from dataclasses import dataclass
from numpy import sqrt, ndarray, float64
from numpy.typing import NDArray


# ==============================================================================
# Dataclasses
# ==============================================================================

@dataclass(slots=True, frozen=True)
class State:
    """
    Container for state data.

    This dataclass represents the metadata associated with a computed
    radial eigenstate for a given set of quantum numbers. It provides a
    lightweight, immutable snapshot of the state returned by the underlying C
    implementation.

    Attributes
    ----------
    N : int
        Total number of discretisation points.
    dim : int
        Dimension of the basis used in the calculation.
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number j = |l - 1/2| or |l + 1/2|.
    t : NDArray[float64]
        Discretisation points. This is a numpy array which contains the
        the points tk, k = 0, ..., N-1, which discretise the interval [0, tmax].
        Here, tmax represents the maximally considered radius rmax = tmax * aB
        and N is the number of discretisation points tk. Note that t0 = 0 and
        tN-1 = tmax.
    h : NDAarray[float64]
        Step sizes. This is a numpy array of length N-1, with entries
        hk = tk - tk-1, k = 1, ..., N-1. These are the (non-uniform) step sizes
        of the grid defined by tk.
    fnlsj : NDArray[float64]
        Values of the radial eigenstate. This is a numpy array which contains
        the values fnlsj(tk), k = 1, ..., dim, where dim = N-2. Note that the
        boundary terms are by construction zero, i.e.,
        fnlsj(t0) = fnlsj(tN-1) = 0.

    Notes
    -----
    The total wave-function is defined as the product of the angular spinor and
    the radial eigenstate Rnlsj. This dataclass represents Rnlsj in terms of
    fnlsj(t) = sqrt(aB) * r * Rnlsj(r), where t = r / aB with Bohr's radius aB.
    Note that fnlsj is dimensionless, because the dimension of Rnlsj is
    1 / sqrt(aB)**3.
    """
    N: int
    dim: int
    n: int
    l: int
    j: float
    t: NDArray[float64]
    h: NDArray[float64]
    fnlsj: NDArray[float64]

@dataclass
class CG:
    """
    Container for Clebsch-Gordan coefficient.

    This dataclass represents the metadata associated with a Clebsch-Gordan
    coefficient. It provides a lightweight, immutable snapshot of the
    Clebsch-Gordan coefficient returned by the underlying C implementation.

    Attributes
    ----------
    sign : int
        Sign of the Clebsch-Gordan coefficient (see notes below).
    numerator : int
        Numerator of the Clebsch-Gordan coefficient (see notes below).
    denominator : int
        Denominator of the Clebsch-Gordan coefficient (see notes below).

    Notes
    -----
    The Clebsch-Gordan coefficient is returned in three parts: its s, its
    numerator, and its denominator. This allows to obtain it analytically, not
    just as a floating-point number. Explicitly, the Clebsch-Gordan coefficient
    C can be computed in the following way:

        C = s * sqrt(numerator / denominator) .

    The Clebsch-Gordan determines the change of basis between the states
    |j1,m1>|j2,m2> and |j,mj>. We choose the Condon-Shortley phase convention,
    that is, we choose the Clebsch-Gordan coefficient with mj = j1 + j2 to be
    positive.
    """
    sign: int
    numerator: int
    denominator: int

@dataclass
class Spinor:
    """
    Container for spinor data.

    The spinor is represented in the basis (|u>, |d>), where |u> and |d> are the
    spin-1/2 eigenstates with spin +1/2 and -1/2, respectively.

    Attributes
    ----------
    u: float
        Component of spinor associated to +1/2.
    d: float
        Component of spinor associated to +1/2.
    """
    u: float
    d: float


# ==============================================================================
# Public Functions
# ==============================================================================

def energy(species: str, n: int, l: int, j: float) -> float:
    """
    Eigenenergy of fine-structure state.

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

    Returns
    -------
    float
        Eigenenergy corresponding to the quantum numbers n, l, s = 1/2, and j in
        units of Hartree.
    """
    species_ascii = species.encode("ascii")
    return _enlsj_core(species_ascii, n, l, j)

def state(species: str, n: int, l: int, j: float, result: str = "f") -> State:
    """
    Radial eigenstate.

    The total wave-function is defined as the product of the angular spinor and
    the radial eigenstate Rnlsj. This function returns the state
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
    result : str
        Character for deciding if discretisation data is returned. If it is
        `"f"` (full), the full result is returned. If result is `"p"` (partial),
        the discretisation data is not returned. In the latter case, the arrays
        `t` and `h` are empty (see documentation of the dataclass: State).

    Returns
    -------
    State
        Radial eigenstate.
    """
    species_ascii = species.encode("ascii")
    result_ascii = result.encode("ascii")
    core = _StateCore(species_ascii, n, l, j, result_ascii)
    core.t.flags.writeable = False
    return State(
            N=core.N,
            dim=core.dim,
            n=core.n,
            l=core.l,
            j=core.j,
            t=core.t,
            h=core.h,
            fnlsj=core.fnlsj,
    )

def radial_matrix_element(
        species: str,
        nb: int,
        lb: int,
        jb: float,
        p: float,
        nk: int,
        lk: int,
        jk: float
) -> float:
    """
    Radial matrix element.

    Parameters
    ----------
    nb : int
        Principal quantum number of bra.
    lb : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1 of bra.
    jb : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2| of bra.
    p : float
        Power of radius operator. For instance, if `p` is zero, the overlap
        between the radial eigenstates is calculated, or if `p` is unity, the
        radial dipole matrix element.
    nk : int
        Principal quantum number of ket.
    lk : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1 of ket.
    jk : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2| of ket.

    Returns
    -------
    float
        Radial matrix element <R', r**p R> of power `p`. Here, R' and R are the
        radial eigenstates (bra and ket) with the associated quantum numbers
        (nb, lb, jb) and (nk, lk, jk).
    """
    species_ascii = species.encode("ascii")
    return _rp_core(species_ascii, nb, lb, jb, p, nk, lk, jk)

def clebsch_gordan_coefficient(
        j1: float,
        m1: float,
        j2: float,
        m2: float,
        j: float,
        mj: float,
        result: str = "analytic"
) -> CG:
    """
    Clebsch-Gordan coefficient.

    Parameters
    ----------
    float : j1
        Angular momentum quantum number of first factor.
    float : m1
        Magnetic quantum number of first factor.
    float : j2
        Angular momentum quantum number of second factor.
    float : m2
        Magnetic quantum number of second factor.
    float : j
        Total angular momentum quantum number.
    float : mj
        Total magnetic quantum number.
    result : str
        If `result` is `"analytic"`, an instance of the dataclass `CG` is
        returned and the Clebsch-Gordan coefficient is returned analytically. If
        `result` is `"numeric"` the Clebsch-Gordan coefficient is returned as
        a real number.

    Returns
    -------
    CG, float
        Clebsch-Gordan coefficient for coupling j1 and j2 to obtain j.
    """
    core = _CGCore(j1, m1, j2, m2, j, mj)
    cg = CG(
          sign=core.sign,
          numerator=core.numerator,
          denominator=core.denominator,
    )
    if result is None or result == "analytic":
        return cg
    elif result == "numeric":
        return cg.sign * sqrt(cg.numerator / cg.denominator)
    else:
        raise ValueError("Invalid value for `result`. Expected `\"analytic\"` "
                         "or `\"numeric\"`.")

def spinor_uncoupled_basis(
        l: int,
        ml: int,
        ms: float,
        theta: float,
        phi: float
) -> Spinor:
    """
    Spinor associated to an uncoupled basis state.

    Parameters
    ----------
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    ml : int
        Magnetic quantum number of obital angular momentum.
    ms : float
        Magnetic quantum number of electron spin.
    theta : float
        Polar angle (Zenitwinkel) in the range [0, pi].
    phi : float
        Azimuthal angle (Azimut) in the range [0, 2pi].

    Returns
    -------
    Spinor
        Spinor associated to an uncoupled basis state. It is computed as
        follows:

            [ Ylml, 0 ] if `ms` is +1/2 (up) ,

            [ 0, Ylml ] if `ms` is -1/2 (down) .

        Here, Ylml is a spherical harmonic. For these, we assume the
        Condon-Shortley phase convention.
    """
    core = _SpinorUncoupledBasis(l, ml, ms, theta, phi)
    return Spinor(
            u=core.u,
            d=core.d,
    )

def spinor_coupled_basis(
        l: int,
        j: float,
        mj: float,
        theta: float,
        phi: float
) -> Spinor:
    """
    Spinor associated to a coupled basis (fine-structure) state.

    Parameters
    ----------
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number.
    mj : float
        Total magnetic quantum number.
    theta : float
        Polar angle (Zenitwinkel) in the range [0, pi].
    phi : float
        Azimuthal angle (Azimut) in the range [0, 2pi].

    Returns
    -------
    Spinor
        Spinor associated to a coupled basis (fine-structure) state. It is
        computed as follows:

            [ C(l,mj-1/2,1/2,+1/2,j,mj) * Ylmj-1/2,
              C(l,mj+1/2,1/2,-1/2,j,mj) * Ylmj+1/2 ] .

        Here, Ylml is a spherical harmonic and C(j1,m1,j2,m2,j,mj) is the
        Clebsch-Gordan coefficient which couples j1 and j2 to yield j. For both,
        we assume the Condon-Shortley phase convention.
    """
    core = _SpinorCoupledBasis(l, j, mj, theta, phi)
    return Spinor(
            u=core.u,
            d=core.d,
    )

def oscillator_strength(
        species: str,
        ni: int,
        li: int,
        ji: float,
        nf: int,
        lf: int,
        jf: float
) -> float:
    """
    Oscillator strength between fine-structure states.

    Parameters
    ----------
    species : str
        String to specify atom/ion species, e.g., 1H for Hydrogen, or 88SR+ for
        the 88Sr+ ion.
    ni : int
        Principal quantum number of initial state.
    li : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1 of initial
        state.
    ji : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2| of final
        state.
    nf : int
        Principal quantum number of initial state.
    lf : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1 of final
        state.
    jf : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2| of
        final state.

    Returns
    -------
    float
        Oscillator strength between the fine-structure states associated to the
        initial and final state with quantum numbers (ni, li, ji) and
        (nf, lf, jf), respectively. Note that this result is obtained after
        summation over mjf. Because of this, the oscillator strength is the same
        for all mji. Therefore, mji is not an argument here.
    """
    species_ascii = species.encode("ascii")
    return _fitof_core(species_ascii, ni, li, ji, nf, lf, jf)

def lifetime(
    T: float,
    species: str,
    n: int,
    dn: int,
    l: int,
    j: float
) -> float:
    """
    Lifetime of fine-structure state.

    Parameters
    ----------
    T : float
        Temperature of black-body excitation spectrum in Kelvin.
    species : str
        String to specify atom/ion species, e.g., 1H for Hydrogen, or 88SR+ for
        the 88Sr+ ion.
    n : int
        Principal quantum number.
    dn : int
        At finite temperature, decay also takes place to states with higher
        energy. The number `dn` specifies how many states, with energy above
        the requested one, are are taken into account for calculating the
        lifetime. For instance, if `dn` is 5, decay up to (and including) states
        with n+5 is considered.
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2|.

    Returns
    -------
    float
        Lifetime of the fine-structure state with quantum numbers n, l, s = 1/2,
        and j in nanoseconds.
    """
    species_ascii = species.encode("ascii")
    return _tau_core(T, species_ascii, n, dn, l, j)
