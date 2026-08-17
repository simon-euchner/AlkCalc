/* -------------------------------------------------------------------------- *
 * AlkCalc: Alkali-metal atom and alkaline-earth-metal ion calculator         *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

(Deutsche Version: "READMEger.txt")


Contact.

    If you have doubts, questions, and/or suggestions, please do not hesitate to
    contact the author.

        Simon Euchner
        Electronic mail: <euchner.se@gmail.com>


Introduction.

        Bound single-electron eigenfunctions of alkali-metal atoms and
    alkaline-earth-metal ions are of significant relevance in quantum computing,
    quantum information processing, quantum simulation, Rydberg physics, quantum
    optics, and more. However, except for the special case of the Hydrogen atom
    and Hydrogenic ions, the time-independent Schrödinger equation cannot be
    solved analytically. For this reason, eigenenergies and eigenstates must be
    obtained numerically. The library AlkCalc is designed for exactly this
    purpose: it allows the user to compute the relevant part of the bound-state
    energy spectrum together with the associated eigenstates.

        AlkCalc itself is split into two parts. The first part computes the
    eigenenergies and eigenstates once and stores the resulting data. The second
    part consists of the library functions, which implement useful functionality
    for reading and processing the precomputed atomic structure data. For
    instance, with these functions one can compute radial transition matrix
    elements and oscillator strengths.

        AlkCalc is very much 'hands-on', in the sense that it is simple enough
    to interact directly with the source code. Further, because it is written in
    a low-level language conforming strictly to the C99 standard, it is
    virtually platform-independent. Since high-level languages typically offer
    easy-to-use protocols for calling C code (e.g. ctypes or Cython for Python),
    it is straightforward to write an interface to the high-level language of
    one's choice (for instance Python, Julia, MATLAB, etc. --- whatever is
    popular these days). AlkCalc should therefore not solely be seen as a
    ready-to-use software package, but more broadly as a basis for developing
    other software packages written in high-level languages. An interface for
    Python is already provided by AlkCalc (see AlkCalc/pyalkcalc).

        To model the atoms and ions, AlkCalc employs so-called parametric model
    potentials, in particular those introduced in 1994 by Marinescu et al. in
    Ref. [Mar1994] and in 1996 by Aymar et al. in Ref. [Aym1996]. These
    parametric model potentials are described in detail in theory/theory.pdf.
    Unlike typical packages for treating Rydberg atoms, such as ARC [Sib2017] or
    PairInteraction [Web2017], AlkCalc does not rely on quantum defect theory.
    Instead, AlkCalc computes everything self-consistently from the parametric
    model potential. This has the major advantage that AlkCalc is not at all
    restricted to high principal quantum numbers, where quantum defect theory
    becomes accurate, but automatically works equally well for low-lying
    electronic states. Further, not relying on quantum defect theory grants
    AlkCalc the power to solve essentially any problem with radial symmetry.
    Besides other things, this flexibility allows AlkCalc to natively treat both
    Rydberg states in neutral atoms and AND atomic ions --- note that ARC and
    PairInteraction are designed for neutral Rydberg atoms. To keep the
    interface clean and simple, all parameters associated with the parametric
    model potentials are collected in a single text file. In fact, all data
    associated with AlkCalc is stored in plain text files, which has the
    advantage of not introducing additional software requirements for reading
    binary data formats and alike.

        Ultimately, AlkCalc is designed for researchers who value full control
    over the data underlying their research. More broadly, AlkCalc is for anyone
    who values transparent software: nothing is 'hidden', everything is laid out
    plainly, no unnecessarily complicated data formats are used (plain text
    files only), and the software remains simple and lightweight while still
    being complete --- complete in the sense that the eigenenergies and
    eigenstates can indeed be computed, so the problem is actually solved.
    AlkCalc is fully open source, and users are explicitly encouraged to inspect
    and interact with the source code. Finally, AlkCalc comes with no external
    dependencies beyond a C compiler, a FORTRAN compiler, and a C library. This
    is because all other required software is built directly into AlkCalc
    itself, using decades-tested code from Netlib [Don1987].

        In addition to this README, the file theory/theory.pdf provides a
    thorough introduction to the physics of the single-atom and single-ion
    Hamiltonians that AlkCalc is designed to diagonalise. It also contains a
    full description of the numerical (B-spline) method AlkCalc uses to reduce
    the radial eigenvalue problem to a generalised matrix eigenvalue problem.
    Finally, theory/theory.pdf includes a comprehensive reference manual for
    AlkCalc's library functions, covering all the information a user may need.


Structure of this README.

        This README is organised into the following sections. The section
    'Software requirements' lists the required external software. The section
    'Installation' provides instructions for correctly installing and setting up
    AlkCalc. The section 'Data generation' explains how the eigenenergies and
    radial eigenstates are computed and stored with AlkCalc. Finally, the
    section 'Important additional information' discusses technical aspects that
    should be taken into consideration before using AlkCalc.


Software requirements.

    -----------------------------------------------------------------------
    Software            Example (tested)                          Reference
    -----------------------------------------------------------------------
    -----------------------------------------------------------------------
    C Library           glibc 2.44+r24+g16be1518495f-1            [GLC]

    C Compiler          gcc 15.2.1+r22+gc4e96a094636-1            [GCC]

    FORTRAN Compiler    gcc-fortran 15.2.1+r22+gc4e96a094636-1    [GCC]
    -----------------------------------------------------------------------

    - AlkCalc conforms to the C99 standard, Ref. [C99], with the additional
      requirement that the fixed-width integer types int8_t, int32_t, and
      int64_t be defined. The C99 standard treats these types as optional (see
      Sec. 7.18.1.1 in Ref. [C99]). However, most modern C libraries (e.g. GNU's
      C library) define them, and the C99 standard guarantees that, where these
      types are defined, they use two's-complement representation without
      padding bits (see Sec. 7.18.1.1 in Ref. [C99]). AlkCalc relies on these
      types and their properties for definiteness and to simplify overflow
      checks in integer arithmetic. If the fixed-width integer types are not
      defined on your system the compiler will raise an error.


Installation.

    The installation process is split into eight steps.

    1. Navigate to the directory GAUSSQ and run the Makefile. This builds the
       program GAUSSQ [Gol1969] (source code taken from Netlib [Don1987]), which
       is used to compute Gaussian quadratures. The Makefile is written for GNU
       make. If GNU make is not available on your system, simply adjust the
       Makefile accordingly, or manually carry out the steps it performs.

    2. Navigate to the directory BSPLINES and run the Makefile (the same
       conditions as in step 1 apply to the Makefile). This builds a version of
       de Boor's algorithm [dBo2001] written by D. E. Amos [Amo1993].

    3. Navigate to the directory EIGLAPACK and run the Makefile (again, the same
       conditions as in step 1 apply). This builds the two eigensolvers DSBGVX
       and DSBEVX, which are part of LAPACK [And1999] (source code taken from
       Netlib [Don1987]).

    4. Navigate to the directory AlkCalc/interface and set the variable
       PATH_TO_STATES, which defines the location to which the data for the
       radial eigenstates will be written. Typically, it is fine to store this
       data directly in AlkCalc's data directory.

--- Note (*)

    5. Navigate to the directory MVMBLAS and run the Makefile (the same
       conditions as in step 1 apply). This builds the routine DSBMV, a BLAS
       level-2 routine. This routine is part of LAPACK [And1999,Don1987] and
       allows matrix-vector products to be computed efficiently.

    6. Navigate to the directory AlkCalc/interface and set the relevant paths in
       the file alkcalk.h. The variable PATH_TO_ALKCALC must be the absolute
       path to the location of the AlkCalc directory. The variable
       PATH_TO_STATES defines where the radial eigenstates will be read from.
       The radial eigenstates can take up a fair amount of disk space (~5 GB),
       which may be an issue on legacy hardware. Being free to choose where the
       radial eigenstates are stored gives more flexibility here; for example,
       they could be stored on an external drive. On modern hardware, of course,
       disk space is usually not a concern.

    7. Run the Makefile in the top-level directory of AlkCalc with the argument
       'lib'. This builds AlkCalc's library functions, which are the ones meant
       for end-user interaction. For a comprehensive reference manual of the
       library functions see theory/theory.pdf.

    8. Before the library functions can be used, the associated data (i.e. the
       eigenenergies and the radial eigenstates) must be computed. This process
       is described in the next section.

(*) This marks the point at which the eigenenergies and radial eigenstates can
    already be computed. The steps that follow are only required if AlkCalc's
    library functions are also needed (see theory/theory.pdf for comprehensive
    documentation of the library functions).


Data generation.

        This section describes how AlkCalc is used to diagonalise the full
    single-atom or single-ion Hamiltonian presented in theory/theory.pdf. This
    part of AlkCalc is intended for generating the data (eigenenergies and
    radial eigenstates) that the library functions (see Sec. Installation) rely
    on. In what follows, the focus is on a generic atom or ion species X. The
    four steps below describe how the data for X is generated using AlkCalc.

    1. Navigate to the directory AlkCalc/interface and open the file
       species.dat. Make sure this file contains the necessary data for the
       species X. When adding a new species, make sure to keep the formatting
       consistent with the existing entries. When adding new species, please
       observe the following rules and assumptions:

          (1) Make sure entries are sorted in ascending order in the orbital
              angular momentum quantum number, e.g., data for P states must be
              placed before data for F states.

          (2) Data for some orbital angular momentum quantum numbers may be
              omitted, e.g., it is valid to provide data only for P, D, and G
              states. This is useful when, for instance, only P states or
              specific circular states are of interest.

          (3) There must be data for at least one orbital angular momentum
              quantum number in species.dat.

          (4) Suppose an orbital angular momentum quantum number l is requested
              in settings.c that is larger than the largest one, l0, specified
              in species.dat. For l > l0, AlkCalc internally uses l, but with
              the data associated with l0 in species.dat. Note that this
              behaviour is consistent with Refs. [Mar1994,Aym1996].

          (5) In species.dat, a minimum principal quantum number, nl, is
              specified for each l. There are two possibilities: either nl
              follows the Hydrogenic law (i.e., nl = l + 1), or nl is anomalous
              in the sense that nl > l + 1. The correct nl can be read off
              directly from the electronic configuration of the atom or ion
              species. For example, Rubidium has the electronic configuration
              [Kr]5s1, meaning that for S states (l = 0), nl = n0 = 5 > 0 + 1.
              This is the anomalous case and must be specified explicitly in
              species.dat. When data for a requested orbital angular momentum
              quantum number is not supplied explicitly in species.dat, AlkCalc
              internally assumes the Hydrogenic law (i.e., nl = l + 1).

    2. Open the file settings.c and set the parameters. The species identifier
       referring to species X is defined in species.dat. One particularly
       important parameter, offset, must be chosen manually and with some care.
       It controls which eigenenergy is associated with the minimal principal
       quantum number nl. This is necessary because, for anomalous nl, the
       lowest computed eigenenergy is not always the correct one: the potential
       can host eigenenergies lower than the actual ground-state energy, and if
       the generalised eigenvalue problem is badly conditioned, extremely
       negative, unphysical eigenvalues can also appear. To eliminate these
       unwanted eigenvalues, proceed as follows: start with offset = -nl. This
       associates the lowest possible eigenvalue hosted by the potential,
       physical or not, with nl (check data/energies-X-...). If the energy
       associated with nl is not the correct one, increment offset by one
       (offset += 1) and run again 'make solve' (see point 3). Do this until all
       unphysical eigenvalues have been cut off. In practice, the unphysical
       eigenvalues are typically far off and easy to identify.

    3. Run the Makefile in the top-level directory AlkCalc with the argument
       'solve'. This generates the eigenenergies and radial eigenstates. The
       radial eigenstates are stored in the user-specified location (see Sec.
       Installation), and the eigenenergies are stored in AlkCalc's data
       directory. When run with the argument 'solve', AlkCalc also prints the
       condition number of the mass matrix of the generalised eigenvalue problem
       (see theory/theory.pdf). This is only meant as an indicator, not an exact
       mathematical bound, of the potential loss of precision, assuming that the
       Hamiltonian itself is well conditioned. Even so, it is useful for
       choosing the order of the B-splines, k, and the number of knots without
       multiplicities, N. As a rule of thumb, one looses roughly log10(kappa)
       digits, where kappa is the condition number of the mass matrix M (see
       theory/theory.pdf). For 64-bit floating-point arithmetic, one should try
       to keep the condition number kappa less than ~1e6, so roughly 10 digits
       of precision are left. Note that all of this is assuming H is
       well-conditioned. Therefore, this calculation is more a rough estimate
       than an exact estimation of the numerical error. In general, all numbers
       should be chosen within reason, so the results can be trusted.

       IMPORTANT: Depending on the maximum desired principal quantum number,
                  choose rmax large enough, and keep in mind that larger core
                  charges result in radial eigenstates with support at smaller
                  distances --- that is, rmax typically needs to be larger for
                  atoms than for ions.

    4. From this point on, KEEP THE PARAMETERS species, k, N, nmax, and rmax in
       interface/settings.c FIXED, and change only the orbital, l, and total
       angular momentum quantum number, j. For each desired pair (l, j),
       generate the eigenenergies and the radial eigenstates by running the
       Makefile with the argument 'solve'. Note that it may be necessary to
       adjust the offset parameter for each pair (l, j) individually.

       IMPORTANT: The knot data (knot vector and step sizes; see
                  theory/theory.pdf) is stored exactly ONCE per species in
                  AlkCalc's data directory. Therefore, after testing different
                  settings for the associated parameters, delete the file
                  data/knotdata-X.dat manually. It will then be regenerated
                  exactly ONCE for the finally chosen parameters. This may sound
                  somewhat convoluted, but can be summarised as follows:

                  - Check whether the eigenenergies are computed to a reasonable
                    precision by testing different values of the parameters k,
                    N, and rmax.

                  - Once satisfied, keep the parameters k, N, and rmax fixed for
                    all pairs (l, j), and, before generating the first set of
                    eigenenergies and radial eigenstates, delete the file
                    containing the knot data in AlkCalc/data. When the
                    eigensolver is subsequently run for the next pair (l, j),
                    the file containing the knot data is NOT overwritten --- it
                    is only regenerated in case AlkCalc cannot locate the
                    knot-data file, i.e., if it was deleted.


Important additional information.

    - By default, the mass correction factor C from theory/theory.pdf is set to
      unity in the source code in src/eigensolver.c, that is, the reduced mass
      is approximated by the electron mass. This is because the parameters for
      the parametric model potential (see Refs. [Mar1994,Aym1996]) were computed
      without the mass correction: the computed ground-state energies fit the
      ideal ionisation-energy values better when the mass correction is omitted.
      If model parameters obtained WITH the mass correction are used instead,
      the mass correction must also be included in the source code. This amounts
      to uncommenting one line (and commenting out another) in the file
      src/potential.c (see also the explanation given there).

    - AlkCalc includes the Hydrogen atom (1H) and the Helium ion (4HE+), both
      with Russell-Saunders (LS) coupling and, by default, NO mass correction.
      These species are provided to benchmark results for both atoms and ions.
      Note that LS coupling is weak, in particular for small orbital angular
      momentum quantum numbers and large n. Consequently, to check whether the
      eigensolver settings in interface/settings.c are reasonable, it is usually
      sufficient to compare the eigenenergies to the analytical result for 1H
      and 4HE+, even without considering LS coupling.

    - By default, AlkCalc uses for the knots for the B-splines step sizes, i.e.,
      distances between consecutive knots, that increase linearly over the
      interval [0, rmax]. This can be adjusted by modifying the function step in
      src/eigensolver.c. More information can be found directly in the
      corresponding source code.

    - All ion masses given in the file interface/species.dat are the full atom's
      mass MINUS the mass of the removed electron(s).


References.

    Note that references from the remaining files of AlkCalc are also listed
    here.

    [Mar1994] M. Marinescu, H. R. Sadeghpour, and A. Dalgarno, 'Dispersion
              coefficients for alkali-metal dimers', Phys. Rev. A 49, 982 (1994)

    [Aym1996] M. Aymar, C. H. Greene, E. Luc-Koenig, 'Multichannel Rydberg
              spectroscopy of complex atoms', Rev. Mod. Phys. 68, 1015 (1996)

    [Sib2017] N. Šibalić, J. D. Pritchard, C. S. Adams, and K. J. Weatherill,
              'ARC: An open-source library for calculating properties of alkali
              Rydberg atoms', Comput. Phys. Commun. 220, 319–-331 (2017)

    [Web2017] S. Weber, C. Tresp, H. Menke, A. Urvoy, O. Firstenberg, H. P.
              Büchler, and S. Hofferberth, 'Calculation of Rydberg interaction
              potentials', J. Phys. B: At. Mol. Opt. Phys. 50 133001 (2017)

    [Don1987] J. J. Dongarra and E. Grosse, 'Distribution of Mathematical
              Software via Electronic Mail', Commun. ACM 30, 403--407 (1987)

    [GLC]     'The GNU C Library (glibc)', Free Software Foundation, URL:
              https://www.gnu.org/software/libc/

    [GCC]     'GCC, The GNU Compiler Collection', Free Software Foundation, URL:
              https://gcc.gnu.org/software/gcc/

    [C99]     'INTERNATIONAL STANDARD ISO/IEC 9899:1999(E)' (American National
              Standard Institute, New York, 1999) 2nd ed.

    [Gol1969] G. H. Golub and J. H. Welsch, 'Calculation of Gauss Quadrature
              Rules', Math. Comp. 23, 221--230 (1969)

    [dBo2001] C. de Boor, 'A Practical Guide to Splines' (Springer, New York,
              2001) 1st ed., ISBN: 978-0-387-95366-3

    [Amo1993] D. E. Amos, 'Implementation of de Boor's algorithm', URL:
              http://www.netlib.org/slatec/src/dbspvd.f

    [And1999] E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel, J.
              Dongarra, J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney,
              and D. Sorensen, 'LAPACK User's Guide' (Society for Industrial and
              Applied Mathematics, Philadelphia, PA, 1999) 3rd ed., ISBN:
              0-89871-447-8

    References which appear in the code, the other READMEs, and the data files
    of AlkCalc.

    [NISTcuu] 'The NIST Reference on Constants, Units, and Uncertainty', URL:
              https://physics.nist.gov/cuu/Constants/

    [NISTaw]  'Atomic Weights and Isotopic Compositions with Relative Atomic
              Masses', URL: https://www.nist.gov/pml/atomic-weights-and-isotopic
              -compositions-relative-atomic-masses

    [NISTie]  'NIST Atomic Spectra Database Ionization Energies Data', URL:
              https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
