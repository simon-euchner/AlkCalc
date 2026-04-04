/* -------------------------------------------------------------------------- *
 * AlkCalc: Alkali-metal atom and alkaline-earth-metal ion calculator         *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

(Deutsche Version: "READMEger.txt")


Contact.

    With doubts, questions, and/or suggestions, do not hesitate to contact the
    author.

        Simon Euchner
        Electronic mail: <euchner.se@gmail.com>


Introduction.

        Physics involving trapped atoms and ions is of great interest across the
    natural sciences. Particularly, quantum computing and quantum simulation
    applications based on ultracold atoms and ions have gained popularity over
    the last decades, as the necessary control over these complex synthetic
    quantum systems became a tangible reality through extraordinary experimental
    advancement.

        To construct and develop the systems of interest, it is indispensable to
    be able to compute atomic physics quantities of the various atom and ion
    species. The atoms and ions are approximated as single-electron systems,
    with all electrons residing in their ground state configuration except for
    one, which defines the single-electron problem in combination with an
    effective nucleus, i.e., the nucleus surrounded by the bound electron cloud.
    Since the mid-1990s, see Refs. [6,8], model potentials are employed
    for treating these problems numerically. Moreover, the model potentials and
    experimentally obtained quantum defects are used in software packages such
    as Ref. [10]. Typically, these software packages are rather extensive, which
    makes it challenging and tedious to gain full understanding of their
    internal workings. Additionally, typical packages are developed with a
    specific high-level programming language, like Python, in mind, which forces
    scientists into learning the 'correct' language. We note that this is of
    course exactly what most people need. However, if systems need to be
    described where it is not sufficient to compute some common quantity, such
    as a C6 dispersion coefficient, the full eigenstates are needed. If these
    eigenstates are only avilable through quantum defects, this highly limits
    the breadth of describable systems, as the defects must be known --- a
    different, more consistent, approach is to also compute the eigenenergies
    from the model potential, without relying on quantum defects.

        As a solution to the problems stated above, we present the library
    AlkCalc (Alkali-metal atom and alkaline-earth-metal ion calculator), which
    is developed with a different philosophy in mind: We want to provide simple
    software that can do everything, but includes only the necessary features.
    AlkCalc is a collection of simple programs written in native C and FORTRAN.
    The software consists of two parts: the library functions, which are used to
    obtain single atom or ion eigenenergies, eigenstates, Clebsch-Gordan
    coefficients, and (virtually) arbitrary radial matrix elements, as well as
    the data generation part, allowing to compute the desired eigenenergies and
    radial eigenstates. The generated data is then used by the library functions
    to provide the data in an easy-to-use fashion to the user.
        AlkCalc is very much 'hands-on' in the sense that it is simple enough to
    directly interact with the source code. Further, because it is written in a
    low-level language conforming to the C99 standard, it is virtually platform
    independent. Importantly, since high-level languages typically come with
    easy-to-use protocols to call C code (e.g. ctypes or Cython for Python), it
    is easy to write an interface to the high-level language of choice (for
    instance Python, Julia, Matlab, etc. --- whatever is popular). Therefore, we
    view AlkCalc not just as a ready-to-use software package, but also as the
    basis for development of other software packages written in high-level
    languages.
        We note that because of the simplicity of AlkCalc, it is very easily
    extendable, e.g., adding external fields to the single-atom or single-ion
    Hamiltonian, without the need of relying on perturbative methods, is no
    problem, as the source code, opposed to typical software packages, consists
    only of a few files without 'hiding' details. Finally, the parameters for
    the model potentials are included in a single text file, i.e., are easy to
    find and can be directly accessed. In this context we note that all data is
    stored in simple text format files as well, having the advantage of not
    introducing additional software requirements to read binary data formats and
    alike.

        In conclusion, the software AlkCalc is for researchers who value full
    control over how data for their research is generated and who want
    transparent software: nothing is 'hidden', all steps are laid out, no
    unnecessarily complicated data formats are used, and the software is simple
    and lightweight while maintaining completeness; that is, the eigenenergies
    and eigenstates are available, i.e., the problem is solved! It is fully open
    source and the user is actually encouraged to look at and interact with the
    source code.


Structure of this README.

        This README file is organised into the following sections. In the
    section 'Software requirements' we state the (very few) required external
    software (most of the employed external software is already included in the
    form of source code [text] files in AlkCalc itself --- only the necessary
    part to avoid unneeded files and reduce complexity). In the section
    'Installation' we provide instructions on how AlkCalc can be installed and
    set up. In the section 'Data generation' it is explained how the
    eigenenergies and radial eigenstates can be computed and stored. Finally, in
    the section 'Important additional information' we discuss (technical)
    aspects which should be considered before employing AlkCalc.
        In addition to this README, we provide the file 'theory.d/theory.pdf',
    which contains an introduction to the physics of the single-atom and
    single-ion Hamiltonian which AlkCalc is designed to diagonalise, as well as
    a full description of the finite-element method we employ to reduce the
    problem to a matrix eigenproblem. Finally, 'theory.d/theory.pdf' includes an
    exhaustive reference manual for the library functions of AlkCalc, including
    all the information the user may need.


Software requirements.

    ----------------------------------------------------------------
    Software             Example (tested)
    ----------------------------------------------------------------
    ----------------------------------------------------------------
    BLAS                 blas-3.12.1-2 [1]

    LAPACK               lapack-3.12.1-2 [2]

    C Compiler           gcc-15.2.1+r22+gc4e96a094636-1 [12]

    FORTRAN Compiler     gcc-fortran-15.2.1+r22+gc4e96a094636-1 [12]

    C Library            glibc-2.42+r17+gd7274d718e6f-1 [13]
    ----------------------------------------------------------------

    - AlkCalc conforms to the C99 standard, Ref. [11], with the additional
      assumption that the fixed-width integer types, int8_t, int32_t, and
      int64_t are defined. The C99 standard defines these types as optional (see
      Sec. 7.18.1.1 in Ref. [11]). However, most modern C libraries (e.g. GNU's
      C library) define these types, and the C99 standard ensures that, if these
      types are defined, they are two's-complement represented and come without
      padding (see Sec. 7.18.1.1 in Ref. [11]). We make use of these types and
      their properties for definiteness and to simplify overflow checks in
      integer arithmetic. We remark that if these types are not defined on the
      system (with a C library conforming to C99 standard) the compiler throws
      an error.


Installation.

    The installation process is split into seven steps.

    1. Make sure that on the system a BLAS (e.g. Ref. [1]) implementation and a
       version of LAPACK (see Ref. [2]) is installed.

    2. Navigate to the directory LANCZOS.D and run the Makefile. This builds the
       relevant part of ARPACK (see Ref. [3]). The Makefile is written for GNU
       make. If GNU make is not available, just adjust the Makefile accordingly,
       or manually execute the steps defined in the Makefile.

    3. Navigate to the directory LUFac.d and run the Makefile (regarding the
       Makefile, the same conditions apply as in step 2). This builds the
       relevant part of the LU-factorisation software SuperLU (see Ref. [7]).

    4. Navigate to the directory interface.d and set the variable
       PATH_TO_STATES, which defines the location to which the data for the
       radial eigenstates will be written.

--- Note (*)

    5. Navigate to the directory interface.d and set the relevant paths in the
       file alkcalk.h. The variable PATH_TO_ALKCALC must be the absolute path to
       the place where the directory AlkCalc is located. The variable
       PATH_TO_STATES defines where the radial eigenstates will be read from.
       The radial eigenstates take a lot of disk space (~10 GB), which allows
       for more flexibility; for example, one could store the radial eigenstates
       on an external drive.

    6. Run the Makefile in the directory AlkCalc with the argument 'lib'. This
       builds the library functions of AlkCalc, which are the ones meant for
       end-user interaction. For a full manual of the library functions see
       'theory.d/theory.pdf'.

    7. To use the library functions, the associated data, i.e., eigenenergies
       and radial eigenstates, must be computed first. This process is described
       in the next section.

(*) This marks the point at which the eigenenergies and radial eigenstates can
    be computed. The following steps are only required if the library functions
    of AlkCalc (see 'theory.d/theory.pdf' for full documentation of the library
    functions) are desired.


Data generation.

        In this section we describe how AlkCalc can be used to (partially)
    diagonalise the full single-atom or single-ion Hamiltonian presented in the
    theory part in 'theory.d/theory.pdf'. This part of AlkCalc is crucial for
    generating the data (eigenenergies and radial eigenstates) which the library
    functions (see Sec. Installation) depend on. In the following we focus on
    the atom or ion species X. The four steps below describe how the data for X
    can be generated.

    1. Navigate to the directory interface.d and open the file species.dat. Make
       sure that in this file the necessary data for the species X is located.
       If a new species is added, make sure to keep the formatting correct (cf.
       already available entries).

    2. Open the file settings.c and set the parameters. The species identifier
       which refers to the species X is defined in species.dat. There are two
       important parameters: offset and shift. These should be chosen in the
       following way: Suppose the ground state energy (ionisation energy) of X
       is -0.5 Hartree. First set 'offset' such that the energies are bunched
       around 1.0, i.e., in our example something like 1.3. The parameter
       'shift' can then be set to 1.0, which means that the algorithm searches
       for eigenenergies around 1.0. Generally, it is most secure to set
       'offset' such that all eigenenergies are larger than 'shift', e.g., 1.6
       in our example. The resulting eigenenergies are automatically corrected
       to the actual eigenenergies in units of Hartree (with the negative
       ionisation energy being the ground state energy). Note that these
       parameters have a significant impact on runtime, so testing different
       choices is encouraged. To test which ones yield the best performance,
       repeat step 3 for different choices.

    3. Run the Makefile in the top-level directory AlkCalc with the argument
       'solve'. This will generate the eigenenergies and radial eigenstates. The
       radial eigenstates are stored in a user-specified location (see
       Sec. Installation) and the eigenenergies are stored in data.d. Check if
       the ground state energy is correctly captured. Otherwise, readjust the
       parameters 'offset' and 'shift'.

       IMPORTANT: Depending on the maximum desired principal quantum number,
                  choose 'rmax' large enough and keep in mind that larger core
                  charges result in radial eigenstates with support at smaller
                  distances, that is, typically, 'rmax' must be larger for atoms
                  compared to ions.

    4. From now on, KEEP THE PARAMETERS species, N, nmax, rmax in
       'interface.d/settings.c' FIXED, and only change the orbital (l) and total
       angular momentum (j). For each desired pair (l, j) generate the
       eigenenergies and radial eigenstates by running the Makefile with the
       argument 'solve'. Importantly note that it might be necessary to readjust
       the parameters 'offset' and 'shift' while generating the data for
       different pairs (l, j). For instance, it could be that numerically an
       eigenenergy is obtained that is lower than the ground state energy;
       readjusting the parameters resolves this problem and yields the correct
       lowest eigenenergy (larger than the ground state energy) of the series
       (l, j).

       IMPORTANT: The discretisation data (positions and step sizes) is stored
                  in data.d for each species exactly ONCE. Therefore, after
                  testing different settings for the associated parameters,
                  delete this file manually. It will then be regenerated ONCE
                  for the final parameters that the user selected. This might
                  sound somewhat convoluted but can be summarized as follows:

                  - Check if eigenenergies are computed efficiently by testing
                    different parameters, as described above (point 3).

                  - Once finished, keep the parameters fixed for all pairs
                    (l, j), and before generating the first set of eigenenergies
                    and radial eigenstates, delete the file containing the
                    discretisation data in 'data.d'. When the solver is run for
                    the next pair (l, j), the file containing the discretisation
                    data is NOT overwritten --- it is only regenerated in case
                    the data file could not be located, i.e., was deleted.


Important additional information.

    - By default, the mass correction factor C from 'theory.d/theory.pdf' is set
      to unity in the source code in 'src.d/eigensolver.c', that is, we make the
      approximation that the reduced mass is equal to the electron's mass. This
      is due to the fact that the model parameters for the model potential we
      use by default (see Refs. [6,8]) have been computed without the mass
      correction. This we conclude from the fact that the computed ground state
      energies fit the ideal values for the ionisation energy better if we omit
      the mass correction. If model parameters obtained WITH the mass correction
      are employed, one must include the mass correction in the source code as
      well. This amounts to uncommenting one line (and commenting out another
      one) in the file 'src.d/potential.c' (please see also the explanation in
      this file).

    - Alkcalc includes the Hydrogen atom (1H) and the Helium ion (4HE+). Both
      including Russell-Saunders (LS) coupling and, by default, NO mass
      correction, i.e., we assume that the reduced mass is the electron's mass.
      We provide these species to benchmark results for both atoms and ions.
      Note that LS-coupling is weak (especially for small orbital angular
      momentum quantum numbers and small n), such that to test if the settings
      for the eigensolver in 'interface.d/settings.c' are justified, it is
      sufficient to compare the eigenenergies to the analytical result for 1H
      and 4HE+ without LS-coupling.

    - By default, AlkCalc uses step sizes that increase linearly over the
      interval [0, rmax]. However, this can be adjusted by altering the function
      'step' in 'src.d/eigensolver.c'. More information can be found directly in
      the corresponding source code.

    - All ion masses in the file 'interface.d/species.dat' are the full atom's
      mass MINUS the electron's mass.


References.

    Note that also the references from the source code are listed here.

    [1] 'OpenBLAS: An optimised BLAS Library',
        URL: http://www.openmathlib.org/OpenBLAS/
    [2] 'LAPACK---Linear Algebra PACKage', URL: https://netlib.org/lapack/
    [3] 'Opencollab, ARPACK-NG', URL: https://github.com/opencollab/arpack-ng
    [4] J. W. P. Wilkinson, K. Bolsmann, T. L. M. Guedes, M. M\"uller, and
        I. Lesanovsky, New J. Phys. 27, 064502 (2025)
    [5] 'NIST: National Institute of Standards and Technology (NIST)',
        URL: https://www.nist.gov
    [6] M. Aymar, C. H. Greene, E. Luc-Koenig, Rev. Mod. Phys. 68, 1015 (1996)
    [7] 'SuperLU 7.0.1',
        URL: https://portal.nersc.gov/project/sparse/superlu/ug.pdf
    [8] M. Marinescu, H. R. Sadeghpour, and A. Dalgarno, Phys. Rev. A 49, 982
        (1994)
    [9] 'Commission on isotopic abundances and atomic weights (CIAAW)',
        URL: https://www.ciaaw.org/lithium.htm
   [10] N. Šibalić, J. D. Pritchard, C. S. Adams, and K. J. Weatherill, Comput.
        Phys. Commun. 220, 319–331 (2017)
   [11] INTERNATIONAL STANDARD ISO/IEC 9899:1999(E) (American National Standard
        Institute, New York, 1999) 2nd ed.
   [12] 'GCC, The GNU Compiler Collection', URL: https://gcc.gnu.org
   [13] 'The GNU C Library', URL: https://www.gnu.org/software/libc
