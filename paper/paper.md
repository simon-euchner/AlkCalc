---
title: "AlkCalc: Alkali-metal atom and alkaline-earth-metal ion calculator"
tags:
  - atomic physics
  - Rydberg ions
  - Rydberg atoms
  - C
  - FORTRAN
  - Python
authors:
  - given-names: Simon
    surname: Euchner
    orcid: 0009-0007-5445-359X
    affiliation: "1"
  - given-names: Weibin
    surname: Li
    orcid: 0000-0001-6731-1311
    affiliation: "2"
  - given-names: Igor
    surname: Lesanovsky
    orcid: 0000-0001-9660-9467
    affiliation: "1,2"
affiliations:
  - index: 1
    name: "Institut für Theoretische Physik, Universität Tübingen, Auf der
           Morgenstelle 14, 72076 Tübingen, Deutschland"
  - index: 2
    name: "School of Physics and Astronomy and Centre for the Mathematics and
           Theoretical Physics of Quantum Non-Equilibrium Systems, The
           University of Nottingham, Nottingham, NG7 2RD, United Kingdom"
date: 10 October 2026
bibliography: paper.bib
---

# Summary

AlkCalc [@AlkCalc] is software designed to numerically compute properties of
alkali-metal atoms and alkaline-earth-metal ions. More specifically, AlkCalc is
designed to compute bound single-electron eigenenergies and eigenstates of the
valence electron together with derived quantities such as transition matrix
elements, oscillator strengths, and lifetimes. The eigenenergies and eigenstates
are computed from the effective, single-electron model potentials given in
[@Marinescu1994;@Aymar1996]. The eigenstates and the derived quantities are of
fundamental relevance in atomic, molecular, and optical (AMO) physics, but are
also central to the design and characterisation of new technologies based on
trapped atoms and ions. Examples include protocols for quantum information
processing [@Sorensen1999;@Sorensen2000;@Saffman2010;@Zhang2020] and quantum
simulations based on the enhanced properties of highly excited electronic
states, so-called Rydberg states [@Weimer2010;@Mueller2008].

AlkCalc is lightweight and portable software. The library itself is written in C
and conforms to the C99 standard [@C99]. The provided Makefiles are written for
GNU make, but can be easily adapted to any other build system. AlkCalc further
builds on decades-long tested FORTRAN routines adapted from Netlib
[@Dongarra1987]. However, this third-party software is built into AlkCalc
directly and therefore does not constitute a dependency. Because of this design
choice, it is straightforward to port AlkCalc to virtually any computing
platform. Its core’s only dependencies are a C library, a C compiler, and a
FORTRAN compiler.

AlkCalc is thoroughly documented in the `README`s and `theory/theory.pdf`, which
constitute its user’s guide. The `README`s contain technical information on how
to install and use AlkCalc. On the other hand, `theory/theory.pdf` provides
comprehensive documentation of the library as well as a detailed introduction to
the problems AlkCalc is designed to address and solve. It allows users who are
new to the atomic and Rydberg physics communities to find all relevant
information for the use of AlkCalc in a single, self-contained document. In AMO
physics there can be different conventions for defining common quantities. One
example is the Condon-Shortley phase for the spherical harmonics
[@Friedrich2017]. AlkCalc explicitly documents all conventions used, eliminating
possible ambiguities.

# State of the field

For a long time, the Rydberg atom and Rydberg ion communities relied on small,
independently developed, unpublished codes for numerically computing properties
of low-lying electronic and Rydberg states. For Rydberg atoms, in particular
the alkali metals, two community-driven open-source software packages unified
these private codes: ARC [@Sibalic2017] and PairInteraction [@Weber2017]. ARC is
designed to compute properties of Rydberg states, while PairInteraction focuses
on the computation of interaction potentials between two Rydberg atoms. Both
of these software packages extract the valence electron’s eigenenergies from
experimentally measured quantum defects that must be externally supplied. The
corresponding radial wave functions are subsequently obtained by numerically
integrating the radial Schrödinger equation at these eigenenergies using an
effective, single-electron model potential to describe the valence electron,
introduced in [@Marinescu1994]. Alternatively, PairInteraction offers analytical
Whittaker functions to describe the eigenstates. A recent development in the
numerical methods of the field is the treatment of divalent atoms, i.e., atoms
with two valence electrons: with version 3.0, ARC was extended to the divalent
atoms Ca, Sr, and Yb [@Robertson2021]. Recently, PairInteraction was likewise
extended to divalent atoms using multichannel quantum defect theory, with states
and matrix elements computed by the separate Python package RydState
[@Moegerle2026].

In contrast to alkali-metal atoms, there is no unified codebase for the
electronic states of the valence electron in singly charged alkaline-earth-metal
ions. Moreover, the available software packages ARC and PairInteraction are
tailored to directly output derived electronic quantities, such as transition
matrix elements, oscillator strengths, and lifetimes. These software packages
are not designed to output and handle large subsets of eigenstates and their
associated eigenenergies. In particular, these software packages do not include
functionality to compute the eigenenergies directly from the model potential:
the eigenenergies are externally provided by experimentally measured quantum
defects. Up to this point there is no open-source software available that
obtains *both* eigenenergies *and* eigenstates for alkali-metal atoms and
alkaline-earth-metal ions directly from a single effective model potential.
AlkCalc computes the eigenenergies and eigenstates directly from the model
potentials given in [@Marinescu1994;@Aymar1996]. Because this approach is
different from ARC, PairInteraction, and RydState, which are all based on
quantum defect theory, AlkCalc is implemented as a stand-alone open-source
software project and not as an extension to the existing libraries.

# Statement of need

The target audiences of the software AlkCalc are the atomic and Rydberg physics
communities, including Rydberg-atom/-ion based quantum computing and simulation,
metrology with trapped ions, sensing with Rydberg atoms, and atomic physics with
energetically low-lying states, relevant for atomic clock experiments and
optical lattices. AlkCalc addresses three key requirements of these communities:

(i) **Unified framework for both atoms and ions.** AlkCalc offers a unified
codebase for the systematic treatment of Rydberg states in both neutral
alkali-metal atoms and singly charged alkaline-earth-metal ions. To date, no
open-source software package appears to be available for the latter.

(ii) **Efficient computation of eigenenergies and eigenstates.** AlkCalc
computes *both* eigenenergies *and* eigenstates directly from the model
potentials for alkali-metal atoms and alkaline-earth-metal ions given in
[@Marinescu1994;@Aymar1996]. This functionality has two benefits: first, in
general, effective model potentials depend on parameters that must be
numerically optimised such that the model potential reproduces the measured
spectroscopic energies correctly. For fitting such parameters, it is therefore
crucial to have a tool that allows one to efficiently compute the eigenenergies
associated with a model potential. Secondly, for the theory side of the
target communities, it is important to have software that allows one to
efficiently compute eigenenergies and eigenstates. This is because in quantum
physics this means the problem is solved: from the eigenenergies and eigenstates
the full dynamics can be inferred and simulated on arbitrary time scales.

(iii) **Transparency and portability.** AlkCalc is self-contained and fully
documented. Every step, from the model Hamiltonian to the stored eigenstates, is
thoroughly described in `theory/theory.pdf`, including all applied conventions
in the formulation of the problems AlkCalc is designed to address. The code
remains portable and usable in the long term, because AlkCalc does not build on
interpreted programming languages. Its full functionality can be operated within
C. This makes AlkCalc both a reliable foundation for future software in the
steadily evolving Rydberg and atomic physics communities and a self-contained
entry point for newcomers to the field.

# Software design

AlkCalc follows a two-stage software architecture, which is visualised in
\autoref{fig:SoftwareArchitecture}.

![**Software architecture.** Visualisation of AlkCalc’s software architecture.
The leftmost box represents the eigensolver (Stage I), which is operated by the
user indicated by the box labelled ‘Input’. The computed eigenenergies and
eigenstates are represented by the box in the centre. The rightmost box
represents the library functions (Stage II). Through the user interface,
represented by the lower right box, the library functions are called.
\label{fig:SoftwareArchitecture}](./figures/software_architecture.png)

**Stage I.** This stage consists of an eigensolver that is tailored to solve
for the eigenenergies and eigenstates of alkali-metal atoms and
alkaline-earth-metal ions, described by a spherically symmetric, single-electron
model potential in combination with the spin-orbit interaction contribution. Due
to symmetry, the eigenstates factorise into an angular and a radial part. The
angular part is solved analytically, as discussed in detail in the documentation
of AlkCalc. The radial part is the solution to a one-dimensional boundary-value
problem. AlkCalc employs B-splines of user-controllable order to approximate the
radial eigenvalue problem with a banded, real-symmetric, generalised matrix
eigenvalue problem. The relevant matrix elements are computed with de Boor’s
algorithm [@Boor2001] in combination with Gaußian quadrature [@Golub1969] for
computing the associated integrals. The generalised matrix eigenvalue problem is
subsequently solved with the LAPACK routine `DSBGVX` [@Anderson1999]. For
evaluating matrix-vector products efficiently, the BLAS routine `DSBMV`
[@Anderson1999] is employed. The resulting eigenenergies and eigenstates are
stored on disk. All data is stored as plain text files, so no additional
software is required to read AlkCalc’s output data, i.e., the eigenenergies and
the eigenstates. It should be emphasised that all FORTRAN routines mentioned in
this section are built into AlkCalc directly.

**Stage II.** This stage consists of the *library functions*. These functions
read and process the stored eigenenergies and eigenstates. For the library
functions, AlkCalc includes an optional Python interface, provided as the Python
package pyalkcalc. Note that pyalkcalc is optional and not required to run
AlkCalc.

**Design choices.** In contrast to shooting methods [@Johnson2007], the B-spline
approach yields all eigenstates within a chosen energy range from a single
diagonalisation, without the need for initial guesses of the eigenenergies. This
provides the full set of eigenstates required for the non-perturbative treatment
of external fields. The compact support of the B-splines renders the matrices
banded, which allows AlkCalc to use the banded LAPACK solver `DSBGVX` and
reduces both memory requirements and computational cost compared to dense
solvers. Since the diagonalisation is the computationally expensive step,
AlkCalc computes the eigenenergies and eigenstates once and stores them on disk,
such that derived quantities can be evaluated repeatedly without solving the
eigenvalue problem again. The choice to store all data as plain text keeps the
output portable and readable without adding software requirements, at the cost
of larger files. However, the data requirements remain moderate: the stored
eigenstates of a whole atom or ion species typically occupy between 10 and 100
MiB and stay below 1 GiB even for large bases (e.g., all states with principal
quantum numbers below 120).

# Research impact statement

The software AlkCalc has already been used in the publications
[@Martins2026;@Euchner2026a;@Euchner2026b]. In [@Martins2026], AlkCalc was used
to compute lifetimes of Rydberg states in ^40^Ca^+^. In [@Euchner2026a], AlkCalc
was used to compute Rydberg-state lifetimes, eigenenergies, and transition
matrix elements (and from these results Rydberg-state polarisabilities) for
^88^Sr^+^.  In [@Euchner2026b], AlkCalc was used to compute Rydberg-state
lifetimes, eigenenergies, and transition matrix elements in the ions ^40^Ca^+^
and ^88^Sr^+^ and the atom ^85^Rb. AlkCalc is actively maintained, and its
development will be continued within further PhD work on Rydberg ions.

# Planned work and future development

A planned extension of AlkCalc is the incorporation of relativistic corrections,
which are essential for the precise treatment of heavy atoms and ions
[@Friedrich2017]. The B-spline method employed by AlkCalc makes this extension
straightforward, as it provides direct access to higher-order derivatives of
numerically represented wave functions, which are required for the relativistic
correction to the kinetic energy.

# AI usage disclosure

For the preparation of the manuscript, Anthropic’s Claude (Opus 5.5) was used
for spell checking, consistency checks, and language and grammar suggestions.
Anthropic’s Claude (Opus 5) helped to search for inconsistencies and bugs in the
code and to write a test suite for AlkCalc. All AI-assisted work was carefully
reviewed by S.E., who takes full responsibility for it. The code for the library
AlkCalc itself was written by S.E. personally.

# Acknowledgements

S.E. is grateful for fruitful discussions with J. W. P. Wilkinson and W. S.
Martins. This work was supported by the European Union through the ERC grant
OPEN-2QS (Grant No. 101164443) and the Deutsche Forschungsgemeinschaft (DFG)
within the Research Unit FOR 5413 (Grant No. 465199066) and the Research Unit
FOR 5522 (Grant No. 499180199). We also acknowledge funding through
JST-DFG 2024: Japanese-German Joint Call for Proposals on ‘Quantum Technologies’
(Japan-JST-DFG-ASPIRE 2024) under JST Grant No. JPMJAP24C2 and DFG Grant No.
554561799.

# References
