/* -------------------------------------------------------------------------- *
 * Python interface for AlkCalc                                               *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */


Description.


    !!! Python bindings for libX (C library for data generation and analysis) !!!


Structure of this README.

        This file is organised as follows: ...


Software requirements.

    ----------------------------------------------------------------
    Software             Example (tested)
    ----------------------------------------------------------------
    Python               v???
    ----------------------------------------------------------------


Installation.

    !!! Requires: libalkcalc.so installed in system path !!!
    !!! must be available to the dynamic linker (ld.so) !!!

    The installation process is split into ???? steps.

    1. Make sure that on the system a BLAS (e.g. Ref. [1]) implementation and a
       version of LAPACK (see Ref. [2]) is installed.

    2. Navigate to the directory LANCZOS and run the Makefile. This builds the
       relevant part of ARPACK (see Ref. [3]). The Makefile is written for GNU
    ...


Structure of this README.

        This file is organised as follows: ...


Software requirements.

    ----------------------------------------------------------------
    Software             Example (tested)
    ----------------------------------------------------------------
    Python               v???
    ----------------------------------------------------------------


Installation.

    The installation process is split into ???? steps.

    1. Make sure that on the system a BLAS (e.g. Ref. [1]) implementation and a
       version of LAPACK (see Ref. [2]) is installed.

    2. Navigate to the directory LANCZOS and run the Makefile. This builds the
       relevant part of ARPACK (see Ref. [3]). The Makefile is written for GNU
    ...

    3. Put a symlink to the lib in the directory where python packages are
       installed
