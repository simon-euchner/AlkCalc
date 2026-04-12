/* -------------------------------------------------------------------------- *
 * PyAlkCalc: Python interface to AlkCalc                                     *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */


Contact.

    With doubts, questions, and/or suggestions, do not hesitate to contact the
    author.

        Simon Euchner
        Electronic mail: <euchner.se@gmail.com>


Structure of this README.

        This README file is organised into the following sections ...


Software requirements.

    ----------------------------------------------------------------
    Software             Example (tested)
    ----------------------------------------------------------------
    ----------------------------------------------------------------
    AlkCalc              ...

    Python               ...
    ----------------------------------------------------------------


Installation.

    The installation process is split into ???????????????????? steps.

    1. Make sure that on the system a BLAS (e.g. Ref. [1]) implementation and a
       version of LAPACK (see Ref. [2]) is installed.

    2. Navigate to the directory LANCZOS.D and run the Makefile. This builds the
       relevant part of ARPACK (see Ref. [3]). The Makefile is written for GNU
       make. If GNU make is not available, just adjust the Makefile accordingly,
       or manually execute the steps defined in the Makefile.
