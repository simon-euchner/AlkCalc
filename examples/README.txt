/* -------------------------------------------------------------------------- *
 * Example programs and tests for the library functions of AlkCalc            *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */


Notes.

    In each example's source file, the header includes the command that is used
    to compile the example into a binary. Please change the compiler to the
    appropriate one if the example shall be compiled on a system not hosting
    GNU's C library.

    Note that before some of the examples can be run, the corresponding data
    (eigenenergies and radial eigenstates) must be generated. Please see the
    examples' source code and the manual in theory/theory.pdf for details.

    The file tests.c collects all tests of AlkCalc's library functions in a
    single source file, and the Makefile in this directory runs all of them at
    once when it is called with the argument 'test'. The tests are performed
    exclusively for the species 1H. Before they can be run, the data for 1H
    must be generated, and the mass correction must not be included. Please see
    the header of tests.c and the section 'Testing' in README.txt for details.

    Of course, make sure the library AlkCalc is compiled before trying to use
    the example programs.
