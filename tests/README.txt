/* -------------------------------------------------------------------------- *
 * Testing of AlkCalc                                                         *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */


Notes.

    The file tests.c contains the code for all tests of AlkCalc and the Makefile
    in this directory runs all of the tests at once, when executed with the
    argument 'test'. The tests are performed exclusively on data for the species
    1H. Before the tests can be run, the data for 1H must be generated, and the
    mass correction must not be included. For detials, please see the header of
    tests.c and the section 'Testing' in the README contained in the top-level
    directory of AlkCalc.

    Of course, make sure the library functions of AlkCalc are compiled before
    trying to test them.
