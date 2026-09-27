/* -------------------------------------------------------------------------- *
 * Testing of AlkCalc                                                         *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */


Notes.

    The file tests.c contains the code for all tests of AlkCalc, and the
    Makefile in this directory compiles and runs all of them at once, when
    executed with the argument 'test'. Before the tests can be run, make sure
    that

        (1) the library functions of AlkCalc are compiled (run the Makefile in
            the top-level directory of AlkCalc with the argument 'lib'),

        (2) the eigenenergies and radial eigenstates of the species 1H are
            generated for the pairs (l, j) = (0, 1/2), (1, 1/2), (1, 3/2),
            (2, 3/2), (2, 5/2), and (3, 7/2), and

        (3) the mass correction is NOT included (default behaviour of AlkCalc).

    The tests are performed exclusively on data for the species 1H (Hydrogen),
    because for the Hydrogen atom every reference value is known in closed
    analytical form. For details, please see the header of tests.c and the
    section 'Testing' in the README contained in the top-level directory of
    AlkCalc.

    Each test prints one line containing the outcome (PASS or FAIL), the tested
    quantity, the computed value (IS), the reference value (SHOULD BE), and the
    relative discrepancy (ERROR). Complex quantities are split into their real
    and imaginary parts (RE and IM), which are reported on two consecutive
    lines. A summary of all tests is printed at the end. Note that quantities
    involving states with l > 0 are tested with looser tolerances, because
    AlkCalc includes the LS coupling, whereas the analytical expressions for
    obtaining the reference values do not (see the header of the file tests.c).
