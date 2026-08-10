/* -------------------------------------------------------------------------- *
 * Eigenenergies and discretisation data (possibly also radial eigenstates)   *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */


Table containing information on the data files. The symbol X is a placeholder
for a species identifier, i.e., X = 6LI, X = 40CA+, etc., l is the orbital
angular momentum quantum number, 2j represents the integer value 2 * j, where j
is the half-integer total angular momentum quantum number, and n is the
principal quantum number.

--------------------------------------------------------------------------------
Quantity              Naming convention        Description
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Knot data             knotdata-X.dat           File including the knots (ti)
                                               (see theory/theory.pdf),
                                               including the multiplicities of
                                               the knots. The file also includes
                                               the non-zero step sizes hi = ti
                                               - tim1, where im1 means i - 1.

Eigenenergies         energies-X-l-2j.dat      Eigenenergies up to the
                                               user-defined maximal principal
                                               quantum number, for fixed l and
                                               j.

Radial eigenstate     state-X-n-l-2j.dat       Radial eigenstate with quantum
                                               numbers n, l, j.
--------------------------------------------------------------------------------


Notes.

    Depending on the user-defined variable PATH_TO_STATES in
    interface/settings.h, the radial eigenstates are stored here as well.

    All data files are simple text files.
