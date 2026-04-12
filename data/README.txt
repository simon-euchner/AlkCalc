/* -------------------------------------------------------------------------- *
 * Eigenenergies and discretisation data (possibly also radial eigenstates)   *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */


Table containing information on the data files. The symbol X is a placeholder
for a species identifier, i.e., X=6LI, X=40CA+, etc., l is the orbital angular
momentum quantum number, 2j represents the value 2*j, where j is the total
angular momentum quantum number, and n is the principal quantum number.

--------------------------------------------------------------------------------
Quantity              Naming convention        Description
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Mesh                  discretisation-X.dat     Mesh points and corresponding
                                               step sizes, used to discretise
                                               the interval [0, rmax].

Eigenenergies         energies-X-l-2j.dat      Eigenenergies up to the
                                               user-defined maximal principal
                                               quantum number, for fixed l and
                                               j.

Radial eigenstate     state-X-n-l-2j.dat       Radial eigenstate with quantum
                                               numbers n, l, j.
--------------------------------------------------------------------------------


Notes.

    Depending on the user-defined variable PATH_TO_STATES in
    'interface.d/settings.h', the radial eigenstates are stored here as well.

    All data files are simple text files.
