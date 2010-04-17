/*****************************************************************************
*
* diffpy.lib        by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Michigan State University.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* shared helpers for loading diffpy structure objects.
*
* $Id$
*
*****************************************************************************/

#include <diffpy/PythonInterface.hpp>
#include "python_helpers.hpp"
#include "tests_dir.hpp"

using namespace std;
using namespace boost;

python::object loadTestStructure(const string& tailname)
{
    using namespace diffpy;
    initializePython();
    python::object Structure =
        importFromPyModule("diffpy.Structure", "Structure");
    python::object stru = Structure();
    string fp = prepend_testdata_dir(tailname);
    stru.attr("read")(fp);
    return stru;
}

// End of file
