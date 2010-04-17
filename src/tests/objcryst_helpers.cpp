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

#include <iostream>
#include <fstream>
#include <sstream>

#include <ObjCryst/ObjCryst/CIF.h>

#include "objcryst_helpers.hpp"
#include "tests_dir.hpp"

using namespace std;


ObjCryst::Crystal* loadTestCrystal(const string& tailname)
{
    string fp = prepend_testdata_dir(tailname);
    ifstream in(fp.c_str());
    ObjCryst::CIF cif(in);
    // redirect cout to hide the CreateCrystalFromCIF chatty output
    streambuf* cout_sb = cout.rdbuf();
    ostringstream cifout;
    cout.rdbuf(cifout.rdbuf());
    ObjCryst::Crystal* crst = ObjCryst::CreateCrystalFromCIF(cif, false);
    // restore cout
    cout.rdbuf(cout_sb);
    return crst;
}

// End of file
