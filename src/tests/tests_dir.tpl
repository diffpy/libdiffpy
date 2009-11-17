/*****************************************************************************
*
* diffpy.lib        by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Michigan State University.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Path definitions of datafiles that are used in testing.
*
* $Id: globals.hpp 2706 2009-02-15 03:49:57Z juhas $
*
*****************************************************************************/

#include "tests_dir.hpp"

std::string prepend_tests_dir(const std::string& f)
{
    using namespace std;
    string rv("%(tests_dir)s");
    rv = rv + '/' + f;
    return rv;
}


std::string prepend_testdata_dir(const std::string& f)
{
    using namespace std;
    string ft = "testdata/";
    ft += f;
    string rv = prepend_tests_dir(ft);
    return rv;
}

// vim:ft=cpp:
// End of tests_dir.cpp
