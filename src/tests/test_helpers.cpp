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
* Definitions of helper functions for other unit tests.
*
*****************************************************************************/

#include "test_helpers.hpp"

// Resolve DIFFPYTESTSDIRPATH ------------------------------------------------

#ifndef DIFFPYTESTSDIRPATH
#error "Compiler must define the DIFFPYTESTSDIRPATH macro."
#endif
#define STRINGIFY(m) STRINGIFY_BRAIN_DAMAGE(m)
#define STRINGIFY_BRAIN_DAMAGE(m) #m

// Path utilities ------------------------------------------------------------

std::string prepend_tests_dir(const std::string& f)
{
    std::string rv = STRINGIFY(DIFFPYTESTSDIRPATH);
    rv = rv + '/' + f;
    return rv;
}


std::string prepend_testdata_dir(const std::string& f)
{
    using std::string;
    string ft = "testdata/";
    ft += f;
    string rv = prepend_tests_dir(ft);
    return rv;
}

// End of tests_dir.cpp
