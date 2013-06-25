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
*****************************************************************************/

#ifndef TESTS_DIR_HPP_INCLUDED
#define TESTS_DIR_HPP_INCLUDED

#include <string>

std::string prepend_tests_dir(const std::string& f);
std::string prepend_testdata_dir(const std::string& f);

#endif  // TESTS_DIR_HPP_INCLUDED
