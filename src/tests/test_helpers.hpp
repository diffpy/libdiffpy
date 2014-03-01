/*****************************************************************************
*
* diffpy.lib        by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Helper functions used in other unit tests.
*
*****************************************************************************/

#ifndef TEST_HELPERS_HPP_INCLUDED
#define TEST_HELPERS_HPP_INCLUDED

#include <string>
#include <diffpy/srreal/forwardtypes.hpp>

std::string prepend_tests_dir(const std::string& f);
std::string prepend_testdata_dir(const std::string& f);

diffpy::srreal::StructureAdapterPtr
    loadTestPeriodicStructure(const std::string& tailname);

#endif  // TEST_HELPERS_HPP_INCLUDED
