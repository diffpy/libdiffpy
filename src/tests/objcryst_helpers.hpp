/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* shared helpers for loading diffpy structure objects.
*
*****************************************************************************/

#ifndef OBJCRYST_HELPERS_HPP_INCLUDED
#define OBJCRYST_HELPERS_HPP_INCLUDED

#include <string>
#include <ObjCryst/ObjCryst/Crystal.h>

ObjCryst::Crystal* loadTestCrystal(const std::string& tailname);

#endif  // OBJCRYST_HELPERS_HPP_INCLUDED
