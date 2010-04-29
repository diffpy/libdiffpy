/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
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
* Macro definitions for
*   DIFFPY_VERSION,
*   DIFFPY_VERSION_STR,
*   DIFFPY_VERSION_DATE
*
* $Id$
*
*****************************************************************************/

#ifndef VERSION_HPP_INCLUDED
#define VERSION_HPP_INCLUDED

// DIFFPY_VERSION % 100000 is the subversion revision number
// DIFFPY_VERSION / 100000 % 100 is the minor version
// DIFFPY_VERSION / 10000000 is the major version

#define DIFFPY_VERSION 10005143

// DIFFPY_VERSION_STR is a string form of DIFFPY_VERSION

#define DIFFPY_VERSION_STR "1.0-r5143"

// DIFFPY_VERSION_DATE is the subversion date of DIFFPY_VERSION

#define DIFFPY_VERSION_DATE "2010-04-28 02:30:48Z"

#endif  // VERSION_HPP_INCLUDED

// vim:ft=cpp:
