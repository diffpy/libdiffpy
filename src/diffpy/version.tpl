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
*   DIFFPY_VERSION_MAJOR,
*   DIFFPY_VERSION_MINOR,
*   DIFFPY_VERSION_STR,
*   DIFFPY_VERSION_DATE
*   DIFFPY_GIT_SHA
*
*****************************************************************************/

#ifndef VERSION_HPP_INCLUDED
#define VERSION_HPP_INCLUDED

#define DIFFPY_VERSION_MAJOR 1
#define DIFFPY_VERSION_MINOR 0

// DIFFPY_VERSION % 100000 is a version age in days since Unix time epoch
// DIFFPY_VERSION / 100000 % 100 is the minor version
// DIFFPY_VERSION / 10000000 is the major version

#define DIFFPY_VERSION ${DIFFPY_VERSION}

// DIFFPY_VERSION_STR is a string form of DIFFPY_VERSION

#define DIFFPY_VERSION_STR "${DIFFPY_VERSION_STR}"

// DIFFPY_VERSION_DATE is the subversion date of DIFFPY_VERSION

#define DIFFPY_VERSION_DATE "${DIFFPY_VERSION_DATE}"

// DIFFPY_GIT_SHA is an abbreviated git commit hash for the current version

#define DIFFPY_GIT_SHA "${DIFFPY_GIT_SHA}"

// Optional Features ---------------------------------------------------------

#if ${DIFFPY_HAS_OBJCRYST}
# define DIFFPY_HAS_OBJCRYST
#endif

#endif  // VERSION_HPP_INCLUDED

// vim:ft=cpp:
