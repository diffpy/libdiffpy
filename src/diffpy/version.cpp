/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2014 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Definitions for the libdiffpy_version_info constants.
*
*****************************************************************************/

#include <diffpy/version.hpp>

const long long libdiffpy_version_info::version = DIFFPY_VERSION;
const char* libdiffpy_version_info::version_str = DIFFPY_VERSION_STR;
const int libdiffpy_version_info::major = DIFFPY_VERSION_MAJOR;
const int libdiffpy_version_info::minor = DIFFPY_VERSION_MINOR;
const int libdiffpy_version_info::micro = DIFFPY_VERSION_MICRO;
const int libdiffpy_version_info::patch = DIFFPY_VERSION_PATCH;
const char* libdiffpy_version_info::date = DIFFPY_VERSION_DATE;
const char* libdiffpy_version_info::git_sha = DIFFPY_GIT_SHA;

// End of file
