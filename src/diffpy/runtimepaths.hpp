/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Functions for resolving paths to static data files at runtime.
*
*****************************************************************************/

#ifndef RUNTIMEPATHS_HPP_INCLUDED
#define RUNTIMEPATHS_HPP_INCLUDED

#include <string>

namespace diffpy {
namespace runtimepaths {

/// Return full path to a data file included with DiffPy library.
///
/// The paths are looked up from the DIFFPYRUNTIME environment variable.
/// If DIFFPYRUNTIME does not exist, the paths are resolved as
/// (1) relative path of the datadir directory with respect to the libdir
///     directory at the installation time and if it does not exist as
/// (2) relative path of the runtime directory in a source tree
///     at the build time.
///
/// Throw runtime_error if the base directory cannot be found.
std::string datapath(const std::string& f);


}   // namespace runtimepaths
}   // namespace diffpy

#endif  // RUNTIMEPATHS_HPP_INCLUDED
