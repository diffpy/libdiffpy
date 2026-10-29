/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2013 Brookhaven Science Associates,
 *                   Brookhaven National Laboratory.
 *                   All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 ******************************************************************************
 *
 * Small routines related to atom properties.
 *
 *****************************************************************************/

#ifndef ATOMUTILS_HPP_INCLUDED
#define ATOMUTILS_HPP_INCLUDED

#include <string>
#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {

/// Return bare element symbol from atom symbol that may be isotope or ion,
/// for example "Ca2+" or "12-C".
DLL_EXPORT std::string atomBareSymbol(const std::string& atomtype);

/// Return valence of possibly ionic symbol such as "S2-" or "Cl-".
DLL_EXPORT int atomValence(const std::string& atomtype);

}  // namespace srreal
}  // namespace diffpy

#endif  // ATOMUTILS_HPP_INCLUDED
