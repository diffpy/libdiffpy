/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
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
* Small routines related to atom properties.
*
*****************************************************************************/

#ifndef ATOMUTILS_HPP_INCLUDED
#define ATOMUTILS_HPP_INCLUDED

#include <string>

namespace diffpy {
namespace srreal {

/// Return bare element symbol from atom symbol that may be isotope or ion,
/// for example "Ca2+" or "12-C".
std::string atomBareSymbol(const std::string& atomtype);

/// Return valence of possibly ionic symbol such as "S2-" or "Cl-".
int atomValence(const std::string& atomtype);

}   // namespace srreal
}   // namespace diffpy

#endif  // ATOMUTILS_HPP_INCLUDED
