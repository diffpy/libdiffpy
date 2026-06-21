/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
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
* Forward declaration of commonly used pointer types.
*
*****************************************************************************/

#ifndef FORWARDTYPES_HPP_INCLUDED
#define FORWARDTYPES_HPP_INCLUDED

#include <vector>
#include <memory>

namespace diffpy {
namespace srreal {

/// shared pointers to structure adapter and bond generator related objects
typedef std::shared_ptr<class BaseBondGenerator> BaseBondGeneratorPtr;
typedef std::shared_ptr<class StructureAdapter> StructureAdapterPtr;
typedef std::shared_ptr<const class StructureAdapter> StructureAdapterConstPtr;

/// type for holding indices for atom sites
typedef std::vector<int> SiteIndices;

}   // namespace srreal
}   // namespace diffpy

#endif  // FORWARDTYPES_HPP_INCLUDED
