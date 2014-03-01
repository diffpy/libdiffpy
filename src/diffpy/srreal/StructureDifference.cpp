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
* class StructureDifference for storing differences between two
* instances of the StructureAdapter class.
*
* StructureDifference is returned by the StructureAdapter::diff method.
*
*****************************************************************************/

#include <diffpy/srreal/StructureDifference.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class StructureDifference
//////////////////////////////////////////////////////////////////////////////

// Constructors --------------------------------------------------------------

StructureDifference::StructureDifference() : diffmethod(Method::NONE)
{ }


StructureDifference::StructureDifference(
        StructureAdapterConstPtr oldstru,
        StructureAdapterConstPtr newstru)
    :
    stru0(oldstru), stru1(newstru), diffmethod(Method::NONE)
{ }

// Public Methods ------------------------------------------------------------

bool StructureDifference::allowsfastupdate() const
{
    int N0 = stru0 ? stru0->countSites() : 0;
    double popbound = (1 - sqrt(0.5)) * N0;
    return int(pop0.size()) < popbound;
}

}   // namespace srreal
}   // namespace diffpy
