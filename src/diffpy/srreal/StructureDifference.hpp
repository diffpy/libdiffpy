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
* class StructureDifference for storing differences between two
* instances of the StructureAdapter class.

* StructureDifference is returned by the StructureAdapter::diff method.
*
*****************************************************************************/

#ifndef STRUCTUREDIFFERENCE_HPP_INCLUDED
#define STRUCTUREDIFFERENCE_HPP_INCLUDED

#include <vector>
#include <diffpy/srreal/forwardtypes.hpp>

namespace diffpy {
namespace srreal {

class StructureDifference {

    public:

        // constructors
        StructureDifference()  { }
        StructureDifference(StructureAdapterConstPtr, StructureAdapterConstPtr);

        // data

        /// pointer to the original StructureAdapter
        StructureAdapterConstPtr stru0;
        /// pointer to the new StructureAdapter
        StructureAdapterConstPtr stru1;
        /// indices of atoms in stru0 that are not in stru1.
        /// These atoms need to be removed for a fast update of PairQuantity.
        std::vector<int> pop0;
        /// indices of atoms in stru1 that are not in stru0
        /// These atoms need to be added in a fast update of PairQuantity.
        std::vector<int> add1;

        // methods

        /// Return true if PairQuantity can be fast-updated after changing
        /// structure from stru0 to stru1.
        bool allowsfastupdate() const;

};

}   // namespace srreal
}   // namespace diffpy

#endif  // STRUCTUREDIFFERENCE_HPP_INCLUDED
