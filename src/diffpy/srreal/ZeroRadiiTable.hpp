/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class ZeroRadiiTable -- concrete AtomRadiiTable that always returns zero
*
* $Id$
*
*****************************************************************************/

#ifndef ZERORADIITABLE_HPP_INCLUDED
#define ZERORADIITABLE_HPP_INCLUDED

#include <diffpy/srreal/AtomRadiiTable.hpp>

namespace diffpy {
namespace srreal {

class ZeroRadiiTable : public AtomRadiiTable
{
    public:

        // methods - HasClassRegistry
        SharedPtr create() const;
        SharedPtr clone() const;
        const std::string& type() const;
        // own methods
        double tableLookup(const std::string& smbl) const;
};

}   // namespace srreal
}   // namespace diffpy

#endif  // ZERORADIITABLE_HPP_INCLUDED
