/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class SFTXray
*
* X-ray ScatteringFactorTable using Waasmaier and Kirfel approximation.
*
*****************************************************************************/

#ifndef SFTXRAY_HPP_INCLUDED
#define SFTXRAY_HPP_INCLUDED

#include <diffpy/srreal/ScatteringFactorTable.hpp>

namespace diffpy {
namespace srreal {

class SFTXray : public ScatteringFactorTable
{
    public:

        // HasClassRegistry methods
        ScatteringFactorTablePtr create() const;
        ScatteringFactorTablePtr clone() const;
        const std::string& type() const;
        // own methods
        const std::string& radiationType() const;
        // method overloads
        double standardLookup(const std::string& smbl, double q) const;

};  // class SFTXray


}   // namespace srreal
}   // namespace diffpy

#endif  // SFTXRAY_HPP_INCLUDED
