/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class SFTElectronNumber
*
* ScatteringFactorTable with scattering power equal number
* of valence electrons without any Q-dependence
*
*****************************************************************************/

#ifndef SFTELECTRONNUMBER_HPP_INCLUDED
#define SFTELECTRONNUMBER_HPP_INCLUDED

#include <diffpy/srreal/ScatteringFactorTable.hpp>

namespace diffpy {
namespace srreal {

class SFTElectronNumber : public ScatteringFactorTable
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

};  // class SFTElectronNumber


}   // namespace srreal
}   // namespace diffpy

#endif  // SFTELECTRONNUMBER_HPP_INCLUDED
