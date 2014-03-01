/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
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
* class SFTElectron
*
* SFTElectron gives Q-dependent electron scattering factor according to
* the formula in International Tables Volume C, page 224.
* The formula diverges at Q = 0, where SFTElectron returns DOUBLE_MAX.
*
*****************************************************************************/

#ifndef SFTELECTRON_HPP_INCLUDED
#define SFTELECTRON_HPP_INCLUDED

#include <diffpy/srreal/ScatteringFactorTable.hpp>

namespace diffpy {
namespace srreal {

class SFTElectron : public ScatteringFactorTable
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

};  // class SFTElectron


}   // namespace srreal
}   // namespace diffpy

#endif  // SFTELECTRON_HPP_INCLUDED
