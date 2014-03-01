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
* class SFTElectron
*
* Implementation of electron ScatteringFactorTable using the formula in
* International Tables Volume C, page 224.
* The instance can be also created by calling the
* createByType("electron") factory.
*
*****************************************************************************/

#include <diffpy/srreal/SFTElectron.hpp>
#include <diffpy/srreal/scatteringfactordata.hpp>

namespace diffpy {
namespace srreal {

using namespace std;

// Public Methods ------------------------------------------------------------

// HasClassRegistry methods

ScatteringFactorTablePtr SFTElectron::create() const
{
    ScatteringFactorTablePtr rv(new SFTElectron());
    return rv;
}

ScatteringFactorTablePtr SFTElectron::clone() const
{
    ScatteringFactorTablePtr rv(new SFTElectron(*this));
    return rv;
}

const string& SFTElectron::type() const
{
    static string rv = "electron";
    return rv;
}

// own methods - overloads

const string& SFTElectron::radiationType() const
{
    static string rv = "E";
    return rv;
}


double SFTElectron::standardLookup(const string& smbl, double q) const
{
    return felectronatq(smbl, q);
}

// Registration --------------------------------------------------------------

bool reg_SFTElectron = (
        SFTElectron().registerThisType() &&
        ScatteringFactorTable::aliasType("electron", "E")
        );

}   // namespace srreal
}   // namespace diffpy
