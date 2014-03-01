/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class SFTElectronNumber
*
* ScatteringFactorTable implementation where scattering power equals number
* of valence electrons and is constant with Q.
*
*****************************************************************************/

#include <diffpy/srreal/SFTElectronNumber.hpp>
#include <diffpy/srreal/scatteringfactordata.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Public Methods ------------------------------------------------------------

// HasClassRegistry methods

ScatteringFactorTablePtr SFTElectronNumber::create() const
{
    ScatteringFactorTablePtr rv(new SFTElectronNumber());
    return rv;
}


ScatteringFactorTablePtr SFTElectronNumber::clone() const
{
    ScatteringFactorTablePtr rv(new SFTElectronNumber(*this));
    return rv;
}

const string& SFTElectronNumber::type() const
{
    static string rv = "electronnumber";
    return rv;
}

// own methods - overloads

const string& SFTElectronNumber::radiationType() const
{
    static string rv = "EN";
    return rv;
}


double SFTElectronNumber::standardLookup(const string& smbl, double q) const
{
    return electronnumber(smbl);
}

// Registration --------------------------------------------------------------

bool reg_SFTElectronNumber = SFTElectronNumber().registerThisType() &&
        ScatteringFactorTable::aliasType("electronnumber", "EN");

}   // namespace srreal
}   // namespace diffpy
