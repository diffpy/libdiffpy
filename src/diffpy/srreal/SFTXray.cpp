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
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class SFTXray
*
* Implementation of x-ray ScatteringFactorTable using Waasmaier and Kirfel
* approximation.  The instance can be also created by calling the
* createByType("xray") factory.
*
*****************************************************************************/

#include <diffpy/srreal/SFTXray.hpp>
#include <diffpy/srreal/scatteringfactordata.hpp>
#include <diffpy/serialization.ipp>

namespace diffpy {
namespace srreal {

using namespace std;

// Public Methods ------------------------------------------------------------

// HasClassRegistry methods

ScatteringFactorTablePtr SFTXray::create() const
{
    ScatteringFactorTablePtr rv(new SFTXray());
    return rv;
}

ScatteringFactorTablePtr SFTXray::clone() const
{
    ScatteringFactorTablePtr rv(new SFTXray(*this));
    return rv;
}

const string& SFTXray::type() const
{
    static string rv = "xray";
    return rv;
}

// own methods - overloads

const string& SFTXray::radiationType() const
{
    static string rv = "X";
    return rv;
}


double SFTXray::standardLookup(const string& smbl, double q) const
{
    return fxrayatq(smbl, q);
}

// Registration --------------------------------------------------------------

bool reg_SFTXray = SFTXray().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::SFTXray)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::SFTXray)
