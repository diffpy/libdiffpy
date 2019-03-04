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
* class SFTNeutron
*
* Implementation of neutron ScatteringFactorTable using coherent cross
* sections data from Paul Kienzle periodictable library for Python.
* approximation.  The instance can be also created by calling the
* createByType("neutron") factory.
*
*****************************************************************************/

#include <diffpy/srreal/SFTNeutron.hpp>
#include <diffpy/srreal/scatteringfactordata.hpp>
#include <diffpy/serialization.ipp>

namespace diffpy {
namespace srreal {

using namespace std;

// Public Methods ------------------------------------------------------------

// HasClassRegistry methods

ScatteringFactorTablePtr SFTNeutron::create() const
{
    ScatteringFactorTablePtr rv(new SFTNeutron());
    return rv;
}

ScatteringFactorTablePtr SFTNeutron::clone() const
{
    ScatteringFactorTablePtr rv(new SFTNeutron(*this));
    return rv;
}

const string& SFTNeutron::type() const
{
    static string rv = "neutron";
    return rv;
}

// own methods - overloads

const string& SFTNeutron::radiationType() const
{
    static string rv = "N";
    return rv;
}


double SFTNeutron::standardLookup(const string& smbl, double q) const
{
    return bcneutron(smbl);
}

// Registration --------------------------------------------------------------

bool reg_SFTNeutron = SFTNeutron().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::SFTNeutron)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::SFTNeutron)
