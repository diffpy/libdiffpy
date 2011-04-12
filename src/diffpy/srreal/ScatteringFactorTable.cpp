/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class ScatteringFactorTable -- base class for looking up scattering factors
*
* $Id$
*
*****************************************************************************/

#include <boost/serialization/export.hpp>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/HasClassRegistry.ipp>

using namespace std;
using diffpy::srreal::ScatteringFactorTable;
using boost::unordered_map;

// Unique instantiation of the template registry base class.
template class HasClassRegistry<ScatteringFactorTable>;

namespace diffpy {
namespace srreal {

// class ScatteringFactorTable -----------------------------------------------

// public methods

double ScatteringFactorTable::lookup(const string& smbl, double q) const
{
    if (0.0 == q)
    {
        boost::unordered_map<std::string,double>::const_iterator ii =
            mqzerocache.find(smbl);
        if (ii != mqzerocache.end())   return ii->second;
    }
    CustomDataStorage::const_iterator csft = mcustom.find(smbl);
    double rv;
    if (csft != mcustom.end())
    {
        const string& srcsmbl = csft->second.first;
        const double& scale = csft->second.second;
        rv = this->standardLookup(srcsmbl, q) * scale;
    }
    else {
        rv = this->standardLookup(smbl, q);
    }
    if (0.0 == q)  mqzerocache[smbl] = rv;
    return rv;
}


void ScatteringFactorTable::setCustomFrom(
        const string& smbl, const string& srcsmbl,
        double value, double q)
{
    double fsrc = this->standardLookup(srcsmbl, q);
    double scale = value / fsrc;
    mcustom[smbl] = make_pair(srcsmbl, scale);
    mqzerocache.erase(smbl);
}


void ScatteringFactorTable::resetCustom(const string& smbl)
{
    mcustom.erase(smbl);
    mqzerocache.erase(smbl);
}


void ScatteringFactorTable::resetAll()
{
    mcustom.clear();
    mqzerocache.clear();
}


boost::unordered_set<string> ScatteringFactorTable::getCustomSymbols() const
{
    boost::unordered_set<string> rv;
    CustomDataStorage::const_iterator csft;
    for (csft = mcustom.begin(); csft != mcustom.end(); ++csft)
    {
        rv.insert(csft->first);
    }
    return rv;
}

// class ScatteringFactorTableOwner ------------------------------------------

void ScatteringFactorTableOwner::setScatteringFactorTable(
        ScatteringFactorTablePtr sft)
{
    msftable = sft;
}


void ScatteringFactorTableOwner::setScatteringFactorTableByType(
        const string& tp)
{
    msftable = ScatteringFactorTable::createByType(tp);
}


ScatteringFactorTablePtr&
ScatteringFactorTableOwner::getScatteringFactorTable()
{
    return msftable;
}


const ScatteringFactorTablePtr&
ScatteringFactorTableOwner::getScatteringFactorTable() const
{
    return msftable;
}


const string& ScatteringFactorTableOwner::getRadiationType() const
{
    static string empty;
    const string& tp = msftable.get() ? msftable->radiationType() : empty;
    return tp;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT(diffpy::srreal::ScatteringFactorTablePtr)
BOOST_CLASS_EXPORT(diffpy::srreal::ScatteringFactorTableOwner)

// End of file
