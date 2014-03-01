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
* class ScatteringFactorTable -- base class for looking up scattering factors
*
*****************************************************************************/

#include <boost/serialization/export.hpp>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/HasClassRegistry.ipp>
#include <diffpy/validators.hpp>
#include <diffpy/serialization.ipp>

using namespace std;
using diffpy::validators::ensureNonNull;

namespace diffpy {

// Unique instantiation of the template registry base class.
template class HasClassRegistry<srreal::ScatteringFactorTable>;

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
    double rv;
    CustomDataStorage::const_iterator csft = mcustom.find(smbl);
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


void ScatteringFactorTable::setCustomAs(
        const string& smbl, const string& srcsmbl)
{
    const double scale = 1.0;
    mcustom[smbl] = make_pair(srcsmbl, scale);
    mqzerocache.erase(smbl);
    mticker.click();
}


void ScatteringFactorTable::setCustomAs(
        const string& smbl, const string& srcsmbl,
        double value, double q)
{
    double fsrc = this->standardLookup(srcsmbl, q);
    double scale = value / fsrc;
    mcustom[smbl] = make_pair(srcsmbl, scale);
    mqzerocache.erase(smbl);
    mticker.click();
}


void ScatteringFactorTable::resetCustom(const string& smbl)
{
    if (mcustom.count(smbl))  mticker.click();
    mcustom.erase(smbl);
    mqzerocache.erase(smbl);
}


void ScatteringFactorTable::resetAll()
{
    if (!mcustom.empty())  mticker.click();
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
    ensureNonNull("ScatteringFactorTable", sft);
    if (msftable != sft)  sft->ticker().click();
    msftable = sft;
}


void ScatteringFactorTableOwner::setScatteringFactorTableByType(
        const string& tp)
{
    msftable = ScatteringFactorTable::createByType(tp);
    msftable->ticker().click();
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

DIFFPY_INSTANTIATE_PTR_SERIALIZATION(diffpy::srreal::ScatteringFactorTable)
DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::ScatteringFactorTableOwner)

// End of file
