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

double ScatteringFactorTable::lookup(const string& smbl) const
{
    using namespace std;
    unordered_map<string, double>::const_iterator isft;
    isft = mtable.find(smbl);
    if (isft == mtable.end())
    {
        double value = this->lookupatq(smbl, 0.0);
        mtable[smbl] = value;
        isft = mtable.find(smbl);
    }
    return isft->second;
}


void ScatteringFactorTable::setCustom(const string& smbl, double value)
{
    mtable[smbl] = value;
    mcustomsymbols.insert(smbl);
}


void ScatteringFactorTable::resetCustom(const string& smbl)
{
    mtable.erase(smbl);
    mcustomsymbols.erase(smbl);
}


unordered_map<string,double> ScatteringFactorTable::getAllCustom() const
{
    unordered_map<string,double> rv;
    boost::unordered_set<string>::const_iterator smbl;
    for (smbl = mcustomsymbols.begin(); smbl != mcustomsymbols.end(); ++smbl)
    {
        rv[*smbl] = this->lookup(*smbl);
    }
    return rv;
}


void ScatteringFactorTable::resetAll()
{
    mtable.clear();
    mcustomsymbols.clear();
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
