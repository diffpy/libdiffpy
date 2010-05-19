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

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/HasClassRegistry.ipp>

using namespace std;
using diffpy::srreal::ScatteringFactorTable;

// Unique instantiation of the template registry base class.
template class HasClassRegistry<ScatteringFactorTable>;

namespace diffpy {
namespace srreal {

// class ScatteringFactorTable -----------------------------------------------

// public methods

double ScatteringFactorTable::lookup(const string& smbl) const
{
    using namespace std;
    map<string, double>::const_iterator isft;
    isft = mtable.find(smbl);
    if (isft == mtable.end())
    {
        double value = this->lookupatq(smbl, 0.0);
        mtable[smbl] = value;
        isft = mtable.find(smbl);
    }
    return isft->second;
}


double ScatteringFactorTable::lookupatq(const std::string& smbl,
        double q) const
{
    return this->lookup(smbl);
}


void ScatteringFactorTable::setCustom(const string& smbl, double value)
{
    mtable[smbl] = value;
}


void ScatteringFactorTable::resetCustom(const string& smbl)
{
    mtable.erase(smbl);
}


void ScatteringFactorTable::resetAll()
{
    mtable.clear();
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


ScatteringFactorTablePtr
ScatteringFactorTableOwner::getScatteringFactorTable()
{
    assert(msftable.get());
    return msftable;
}


const ScatteringFactorTablePtr
ScatteringFactorTableOwner::getScatteringFactorTable() const
{
    assert(msftable.get());
    return msftable;
}


const string& ScatteringFactorTableOwner::getRadiationType() const
{
    const string& tp = this->getScatteringFactorTable()->radiationType();
    return tp;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
