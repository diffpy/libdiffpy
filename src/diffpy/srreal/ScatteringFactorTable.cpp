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
#include <diffpy/ClassRegistry.hpp>

using namespace std;
using diffpy::ClassRegistry;

namespace diffpy {
namespace srreal {

// class ScatteringFactorTable -----------------------------------------------

// public methods

const double& ScatteringFactorTable::lookup(const string& smbl) const
{
    using namespace std;
    map<string, double>::const_iterator isft;
    isft = mtable.find(smbl);
    if (isft == mtable.end())
    {
        double value = this->fetch(smbl);
        mtable[smbl] = value;
        isft = mtable.find(smbl);
    }
    return isft->second;
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
        const ScatteringFactorTable& sft)
{
    if (msftable.get() == &sft)   return;
    msftable.reset(sft.clone());
}


void ScatteringFactorTableOwner::setScatteringFactorTable(const string& tp)
{
    auto_ptr<ScatteringFactorTable> sft(createScatteringFactorTable(tp));
    this->setScatteringFactorTable(*sft);
}


ScatteringFactorTable&
ScatteringFactorTableOwner::getScatteringFactorTable()
{
    assert(msftable.get());
    return *msftable;
}


const ScatteringFactorTable&
ScatteringFactorTableOwner::getScatteringFactorTable() const
{
    assert(msftable.get());
    return *msftable;
}


const string& ScatteringFactorTableOwner::getRadiationType() const
{
    const string& tp = this->getScatteringFactorTable().radiationType();
    return tp;
}

// Factory Functions ---------------------------------------------------------

ScatteringFactorTable* createScatteringFactorTable(const string& tp)
{
    return ClassRegistry<ScatteringFactorTable>::create(tp);
}


bool registerScatteringFactorTable(const ScatteringFactorTable& ref)
{
    return ClassRegistry<ScatteringFactorTable>::add(ref);
}


bool aliasScatteringFactorTable(const string& tp, const string& al)
{
    return ClassRegistry<ScatteringFactorTable>::alias(tp, al);
}


set<string> getScatteringFactorTableTypes()
{
    return ClassRegistry<ScatteringFactorTable>::getTypes();
}

}   // namespace srreal
}   // namespace diffpy

// End of file
