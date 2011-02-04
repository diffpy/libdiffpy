/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class BVParametersTable -- table of bond valence sum parameters
*
* $Id$
*
*****************************************************************************/

#include <cassert>
#include <boost/serialization/export.hpp>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/BVParametersTable.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Declarations --------------------------------------------------------------

extern const char* bvparm2009cif_lines[];

// Static Methods ------------------------------------------------------------

const BVParam& BVParametersTable::none()
{
    static const BVParam bpnone;
    return bpnone;
}

// Constructor ---------------------------------------------------------------

BVParametersTable::BVParametersTable()
{
    mstandardtable = this->getStandardSetOfBVParam();
}


// Public Methods ------------------------------------------------------------

const BVParam& BVParametersTable::lookup(const BVParam& bpk) const
{
    assert(mstandardtable);
    SetOfBVParam::const_iterator bpit;
    bpit = mcustomtable.find(bpk);
    if (bpit != mcustomtable.end())  return *bpit;
    bpit = mstandardtable->find(bpk);
    if (bpit != mstandardtable->end())  return *bpit;
    // unspecified cations are marked with valence 9
    // perform this search only if the first atom is indeed a cation
    if (!(bpk.mvalence0 > 0))  return this->none();
    BVParam bpk9 = bpk;
    bpk9.mvalence0 = 9;
    bpit = mcustomtable.find(bpk9);
    if (bpit != mcustomtable.end())  return *bpit;
    bpit = mstandardtable->find(bpk9);
    if (bpit != mstandardtable->end())  return *bpit;
    // not found - return blank BVParam
    return this->none();
}


const BVParam& BVParametersTable::lookup(const string& atom0, int valence0,
        const string& atom1, int valence1) const
{
    BVParam bpk(atom0, valence0, atom1, valence1);
    const BVParam& rv = this->lookup(bpk);
    return rv;
}


void BVParametersTable::setCustom(const BVParam& bp)
{
    mcustomtable.insert(bp);
}


void BVParametersTable::resetCustom(const BVParam& bp)
{
    mcustomtable.erase(bp);
}


void BVParametersTable::resetCustom(const string& atom0, int valence0,
        const string& atom1, int valence1)
{
    BVParam bp(atom0, valence0, atom1, valence1);
    this->resetCustom(bp);
}


void BVParametersTable::resetAll()
{
    mcustomtable.clear();
}


BVParametersTable::SetOfBVParam
BVParametersTable::getAll() const
{
    SetOfBVParam rv = *mstandardtable;
    rv.insert(mcustomtable.begin(), mcustomtable.end());
    return rv;
}

// Private Methods -----------------------------------------------------------

BVParametersTable::SetOfBVParam*
BVParametersTable::getStandardSetOfBVParam() const
{
    static SetOfBVParam the_set;
    if (the_set.empty())
    {
        const char** cifline;
        for (cifline = bvparm2009cif_lines; *cifline != NULL; ++cifline)
        {
            BVParam bp;
            bp.setFromCifLine(*cifline);
            assert(!the_set.count(bp));
            the_set.insert(bp);
        }
    }
    return &the_set;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::BVParametersTable)
BOOST_CLASS_EXPORT(diffpy::srreal::BVParametersTable)

// End of file
