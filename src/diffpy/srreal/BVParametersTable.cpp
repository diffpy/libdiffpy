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
#include <string>

#include <diffpy/srreal/BVParametersTable.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Declarations --------------------------------------------------------------

extern const char* bvparm2009cif_lines[];

// Constructor ---------------------------------------------------------------

BVParametersTable::BVParametersTable()
{
    mstandardtable = this->getStandardSetOfBVParam();
}


BVParametersTable* BVParametersTable::copy() const
{
    BVParametersTable* rv = new BVParametersTable(*this);
    return rv;
}

// Public Methods ------------------------------------------------------------

const BVParam& BVParametersTable::lookup(const BVParam& bpk) const
{
    assert(mstandardtable);
    SetOfBVParam::const_iterator bpit;
    const BVParam* rv = NULL;
    if (!rv)
    {
        bpit = mcustomtable.find(bpk);
        if (bpit != mcustomtable.end())  rv = &(*bpit);
    }
    if (!rv)
    {
        bpit = mstandardtable->find(bpk);
        if (bpit != mstandardtable->end())  rv = &(*bpit);
    }
    if (!rv)
    {
        static BVParam bpnone;
        rv = &bpnone;
    }
    return *rv;
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

// End of file
