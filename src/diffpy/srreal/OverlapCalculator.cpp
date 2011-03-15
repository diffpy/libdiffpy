/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class OverlapCalculator -- calculator of atom radii overlaps
*
* $Id$
*
*****************************************************************************/

#include <algorithm>
#include <functional>
#include <diffpy/srreal/OverlapCalculator.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/mathutils.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Local Constants -----------------------------------------------------------

namespace {

enum {
    OVERLAP_OFFSET,
    DISTANCE_OFFSET,
    DIRECTION0_OFFSET,
    DIRECTION1_OFFSET,
    DIRECTION2_OFFSET,
    SITE0_OFFSET,
    SITE1_OFFSET,
    CHUNK_SIZE,
};

}   // namespace

// Constructor ---------------------------------------------------------------

OverlapCalculator::OverlapCalculator()
{
    // attributes
    this->registerDoubleAttribute("rmaxused", this,
            &OverlapCalculator::getRmaxUsed);
};

// Public Methods ------------------------------------------------------------

QuantityType OverlapCalculator::overlaps() const
{
    return this->subvalue(OVERLAP_OFFSET);
}


QuantityType OverlapCalculator::distances() const
{
    return this->subvalue(DISTANCE_OFFSET);
}


vector<R3::Vector> OverlapCalculator::directions() const
{
    vector<R3::Vector> rv;
    rv.reserve(this->count());
    QuantityType::const_iterator v = mvalue.begin();
    for (; v < mvalue.end(); v += CHUNK_SIZE)
    {
        R3::Vector dir(
                *(v + DIRECTION0_OFFSET),
                *(v + DIRECTION1_OFFSET),
                *(v + DIRECTION2_OFFSET));
        rv.push_back(dir);
    }
    return rv;
}


vector<int> OverlapCalculator::sites0() const
{
    QuantityType rv0 = this->subvalue(SITE0_OFFSET);
    vector<int> rv1(rv0.begin(), rv0.end());
    return rv1;
}


vector<int> OverlapCalculator::sites1() const
{
    QuantityType rv0 = this->subvalue(SITE1_OFFSET);
    vector<int> rv1(rv0.begin(), rv0.end());
    return rv1;
}


QuantityType OverlapCalculator::siteSquareOverlaps() const
{
    int cntsites = this->countSites();
    QuantityType rv(cntsites, 0.0);
    QuantityType olps = this->overlaps();
    vector<int> sts0 = this->sites0();
    vector<int> sts1 = this->sites1();
    assert(olps.size() == sts0.size());
    assert(olps.size() == sts1.size());
    QuantityType::const_iterator olpii = olps.begin();
    vector<int>::const_iterator ii0 = sts0.begin(), ii1= sts1.begin();
    for (; olpii != olps.end(); ++olpii, ++ii0, ++ii1)
    {
        double halfsqoverlap = 0.5 * (*olpii) * (*olpii);
        rv[*ii0] += halfsqoverlap * mstructure->siteOccupancy(*ii1);
        rv[*ii1] += halfsqoverlap * mstructure->siteOccupancy(*ii0);
    }
    return rv;
}


double OverlapCalculator::totalSquareOverlap() const
{
    int cntsites = this->countSites();
    QuantityType sqoverlap = this->siteSquareOverlaps();
    double rv = 0.0;
    for (int i = 0; i < cntsites; ++i)
    {
        rv += sqoverlap[i] * 
            mstructure->siteMultiplicity(i) * mstructure->siteOccupancy(i);
    }
    return rv;
}


double OverlapCalculator::totalFlipDiff(int i, int j) const
{
    // FIXME
    return 0.0;
}


vector<R3::Vector> OverlapCalculator::gradients() const
{
    // FIXME
    return vector<R3::Vector>();
}


double OverlapCalculator::msoverlap() const
{
    double totocc = mstructure->totalOccupancy();
    double rv = (totocc > 0) ? (this->totalSquareOverlap() / totocc) : 0.0;
    return rv;
}


double OverlapCalculator::rmsoverlap() const
{
    double rv = sqrt(this->msoverlap());
    return rv;
}


void OverlapCalculator::setAtomRadiiTable(AtomRadiiTablePtr table)
{
    matomradiitable = table;
}


AtomRadiiTablePtr& OverlapCalculator::getAtomRadiiTable()
{
    return matomradiitable;
}


const AtomRadiiTablePtr& OverlapCalculator::getAtomRadiiTable() const
{
    return matomradiitable;
}


double OverlapCalculator::getRmaxUsed() const
{
    int cntsites = this->countSites();
    assert(cntsites == int(mstructure_cache.siteradii.size()));
    double maxradius = mstructure_cache.siteradii.empty() ? 0.0 :
        *max_element(mstructure_cache.siteradii.begin(),
                mstructure_cache.siteradii.end());
    double rv = min(this->getRmax(), maxradius);
    return rv;
}



// Protected Methods ---------------------------------------------------------

void OverlapCalculator::resetValue()
{
    // FIXME
}


void OverlapCalculator::addPairContribution(const BaseBondGenerator&, int)
{
    // FIXME
}


void OverlapCalculator::executeParallelMerge(const QuantityType&)
{
    // FIXME
}


void OverlapCalculator::finishValue()
{
    // FIXME
}

// Private Methods -----------------------------------------------------------

int OverlapCalculator::count() const
{
    int rv = int(mvalue.size()) / CHUNK_SIZE;
    return rv;
}


QuantityType OverlapCalculator::subvalue(int offset) const
{
    assert(0 <= offset && offset < CHUNK_SIZE);
    QuantityType rv;
    rv.reserve(this->count());
    QuantityType::const_iterator v = mvalue.begin() + offset;
    for (; v < mvalue.end(); v += CHUNK_SIZE)  rv.push_back(*v);
    return rv;
}


void OverlapCalculator::cacheStructureData()
{
    int cntsites = this->countSites();
    mstructure_cache.siteradii.resize(cntsites);
    const AtomRadiiTable& table = *(this->getAtomRadiiTable());
    for (int i = 0; i < cntsites; ++i)
    {
        const string& smbl = mstructure->siteAtomType(i);
        mstructure_cache.siteradii[i] = table.lookup(smbl);
    }
}

}   // namespace srreal
}   // namespace diffpy

// End of file
