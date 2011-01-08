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
* class BondDistanceCalculator -- bond distance calculator
*
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/BondDistanceCalculator.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Local Constants -----------------------------------------------------------

namespace {

const double DEFAULT_BONDDISTANCECALCULATOR_RMAX = 5.0;

enum {
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

BondDistanceCalculator::BondDistanceCalculator()
{
    this->setRmax(DEFAULT_BONDDISTANCECALCULATOR_RMAX);
}

// Public Methods ------------------------------------------------------------

void BondDistanceCalculator::mergeParallelValue(const QuantityType& pvalue)
{
    mvalue.insert(mvalue.end(), pvalue.begin(), pvalue.end());
    this->finishValue();
}


QuantityType BondDistanceCalculator::distances() const
{
    QuantityType rv;
    rv.reserve(this->count());
    QuantityType::const_iterator v = mvalue.begin() + DISTANCE_OFFSET;
    for (; v < mvalue.end(); v += CHUNK_SIZE)
    {
        rv.push_back(*v);
    }
    return rv;
}


vector<R3::Vector> BondDistanceCalculator::directions() const
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


vector<int> BondDistanceCalculator::sites0() const
{
    vector<int> rv;
    rv.reserve(this->count());
    QuantityType::const_iterator v = mvalue.begin() + SITE0_OFFSET;
    for (; v < mvalue.end(); v += CHUNK_SIZE)
    {
        rv.push_back(int(*v));
    }
    return rv;
}


vector<int> BondDistanceCalculator::sites1() const
{
    vector<int> rv;
    rv.reserve(this->count());
    QuantityType::const_iterator v = mvalue.begin() + SITE1_OFFSET;
    for (; v < mvalue.end(); v += CHUNK_SIZE)
    {
        rv.push_back(int(*v));
    }
    return rv;
}

// Protected Methods ---------------------------------------------------------

void BondDistanceCalculator::resetValue()
{
    mvalue.clear();
}


void BondDistanceCalculator::addPairContribution(
        const BaseBondGenerator& bnds,
        int summationscale)
{
    int baseidx = mvalue.size();
    mvalue.insert(mvalue.end(), CHUNK_SIZE, 0.0);
    mvalue[baseidx + DISTANCE_OFFSET] = bnds.distance();
    const R3::Vector r01 = bnds.r01();
    mvalue[baseidx + DIRECTION0_OFFSET] = r01[0];
    mvalue[baseidx + DIRECTION1_OFFSET] = r01[1];
    mvalue[baseidx + DIRECTION2_OFFSET] = r01[2];
    mvalue[baseidx + SITE0_OFFSET] = bnds.site0();
    mvalue[baseidx + SITE1_OFFSET] = bnds.site1();
}


void BondDistanceCalculator::finishValue()
{
    vector<QuantityType> chunks;
    QuantityType::const_iterator v0 = mvalue.begin();
    QuantityType w;
    for (; v0 < mvalue.end(); v0 += CHUNK_SIZE)
    {
        w.assign(v0, v0 + CHUNK_SIZE);
        chunks.push_back(w);
    }
    stable_sort(chunks.begin(), chunks.end());
    vector<QuantityType>::const_iterator wi = chunks.begin();
    QuantityType::iterator v2 = mvalue.begin();
    for (; wi != chunks.end(); ++wi, v2 += CHUNK_SIZE)
    {
        copy(wi->begin(), wi->end(), v2);
    }
}


int BondDistanceCalculator::count() const
{
    int rv = int(mvalue.size()) / CHUNK_SIZE;
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
