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

#include <cassert>
#include <cmath>

#include <diffpy/srreal/BondDistanceCalculator.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/mathutils.hpp>

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


void BondDistanceCalculator::filterCone(R3::Vector coneaxis, double degrees)
{
    using namespace diffpy::validators;
    double nmconeaxis = R3::norm(coneaxis);
    ensureEpsilonPositive("magnitude of cone vector", nmconeaxis);
    coneaxis /= nmconeaxis;
    mfilter_directions.push_back(coneaxis);
    mfilter_degrees.push_back(degrees);
}


void BondDistanceCalculator::filterOff()
{
    mfilter_directions.clear();
    mfilter_degrees.clear();
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
    using diffpy::mathutils::eps_eq;
    const R3::Vector& r01 = bnds.r01();
    if (eps_eq(0.0, bnds.distance()))    return;
    static R3::Vector ru01;
    ru01 = r01 / bnds.distance();
    if (!(this->checkConeFilters(ru01)))  return;
    int baseidx = mvalue.size();
    mvalue.insert(mvalue.end(), CHUNK_SIZE, 0.0);
    mvalue[baseidx + DISTANCE_OFFSET] = bnds.distance();
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
    sort(chunks.begin(), chunks.end());
    vector<QuantityType>::const_iterator wi = chunks.begin();
    QuantityType::iterator v2 = mvalue.begin();
    for (; wi != chunks.end(); ++wi, v2 += CHUNK_SIZE)
    {
        copy(wi->begin(), wi->end(), v2);
    }
}

// Private Methods -----------------------------------------------------------

int BondDistanceCalculator::count() const
{
    int rv = int(mvalue.size()) / CHUNK_SIZE;
    return rv;
}


bool BondDistanceCalculator::checkConeFilters(const R3::Vector& ru01) const
{
    using diffpy::mathutils::eps_eq;
    if (mfilter_directions.empty())     return true;
    assert(eps_eq(1.0, R3::norm(ru01)));
    vector<R3::Vector>::const_iterator coneax = mfilter_directions.begin();
    vector<double>::const_iterator deg = mfilter_degrees.begin();
    for (;  coneax != mfilter_directions.end(); ++coneax, ++deg)
    {
        double angledegrees = 180.0 / M_PI * acos(R3::dot(ru01, *coneax));
        if (angledegrees <= *deg)    return true;
    }
    return false;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
