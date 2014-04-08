/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class BondCalculator -- bond distance calculator
*
*****************************************************************************/

#include <cassert>
#include <cmath>
#include <sstream>

#include <diffpy/srreal/BondCalculator.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/serialization.ipp>

using namespace std;

namespace diffpy {
namespace srreal {

// Local Constants -----------------------------------------------------------

namespace {

const double DEFAULT_BONDCALCULATOR_RMAX = 5.0;

enum {
    DISTANCE_OFFSET,
    SITE0_OFFSET,
    SITE1_OFFSET,
    DIRECTION0_OFFSET,
    DIRECTION1_OFFSET,
    DIRECTION2_OFFSET,
    CHUNK_SIZE,
};

bool pchunks_compare(QuantityType::const_iterator p0,
        QuantityType::const_iterator p1)
{
    bool rv = lexicographical_compare(
            p0, p0 + CHUNK_SIZE, p1, p1 + CHUNK_SIZE);
    return rv;
}

typedef vector<QuantityType::const_iterator> ChunkPtrVector;


ChunkPtrVector chunks_split(const QuantityType& pv)
{
    assert(0 == pv.size() % CHUNK_SIZE);
    ChunkPtrVector rv(pv.size() / CHUNK_SIZE);
    QuantityType::const_iterator src = pv.begin();
    ChunkPtrVector::iterator dst = rv.begin();
    for (; src != pv.end(); src += CHUNK_SIZE, ++dst)  *dst = src;
    return rv;
}


void
chunks_merge(const ChunkPtrVector& chunks, QuantityType& datadest)
{
    ChunkPtrVector::const_iterator xc = chunks.begin();
    QuantityType tdest(CHUNK_SIZE * chunks.size());
    QuantityType::iterator dst = tdest.begin();
    for (; xc != chunks.end(); ++xc, dst += CHUNK_SIZE)
    {
        copy(*xc, (*xc) + CHUNK_SIZE, dst);
    }
    datadest.swap(tdest);
}


void
chunks_merge(const QuantityType& data0, const QuantityType& data1,
        QuantityType& datadest)
{
    ChunkPtrVector p0 = chunks_split(data0);
    ChunkPtrVector p1 = chunks_split(data1);
    ChunkPtrVector ptot(p0.size() + p1.size());
    merge(p0.begin(), p0.end(), p1.begin(), p1.end(), ptot.begin());
    chunks_merge(ptot, datadest);
}


}   // namespace

// Constructor ---------------------------------------------------------------

BondCalculator::BondCalculator()
{
    this->setRmax(DEFAULT_BONDCALCULATOR_RMAX);
    mevaluator->setFlag(USEFULLSUM, true);
    mevaluator->setFlag(FIXEDSITEINDEX, true);
}

// Public Methods ------------------------------------------------------------

QuantityType BondCalculator::distances() const
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


vector<R3::Vector> BondCalculator::directions() const
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


vector<int> BondCalculator::sites0() const
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


vector<int> BondCalculator::sites1() const
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


vector<string> BondCalculator::types0() const
{
    return siteIndicesToTypes(mstructure, this->sites0());
}


vector<string> BondCalculator::types1() const
{
    return siteIndicesToTypes(mstructure, this->sites1());
}


void BondCalculator::filterCone(R3::Vector coneaxis, double degrees)
{
    using namespace diffpy::validators;
    double nmconeaxis = R3::norm(coneaxis);
    ensureEpsilonPositive("magnitude of cone vector", nmconeaxis);
    coneaxis /= nmconeaxis;
    mfilter_directions.push_back(coneaxis);
    mfilter_degrees.push_back(degrees);
}


void BondCalculator::filterOff()
{
    mfilter_directions.clear();
    mfilter_degrees.clear();
}

// PairQuantity overloads

string BondCalculator::getParallelData() const
{
    ostringstream storage(ios::binary);
    diffpy::serialization::oarchive oa(storage, ios::binary);
    oa << mpairspop << mpairsadd;
    return storage.str();
}

// Protected Methods ---------------------------------------------------------

void BondCalculator::resetValue()
{
    mvalue.clear();
    mpairspop.clear();
    mpairsadd.clear();
    this->PairQuantity::resetValue();
}


void BondCalculator::addPairContribution(
        const BaseBondGenerator& bnds,
        int summationscale)
{
    assert(summationscale == +1 || summationscale == -1);
    static R3::Vector ru01;
    const R3::Vector& r01 = bnds.r01();
    ru01 = r01 / bnds.distance();
    if (!(this->checkConeFilters(ru01)))  return;
    double entry[CHUNK_SIZE];
    entry[DISTANCE_OFFSET] = bnds.distance();
    entry[SITE0_OFFSET] = bnds.site0();
    entry[SITE1_OFFSET] = bnds.site1();
    entry[DIRECTION0_OFFSET] = r01[0];
    entry[DIRECTION1_OFFSET] = r01[1];
    entry[DIRECTION2_OFFSET] = r01[2];
    QuantityType& pv = (summationscale > 0) ? mpairsadd : mpairspop;
    pv.insert(pv.end(), entry, entry + CHUNK_SIZE);
}


void BondCalculator::executeParallelMerge(const std::string& pdata)
{
    istringstream storage(pdata, ios::binary);
    diffpy::serialization::iarchive ia(storage, ios::binary);
    QuantityType parpop, paradd;
    ia >> parpop >> paradd;
    chunks_merge(mpairspop, parpop, mpairspop);
    chunks_merge(mpairsadd, paradd, mpairsadd);
}


void BondCalculator::finishValue()
{
    // filter-out entries marked for removal
    assert(mpairspop.size() <= mvalue.size());
    ChunkPtrVector ppop = chunks_split(mpairspop);
    sort(ppop.begin(), ppop.end(), pchunks_compare);
    ChunkPtrVector padd = chunks_split(mpairsadd);
    sort(padd.begin(), padd.end(), pchunks_compare);
    if (mevaluator->isParallel())
    {
        chunks_merge(ppop, mpairspop);
        chunks_merge(padd, mpairsadd);
        return;
    }
    ChunkPtrVector pval = chunks_split(mvalue);
    ChunkPtrVector::iterator last;
    last = set_difference(pval.begin(), pval.end(),
            ppop.begin(), ppop.end(), pval.begin(), pchunks_compare);
    pval.erase(last, pval.end());
    ChunkPtrVector ptot(pval.size() + padd.size());
    merge(pval.begin(), pval.end(), padd.begin(), padd.end(), ptot.begin());
    chunks_merge(ptot, mvalue);
    mpairspop.clear();
    mpairsadd.clear();
}


void BondCalculator::stashPartialValue()
{
    mstashedvalue.swap(mvalue);
}


void BondCalculator::restorePartialValue()
{
    mvalue.swap(mstashedvalue);
    mstashedvalue.clear();
}

// Private Methods -----------------------------------------------------------

int BondCalculator::count() const
{
    int rv = int(mvalue.size()) / CHUNK_SIZE;
    return rv;
}


bool BondCalculator::checkConeFilters(const R3::Vector& ru01) const
{
    using diffpy::mathutils::eps_eq;
    if (mfilter_directions.empty())     return true;
    assert(eps_eq(1.0, R3::norm(ru01)));
    vector<R3::Vector>::const_iterator coneax = mfilter_directions.begin();
    vector<double>::const_iterator deg = mfilter_degrees.begin();
    for (;  coneax != mfilter_directions.end(); ++coneax, ++deg)
    {
        if (180.0 <= *deg)    return true;
        double angledegrees = 180.0 / M_PI * acos(R3::dot(ru01, *coneax));
        if (angledegrees <= *deg)    return true;
    }
    return false;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::BondCalculator)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::BondCalculator)

// End of file
