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

}   // namespace

class BondOp {

    public:

        static bool compare(
                const BondCalculator::BondEntry& be0,
                const BondCalculator::BondEntry& be1)
        {
            if (be0.distance < be1.distance)  return true;
            if (be0.distance > be1.distance)  return false;
            if (be0.site0 < be1.site0)  return true;
            if (be0.site0 > be1.site0)  return false;
            if (be0.site1 < be1.site1)  return true;
            if (be0.site1 > be1.site1)  return false;
            if (be0.direction0 < be1.direction0)  return true;
            if (be0.direction0 > be1.direction0)  return false;
            if (be0.direction1 < be1.direction1)  return true;
            if (be0.direction1 > be1.direction1)  return false;
            if (be0.direction2 < be1.direction2)  return true;
            if (be0.direction2 > be1.direction2)  return false;
            return false;
        }


        static bool reverse_compare(
                const BondCalculator::BondEntry& be0,
                const BondCalculator::BondEntry& be1)
        {
            return compare(be1, be0);
        }


        static void bmerge(
                BondCalculator::BondDataStorage& dstbonds,
                const BondCalculator::BondDataStorage& srcbonds)
        {
            dstbonds.resize(dstbonds.size() + srcbonds.size());
            merge(dstbonds.rbegin() + srcbonds.size(), dstbonds.rend(),
                    srcbonds.rbegin(), srcbonds.rend(),
                    dstbonds.rbegin(), reverse_compare);
        }


        static BondCalculator::BondEntry entryFrom(
                const BaseBondGenerator& bnds)
        {
            BondCalculator::BondEntry rv;
            rv.distance = bnds.distance();
            rv.site0 = bnds.site0();
            rv.site1 = bnds.site1();
            rv.direction0 = bnds.r01()[0];
            rv.direction1 = bnds.r01()[1];
            rv.direction2 = bnds.r01()[2];
            return rv;
        }


};  // class BondOp

// Constructor ---------------------------------------------------------------

BondCalculator::BondCalculator()
{
    this->setRmax(DEFAULT_BONDCALCULATOR_RMAX);
    this->setEvaluatorType(OPTIMIZED);
    mevaluator->setFlag(USEFULLSUM, true);
    mevaluator->setFlag(FIXEDSITEINDEX, true);
}

// Public Methods ------------------------------------------------------------

QuantityType BondCalculator::distances() const
{
    QuantityType rv;
    rv.reserve(this->count());
    BondDataStorage::const_iterator bi = mbonds.begin();
    for (; bi != mbonds.end(); bi++)  rv.push_back(bi->distance);
    return rv;
}


vector<R3::Vector> BondCalculator::directions() const
{
    vector<R3::Vector> rv;
    rv.reserve(this->count());
    BondDataStorage::const_iterator bi = mbonds.begin();
    for (; bi != mbonds.end(); ++bi)
    {
        R3::Vector dir(bi->direction0, bi->direction1, bi->direction2);
        rv.push_back(dir);
    }
    return rv;
}


SiteIndices BondCalculator::sites0() const
{
    SiteIndices rv;
    rv.reserve(this->count());
    BondDataStorage::const_iterator bi = mbonds.begin();
    for (; bi != mbonds.end(); ++bi)  rv.push_back(bi->site0);
    return rv;
}


SiteIndices BondCalculator::sites1() const
{
    SiteIndices rv;
    rv.reserve(this->count());
    BondDataStorage::const_iterator bi = mbonds.begin();
    for (; bi != mbonds.end(); ++bi)  rv.push_back(bi->site1);
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
    mticker.click();
}


void BondCalculator::filterOff()
{
    if (!mfilter_directions.empty())  mticker.click();
    mfilter_directions.clear();
    mfilter_degrees.clear();
}

// PairQuantity overloads

string BondCalculator::getParallelData() const
{
    ostringstream storage(ios::binary);
    diffpy::serialization::oarchive oa(storage, ios::binary);
    oa << mpopbonds << maddbonds;
    return storage.str();
}

// Protected Methods ---------------------------------------------------------

void BondCalculator::resetValue()
{
    mvalue.clear();
    mbonds.clear();
    maddbonds.clear();
    mpopbonds.clear();
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
    BondDataStorage& bes = (summationscale > 0) ? maddbonds : mpopbonds;
    bes.push_back(BondOp::entryFrom(bnds));
}


void BondCalculator::executeParallelMerge(const std::string& pdata)
{
    istringstream storage(pdata, ios::binary);
    diffpy::serialization::iarchive ia(storage, ios::binary);
    BondDataStorage bpop, badd;
    ia >> bpop >> badd;
    BondOp::bmerge(mpopbonds, bpop);
    BondOp::bmerge(maddbonds, badd);
}


void BondCalculator::finishValue()
{
    // filter-out entries marked for removal
    assert(mpopbonds.size() <= mbonds.size());
    sort(mpopbonds.begin(), mpopbonds.end(), BondOp::compare);
    sort(maddbonds.begin(), maddbonds.end(), BondOp::compare);
    if (mevaluator->isParallel())  return;
    BondDataStorage::iterator last;
    last = set_difference(mbonds.begin(), mbonds.end(),
            mpopbonds.begin(), mpopbonds.end(),
            mbonds.begin(), BondOp::compare);
    mbonds.erase(last, mbonds.end());
    if (mbonds.empty())  mbonds.swap(maddbonds);
    else  BondOp::bmerge(mbonds, maddbonds);
    mvalue = this->distances();
    mpopbonds.clear();
    maddbonds.clear();
}


void BondCalculator::stashPartialValue()
{
    mstashedvalue.bonds.swap(mbonds);
    mstashedvalue.popbonds.swap(mpopbonds);
    // No need to stash maddbonds as they are evaluated after partial value.
}


void BondCalculator::restorePartialValue()
{
    mbonds.swap(mstashedvalue.bonds);
    mpopbonds.swap(mstashedvalue.popbonds);
    mstashedvalue.bonds.clear();
    mstashedvalue.popbonds.clear();
}

// Private Methods -----------------------------------------------------------

int BondCalculator::count() const
{
    return mbonds.size();
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
