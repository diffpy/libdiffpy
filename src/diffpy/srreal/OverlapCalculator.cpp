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

#include <diffpy/srreal/OverlapCalculator.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/serialization.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Local Constants -----------------------------------------------------------

namespace {

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

OverlapCalculator::OverlapCalculator()
{
    // default configuration
    AtomRadiiTablePtr table(new AtomRadiiTable);
    this->setAtomRadiiTable(table);
    // attributes
    this->registerDoubleAttribute("rmaxused", this,
            &OverlapCalculator::getRmaxUsed);
};

// Public Methods ------------------------------------------------------------

QuantityType OverlapCalculator::overlaps() const
{
    int n = this->count();
    QuantityType rv;
    rv.reserve(n);
    for (int idx = 0; idx < n; ++idx)
    {
        double olp = this->suboverlap(idx);
        if (olp >= 0.0)  rv.push_back(olp);
    }
    return rv;
}


QuantityType OverlapCalculator::distances() const
{
    return this->subvector(DISTANCE_OFFSET, OVERLAPPING);
}


vector<R3::Vector> OverlapCalculator::directions() const
{
    int n = this->count();
    vector<R3::Vector> rv;
    rv.reserve(n);
    for (int index = 0; index < n; ++index)
    {
        if (this->suboverlap(index) <= 0.0)  continue;
        rv.push_back(this->subdirection(index));
    }
    return rv;
}


vector<int> OverlapCalculator::sites0() const
{
    QuantityType rv0 = this->subvector(SITE0_OFFSET, OVERLAPPING);
    vector<int> rv1(rv0.begin(), rv0.end());
    return rv1;
}


vector<int> OverlapCalculator::sites1() const
{
    QuantityType rv0 = this->subvector(SITE1_OFFSET, OVERLAPPING);
    vector<int> rv1(rv0.begin(), rv0.end());
    return rv1;
}


QuantityType OverlapCalculator::siteSquareOverlaps() const
{
    int n = this->count();
    int cntsites = this->countSites();
    QuantityType rv(cntsites, 0.0);
    for (int index = 0; index < n; ++index)
    {
        double olp = this->suboverlap(index);
        if (olp <= 0.0)  continue;
        double sqoverlap = olp * olp;
        int i = int(this->subvalue(SITE0_OFFSET, index));
        int j = int(this->subvalue(SITE1_OFFSET, index));
        int sumscale = (i == j) ? 1 : 2;
        rv[i] += sumscale * sqoverlap * mstructure->siteOccupancy(j);
        rv[j] += sumscale * sqoverlap * mstructure->siteOccupancy(i);
    }
    // Overlaps have been summed twice and need to be divided among 2 atoms.
    QuantityType::iterator xi;
    for (xi = rv.begin(); xi != rv.end(); ++xi)  *xi /= 4;
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


double OverlapCalculator::meanSquareOverlap() const
{
    double totocc = mstructure->totalOccupancy();
    double rv = (totocc > 0) ? (this->totalSquareOverlap() / totocc) : 0.0;
    return rv;
}


double OverlapCalculator::flipDiffTotal(int i, int j) const
{
    int cntsites = this->countSites();
    if (i < 0 || i >= cntsites || j < 0 || j >= cntsites)
    {
        const char* emsg = "Index out of range.";
        throw invalid_argument(emsg);
    }
    bool sameradii = (i == j) ||
        (mstructure_cache.siteradii[i] == mstructure_cache.siteradii[j]);
    if (sameradii)  return 0.0;
    // here we have to remove the overlap contributions for i and j
    double rv = 0.0;
    list<int> allids = this->getNeighborIds(i);
    const list<int>& jneighbors = this->getNeighborIds(j);
    allids.insert(allids.end(), jneighbors.begin(), jneighbors.end());
    list<int>::const_iterator idx;
    for (idx = allids.begin(); idx != allids.end(); ++idx)
    {
        int i1 = int(this->subvalue(SITE0_OFFSET, *idx));
        int j1 = int(this->subvalue(SITE1_OFFSET, *idx));
        if (i1 == j1)  continue;
        double sqscale =
            mstructure->siteOccupancy(i1) * mstructure->siteOccupancy(j1) *
            0.5 * (mstructure->siteMultiplicity(i1) +
                    mstructure->siteMultiplicity(j1));
        double olp0 = this->suboverlap(*idx);
        double olp1 = this->suboverlap(*idx, i, j);
        rv -= sqscale * olp0 * olp0;
        rv += sqscale * olp1 * olp1;
    }
    return rv;
}


double OverlapCalculator::flipDiffMean(int i, int j) const
{
    double totocc = mstructure->totalOccupancy();
    double rv = (totocc > 0) ? (this->flipDiffTotal(i, j) / totocc) : 0.0;
    return rv;
}


vector<R3::Vector> OverlapCalculator::gradients() const
{
    using diffpy::mathutils::eps_gt;
    int n = this->count();
    int cntsites = this->countSites();
    vector<R3::Vector> rv(cntsites, R3::Vector(0.0, 0.0, 0.0));
    R3::Vector gij;
    for (int index = 0; index < n; ++index)
    {
        double olp = this->suboverlap(index);
        if (olp <= 0.0)  continue;
        const double& dst = this->subvalue(DISTANCE_OFFSET, index);
        assert(eps_gt(dst, 0.0));
        int i = int(this->subvalue(SITE0_OFFSET, index));
        int j = int(this->subvalue(SITE1_OFFSET, index));
        gij = -2.0 * olp / dst * this->subdirection(index);
        rv[j] += gij;
        rv[i] -= gij;
    }
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
    assert(this->countSites() == int(mstructure_cache.siteradii.size()));
    double rv = min(this->getRmax(), mstructure_cache.maxseparation);
    return rv;
}

// Protected Methods ---------------------------------------------------------

void OverlapCalculator::resetValue()
{
    mvalue.clear();
    mneighborids.clear();
    this->cacheStructureData();
    this->PairQuantity::resetValue();
}


void OverlapCalculator::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmax(this->getRmaxUsed());
}


void OverlapCalculator::addPairContribution(
        const BaseBondGenerator& bnds, int summationscale)
{
    using diffpy::mathutils::eps_eq;
    assert(bnds.distance() <= mstructure_cache.maxseparation);
    if (eps_eq(0.0, bnds.distance()))    return;
    const R3::Vector& r01 = bnds.r01();
    int baseidx = mvalue.size();
    mvalue.insert(mvalue.end(), CHUNK_SIZE, 0.0);
    mvalue[baseidx + DISTANCE_OFFSET] = bnds.distance();
    mvalue[baseidx + DIRECTION0_OFFSET] = r01[0];
    mvalue[baseidx + DIRECTION1_OFFSET] = r01[1];
    mvalue[baseidx + DIRECTION2_OFFSET] = r01[2];
    mvalue[baseidx + SITE0_OFFSET] = bnds.site0();
    mvalue[baseidx + SITE1_OFFSET] = bnds.site1();
}


void OverlapCalculator::executeParallelMerge(const std::string& pdata)
{
    istringstream storage(pdata, ios::binary);
    diffpy::serialization::iarchive ia(storage, ios::binary);
    QuantityType pvalue;
    ia >> pvalue;
    mvalue.insert(mvalue.end(), pvalue.begin(), pvalue.end());
}

// Private Methods -----------------------------------------------------------

int OverlapCalculator::count() const
{
    int rv = int(mvalue.size()) / CHUNK_SIZE;
    return rv;
}


QuantityType OverlapCalculator::subvector(int offset,
        OverlapCalculator::OverlapFlag flag) const
{
    assert(0 <= offset && offset < CHUNK_SIZE);
    assert(flag == ALLVALUES || flag == OVERLAPPING);
    int n = this->count();
    QuantityType rv;
    rv.reserve(n);
    for (int index = 0; index < n; ++index)
    {
        if (OVERLAPPING == flag && this->suboverlap(index) <= 0.0)  continue;
        rv.push_back(this->subvalue(offset, index));
    }
    return rv;
}


const double& OverlapCalculator::subvalue(int offset, int index) const
{
    assert(0 <= offset && offset < CHUNK_SIZE);
    assert(0 <= index && index < this->count());
    return mvalue[offset + CHUNK_SIZE * index];
}


const R3::Vector& OverlapCalculator::subdirection(int index) const
{
    static R3::Vector rv;
    double* p = rv.data();
    *p++ = this->subvalue(DIRECTION0_OFFSET, index);
    *p++ = this->subvalue(DIRECTION1_OFFSET, index);
    *p++ = this->subvalue(DIRECTION2_OFFSET, index);
    return rv;
}


double OverlapCalculator::suboverlap(int index, int flipi, int flipj) const
{
    int cntsites = this->countSites();
    assert(0 <= flipi && flipi < cntsites);
    assert(0 <= flipj && flipj < cntsites);
    int i = int(this->subvalue(SITE0_OFFSET, index));
    int j = int(this->subvalue(SITE1_OFFSET, index));
    const double& radiusi = (flipi == flipj) ? mstructure_cache.siteradii[i] :
        (i == flipi) ? mstructure_cache.siteradii[flipj] :
        (i == flipj) ? mstructure_cache.siteradii[flipi] :
        mstructure_cache.siteradii[i];
    const double& radiusj = (flipi == flipj) ? mstructure_cache.siteradii[j] :
        (j == flipi) ? mstructure_cache.siteradii[flipj] :
        (j == flipj) ? mstructure_cache.siteradii[flipi] :
        mstructure_cache.siteradii[j];
    const double& dij = this->subvalue(DISTANCE_OFFSET, index);
    double sepij = radiusi + radiusj;
    double rv = (dij < sepij) ? (sepij - dij) : 0.0;
    return rv;
}


void OverlapCalculator::cacheStructureData()
{
    int cntsites = this->countSites();
    mstructure_cache.siteradii.resize(cntsites);
    const AtomRadiiTablePtr& table = this->getAtomRadiiTable();
    for (int i = 0; i < cntsites; ++i)
    {
        const string& smbl = mstructure->siteAtomType(i);
        mstructure_cache.siteradii[i] = table->lookup(smbl);
    }
    double maxradius = mstructure_cache.siteradii.empty() ?
        0.0 : *max_element(mstructure_cache.siteradii.begin(),
                mstructure_cache.siteradii.end());
    mstructure_cache.maxseparation = 2 * maxradius;
}


const std::list<int>&
OverlapCalculator::getNeighborIds(int k) const
{
    typedef std::list<int> NbList;
    int cntsites = this->countSites();
    assert(0 <= k && k < cntsites);
    if (mneighborids.empty())
    {
        mneighborids.resize(cntsites);
        int n = this->count();
        for (int idx = 0; idx < n; ++idx)
        {
            int i = int(this->subvalue(SITE0_OFFSET, idx));
            int j = int(this->subvalue(SITE1_OFFSET, idx));
            if (!mneighborids[i].get())  mneighborids[i].reset(new NbList);
            mneighborids[i]->push_back(idx);
            if (i == j)  continue;
            if (!mneighborids[j].get())  mneighborids[j].reset(new NbList);
            mneighborids[j]->push_back(idx);
        }
    }
    assert(cntsites == int(mneighborids.size()));
    static NbList noneighbors;
    const NbList& rv =
        mneighborids[k].get() ? *mneighborids[k] : noneighbors;
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
