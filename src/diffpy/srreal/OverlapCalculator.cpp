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
* class OverlapCalculator -- calculator of atom radii overlaps
*
*****************************************************************************/

#include <algorithm>

#include <diffpy/srreal/OverlapCalculator.hpp>
#include <diffpy/srreal/ConstantRadiiTable.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/serialization.ipp>

using namespace std;
using diffpy::validators::ensureNonNull;

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


template <class T>
vector<T> vectorconvert(const QuantityType& v)
{
    vector<T> rv;
    rv.reserve(v.size());
    QuantityType::const_iterator x = v.begin();
    for (; x != v.end(); ++x)  rv.push_back(T(*x));
    return rv;
}


}   // namespace

// Constructor ---------------------------------------------------------------

OverlapCalculator::OverlapCalculator()
{
    mevaluator->setFlag(USEFULLSUM, true);
    // default configuration
    AtomRadiiTablePtr table(new ConstantRadiiTable);
    this->setAtomRadiiTable(table);
    this->cacheStructureData();
    // use very large rmax, it will be cropped by rmaxused
    this->setRmax(100);
    mneighborids_cached = false;
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
        if (olp <= 0.0)  continue;
        rv.push_back(olp);
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


SiteIndices OverlapCalculator::sites0() const
{
    QuantityType rv0 = this->subvector(SITE0_OFFSET, OVERLAPPING);
    return vectorconvert<int>(rv0);
}


SiteIndices OverlapCalculator::sites1() const
{
    QuantityType rv0 = this->subvector(SITE1_OFFSET, OVERLAPPING);
    return vectorconvert<int>(rv0);
}


vector<string> OverlapCalculator::types0() const
{
    return siteIndicesToTypes(mstructure, this->sites0());
}


vector<string> OverlapCalculator::types1() const
{
    return siteIndicesToTypes(mstructure, this->sites1());
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
        rv[i] += sqoverlap * mstructure->siteOccupancy(j);
    }
    // overlaps are shared by 2 atoms
    QuantityType::iterator xi;
    for (xi = rv.begin(); xi != rv.end(); ++xi)  *xi /= 2;
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
    const list<int>& ineighbors = this->getNeighborIds(i);
    const list<int>& jneighbors = this->getNeighborIds(j);
    unordered_set<int> allids;
    allids.insert(ineighbors.begin(), ineighbors.end());
    allids.insert(jneighbors.begin(), jneighbors.end());
    unordered_set<int>::const_iterator idx;
    for (idx = allids.begin(); idx != allids.end(); ++idx)
    {
        int i1 = int(this->subvalue(SITE0_OFFSET, *idx));
        int j1 = int(this->subvalue(SITE1_OFFSET, *idx));
        double sqscale =
            ((i1 == j1) ? 1 : 2) *
            mstructure->siteOccupancy(i1) * mstructure->siteOccupancy(j1) *
            mstructure->siteMultiplicity(i1) / 2;
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
        int j = int(this->subvalue(SITE1_OFFSET, index));
        gij = -2.0 * olp / dst * this->subdirection(index);
        rv[j] += gij;
    }
    return rv;
}


unordered_set<int> OverlapCalculator::getNeighborSites(int i) const
{
    unordered_set<int> rv;
    const list<int>& ineighbors = this->getNeighborIds(i);
    list<int>::const_iterator idx;
    for (idx = ineighbors.begin(); idx != ineighbors.end(); ++idx)
    {
        double olp = this->suboverlap(*idx);
        if (olp <= 0.0)  continue;
        assert(i == int(this->subvalue(SITE0_OFFSET, *idx)));
        int j1 = int(this->subvalue(SITE1_OFFSET, *idx));
        rv.insert(j1);
    }
    return rv;
}


QuantityType OverlapCalculator::coordinations() const
{
    int cntsites = this->countSites();
    QuantityType rv(cntsites, 0.0);
    int n = this->count();
    for (int index = 0; index < n; ++index)
    {
        double olp = this->suboverlap(index);
        if (olp <= 0.0)  continue;
        int j0 = int(this->subvalue(SITE0_OFFSET, index));
        int j1 = int(this->subvalue(SITE1_OFFSET, index));
        rv[j0] += mstructure->siteOccupancy(j1);
    }
    return rv;
}


unordered_map<string,double>
OverlapCalculator::coordinationByTypes(int i) const
{
    unordered_map<string,double> rv;
    const list<int>& ineighbors = this->getNeighborIds(i);
    list<int>::const_iterator idx;
    for (idx = ineighbors.begin(); idx != ineighbors.end(); ++idx)
    {
        double olp = this->suboverlap(*idx);
        if (olp <= 0.0)  continue;
        int j0 = int(this->subvalue(SITE0_OFFSET, *idx));
        int j1 = int(this->subvalue(SITE1_OFFSET, *idx));
        if (j0 == i)
        {
            const string& tp = mstructure->siteAtomType(j1);
            rv[tp] += mstructure->siteOccupancy(j1);
        }
        else
        {
            assert(j1 == i);
            const string& tp = mstructure->siteAtomType(j0);
            rv[tp] += mstructure->siteOccupancy(j0) *
                mstructure->siteMultiplicity(j0) /
                mstructure->siteMultiplicity(j1);
        }
    }
    return rv;
}


vector< unordered_set<int> > OverlapCalculator::neighborhoods() const
{
    int cntsites = this->countSites();
    typedef unordered_set<int> SiteSet;
    typedef std::vector< boost::shared_ptr<SiteSet> > SiteSetPointers;
    SiteSetPointers rvptr(cntsites);
    int n = this->count();
    for (int index = 0; index < n; ++index)
    {
        double olp = this->suboverlap(index);
        if (olp <= 0.0)  continue;
        int j0 = int(this->subvalue(SITE0_OFFSET, index));
        int j1 = int(this->subvalue(SITE1_OFFSET, index));
        if (!rvptr[j0].get())
        {
            rvptr[j0].reset(new SiteSet);
            rvptr[j0]->insert(j0);
        }
        if (!rvptr[j1].get())
        {
            rvptr[j0]->insert(j1);
            rvptr[j1] = rvptr[j0];
        }
        assert(rvptr[j0].get() && rvptr[j1].get());
        if (rvptr[j0].get() == rvptr[j1].get())  continue;
        // here we need to unify the j0 and j1 sets
        SiteSetPointers::value_type ss1 = rvptr[j1];
        SiteSet::const_iterator idx1 = ss1->begin();
        for (; idx1 != ss1->end(); ++idx1)
        {
            rvptr[j0]->insert(*idx1);
            rvptr[*idx1] = rvptr[j0];
        }
    }
    // helper function that initializes neighbor set at site i
    auto initsiteset = [&](int i, int& counter) {
        if (!rvptr[i]) {
            rvptr[i].reset(new SiteSet(&i, &i + 1));
            ++counter;
        }
    };
    // count initialized items in rvptr
    int cnt = count_if(rvptr.begin(), rvptr.end(),
            [](SiteSetPointers::value_type& p) { return bool(p); });
    // create self-neighborhoods unless prohibited by mask
    for (int j0 = 0; j0 < cntsites && cnt < cntsites; ++j0)
    {
        for (int j1 = j0; j1 < cntsites; ++j1)
        {
            if (rvptr[j0] && rvptr[j1])  continue;
            if (!this->getPairMask(j0, j1))  continue;
            initsiteset(j0, cnt);
            initsiteset(j1, cnt);
        }
    }
    unordered_set<const SiteSet*> duplicate;
    vector< unordered_set<int> > rv;
    SiteSetPointers::const_iterator ssp;
    for (ssp = rvptr.begin(); ssp != rvptr.end(); ++ssp)
    {
        if (!ssp->get() || duplicate.count(ssp->get()))  continue;
        duplicate.insert(ssp->get());
        rv.push_back(**ssp);
    }
    return rv;
}


void OverlapCalculator::setAtomRadiiTable(AtomRadiiTablePtr table)
{
    ensureNonNull("AtomRadiiTable", table);
    matomradiitable = table;
}


void OverlapCalculator::setAtomRadiiTableByType(const std::string& tp)
{
    matomradiitable = AtomRadiiTable::createByType(tp);
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
    mneighborids_cached = false;
    this->cacheStructureData();
    this->PairQuantity::resetValue();
}


void OverlapCalculator::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmin(this->getRmin());
    bnds.setRmax(this->getRmaxUsed());
}


void OverlapCalculator::addPairContribution(
        const BaseBondGenerator& bnds, int summationscale)
{
    assert(summationscale == 1);
    assert(bnds.distance() <= mstructure_cache.maxseparation);
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
    rv[0] = this->subvalue(DIRECTION0_OFFSET, index);
    rv[1] = this->subvalue(DIRECTION1_OFFSET, index);
    rv[2] = this->subvalue(DIRECTION2_OFFSET, index);
    return rv;
}


double OverlapCalculator::suboverlap(int index, int flipi, int flipj) const
{
    assert(0 <= flipi && flipi < this->countSites());
    assert(0 <= flipj && flipj < this->countSites());
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
    assert(0 <= k && k < this->countSites());
    if (!mneighborids_cached)
    {
        assert(mneighborids.empty());
        int n = this->count();
        for (int idx = 0; idx < n; ++idx)
        {
            int i = int(this->subvalue(SITE0_OFFSET, idx));
            mneighborids[i].push_back(idx);
        }
        mneighborids_cached = true;
    }
    static NbList noneighbors;
    NeighborIdsStorage::const_iterator nbit = mneighborids.find(k);
    const NbList& rv =
        (nbit == mneighborids.end()) ?  noneighbors : nbit->second;
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::OverlapCalculator)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::OverlapCalculator)

// End of file
