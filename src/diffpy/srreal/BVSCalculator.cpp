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
* class BVSCalculator -- concrete counter of pairs in a structure.
*
* $Id$
*
*****************************************************************************/

#include <cmath>
#include <cassert>

#include <diffpy/srreal/BVSCalculator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/srreal/BVParametersTable.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

BVSCalculator::BVSCalculator()
{
    // default configuration
    const double valence_precision = 1e-5;
    BVParametersTable bvtb;
    this->setBVParamTable(bvtb);
    this->setValencePrecision(valence_precision);
}

// Public Methods ------------------------------------------------------------

// results

QuantityType BVSCalculator::valences() const
{
    QuantityType rv(mstructure_cache.valences.size());
    copy(mstructure_cache.valences.begin(), mstructure_cache.valences.end(),
            rv.begin());
    return rv;
}


QuantityType BVSCalculator::bvdiff() const
{
    QuantityType vobs = this->valences();
    assert(vobs.size() == this->value().size());
    int cntsites = this->countSites();
    QuantityType rv(cntsites);
    const QuantityType& vsim = this->value();
    for (int i = 0; i < cntsites; ++i)
    {
        rv[i] = fabs(vobs[i]) - fabs(vsim[i]);
    }
    return rv;
}


double BVSCalculator::bvmsdiff() const
{
    QuantityType bd = this->bvdiff();
    assert(bd.size() == mstructure_cache.multiplicities.size());
    int cntsites = this->countSites();
    double sumofsquares = 0.0;
    for (int i = 0; i < cntsites; ++i)
    {
        sumofsquares += mstructure_cache.multiplicities[i] * bd[i] * bd[i];
    }
    double rv = (mstructure_cache.total_occupancy > 0.0) ?
        (sumofsquares / mstructure_cache.total_occupancy) : 0.0;
    return rv;
}


double BVSCalculator::bvrmsdiff() const
{
    double rvsq = this->bvmsdiff();
    return sqrt(rvsq);
}


void BVSCalculator::setBVParamTable(const BVParametersTable& bvpt)
{
    if (mbvptable.get() == &bvpt)   return;
    mbvptable.reset(bvpt.clone());
}


const BVParametersTable& BVSCalculator::getBVParamTable() const
{
    assert(mbvptable.get());
    return *mbvptable;
}


void BVSCalculator::setValencePrecision(double eps)
{
    mvalenceprecision = eps;
}


double BVSCalculator::getValencePrecision() const
{
    return mvalenceprecision;
}

// Protected Methods ---------------------------------------------------------

// PairQuantity overloads

void BVSCalculator::resetValue()
{
    // calcPoints requires that structure and rlimits data are cached.
    this->cacheStructureData();
    this->resizeValue(this->mstructure->countSites());
    this->PairQuantity::resetValue();
}


void BVSCalculator::configureBondGenerator(BaseBondGenerator& bnds)
{
    double rmaxused = min(this->getRmax(),
            this->rmaxFromPrecision(this->getValencePrecision()));
    bnds.setRmax(rmaxused);
}


void BVSCalculator::addPairContribution(const BaseBondGenerator& bnds)
{
    int summationscale = (bnds.site0() == bnds.site1()) ? 1 : 2;
    const string& a0 = mstructure_cache.baresymbols[bnds.site0()];
    const string& a1 = mstructure_cache.baresymbols[bnds.site1()];
    int v0 = mstructure_cache.valences[bnds.site0()];
    int v1 = mstructure_cache.valences[bnds.site1()];
    const BVParametersTable& bvtb = this->getBVParamTable();
    const BVParam& bp = bvtb.lookup(a0, v0, a1, v1);
    // do nothing if there are no bond parameters for this pair
    if (&bp == &bvtb.none())    return;
    double valencehalf = bp.bondvalence(bnds.distance()) / 2.0;
    int plusminus0 = (v0 >= 0) ? 1 : -1;
    int plusminus1 = (v1 >= 0) ? 1 : -1;
    mvalue[bnds.site0()] += summationscale * plusminus0 * valencehalf;
    mvalue[bnds.site1()] += summationscale * plusminus1 * valencehalf;
}

// Private Methods -----------------------------------------------------------

void BVSCalculator::cacheStructureData()
{
    int cntsites = this->countSites();
    mstructure_cache.baresymbols.resize(cntsites);
    mstructure_cache.valences.resize(cntsites);
    mstructure_cache.multiplicities.resize(cntsites);
    for (int i = 0; i < cntsites; ++i)
    {
        const string& smbl = mstructure->siteAtomType(i);
        mstructure_cache.baresymbols[i] = atomBareSymbol(smbl);
        mstructure_cache.valences[i] = atomValence(smbl);
        mstructure_cache.multiplicities[i] = mstructure->siteMultiplicity(i);
        if (mstructure->siteOccupancy(i) != 1.0)
        {
            const char* emsg = "Non unit occupancy is not supported.";
            throw invalid_argument(emsg);
        }
    }
    mstructure_cache.total_occupancy = mstructure->totalOccupancy();
}


double BVSCalculator::rmaxFromPrecision(double eps) const
{
    const BVParametersTable& bptb = this->getBVParamTable();
    BVParametersTable::SetOfBVParam bpused;
    for (int i0 = 0; i0 < this->countSites(); ++i0)
    {
        for (int i1 = i0; i1 < this->countSites(); ++i1)
        {
            const string& a0 = mstructure_cache.baresymbols[i0];
            const string& a1 = mstructure_cache.baresymbols[i1];
            int v0 = mstructure_cache.valences[i0];
            int v1 = mstructure_cache.valences[i1];
            const BVParam& bp = bptb.lookup(a0, v0, a1, v1);
            bpused.insert(bp);
        }
    }
    double rv = 0.0;
    BVParametersTable::SetOfBVParam::iterator bpit;
    for (bpit = bpused.begin(); bpit != bpused.end(); ++bpit)
    {
        rv = max(rv, bpit->bondvalenceToDistance(eps));
    }
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
