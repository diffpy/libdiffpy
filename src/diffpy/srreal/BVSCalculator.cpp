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

QuantityType BVSCalculator::bvdiff() const
{
    assert(mstructure_cache.valences.size() == this->value().size());
    int cntsites = this->value().size();
    QuantityType rv(cntsites);
    for (int i = 0; i < cntsites; ++i)
    {
        rv[i] = mstructure_cache.valences[i] - this->value()[i];
    }
    return rv;
}


double BVSCalculator::bvmsdiff() const
{
    QuantityType bd = this->bvdiff();
    assert(bd.size() == mstructure_cache.multiplicities.size());
    int cntsites = bd.size();
    double sumofsquares = 0.0;
    for (int i = 0; i < cntsites; ++i)
    {
        sumofsquares += mstructure_cache.multiplicities[i] * pow(bd[i], 2);
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
    mbvptable.reset(bvpt.copy());
}


const BVParametersTable& BVSCalculator::getBVParamTable() const
{
    assert(mbvptable.get());
    return *mbvptable;
}


void BVSCalculator::setValencePrecision(double eps)
{
    BVParametersTable::SetOfBVParam allbp = this->getBVParamTable().getAll();
    double dmx = 0.0;
    BVParametersTable::SetOfBVParam::const_iterator bpit;
    for (bpit = allbp.begin(); bpit != allbp.end(); ++bpit)
    {
        dmx = max(dmx, bpit->bondvalenceToDistance(eps));
    }
    this->setRmax(dmx);
}


double BVSCalculator::getValencePrecision() const
{
    BVParametersTable::SetOfBVParam allbp = this->getBVParamTable().getAll();
    double veps = 0.0;
    BVParametersTable::SetOfBVParam::const_iterator bpit;
    for (bpit = allbp.begin(); bpit != allbp.end(); ++bpit)
    {
        veps = max(veps, bpit->bondvalence(this->getRmax()));
    }
    return veps;
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
{ }


void BVSCalculator::addPairContribution(const BaseBondGenerator& bnds)
{
    int summationscale = (bnds.site0() == bnds.site1()) ? 1 : 2;
    const string& a0 = mstructure_cache.baresymbols[bnds.site0()];
    const string& a1 = mstructure_cache.baresymbols[bnds.site1()];
    int v0 = mstructure_cache.valences[bnds.site0()];
    int v1 = mstructure_cache.valences[bnds.site1()];
    const BVParam& bp = mbvptable->lookup(a0, v0, a1, v1);
    double valence = bp.bondvalence(bnds.distance());
    int plusminus0 = (v0 >= 0) ? 1 : -1;
    int plusminus1 = (v1 >= 0) ? 1 : -1;
    mvalue[bnds.site0()] += summationscale * plusminus0 * valence;
    mvalue[bnds.site1()] += summationscale * plusminus1 * valence;
}

// Private Methods -----------------------------------------------------------

void BVSCalculator::cacheStructureData()
{
    int cntsites = mstructure->countSites();
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

}   // namespace srreal
}   // namespace diffpy

// End of file
