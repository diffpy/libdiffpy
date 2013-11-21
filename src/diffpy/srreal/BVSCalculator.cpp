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
*****************************************************************************/

#include <cmath>
#include <cassert>

#include <diffpy/validators.hpp>
#include <diffpy/serialization.ipp>
#include <diffpy/srreal/AtomUtils.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>

using namespace std;
using namespace diffpy::validators;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

BVSCalculator::BVSCalculator()
{
    // default configuration
    const double valence_precision = 1e-5;
    BVParametersTablePtr bvtb(new BVParametersTable);
    this->setBVParamTable(bvtb);
    this->setValencePrecision(valence_precision);
    this->setStructure(mstructure);
    // attributes
    this->registerDoubleAttribute("valenceprecision", this,
            &BVSCalculator::getValencePrecision,
            &BVSCalculator::setValencePrecision);
    this->registerDoubleAttribute("rmaxused", this,
            &BVSCalculator::getRmaxUsed);
}

// Public Methods ------------------------------------------------------------

// results

QuantityType BVSCalculator::valences() const
{
    QuantityType rv(mstructure_cache.valences.begin(),
            mstructure_cache.valences.end());
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
    int cntsites = this->countSites();
    assert(int(bd.size()) == cntsites);
    double sumofsquares = 0.0;
    for (int i = 0; i < cntsites; ++i)
    {
        sumofsquares += mstructure->siteMultiplicity(i) *
            mstructure->siteOccupancy(i) * bd[i] * bd[i];
    }
    double totocc = mstructure->totalOccupancy();
    double rv = (totocc > 0.0) ? (sumofsquares / totocc) : 0.0;
    return rv;
}


double BVSCalculator::bvrmsdiff() const
{
    double rvsq = this->bvmsdiff();
    return sqrt(rvsq);
}


void BVSCalculator::setBVParamTable(BVParametersTablePtr bvtb)
{
    ensureNonNull("BVParametersTable", bvtb);
    mbvptable = bvtb;
}


BVParametersTablePtr& BVSCalculator::getBVParamTable()
{
    return mbvptable;
}


const BVParametersTablePtr& BVSCalculator::getBVParamTable() const
{
    return mbvptable;
}


void BVSCalculator::setValencePrecision(double eps)
{
    ensureEpsilonPositive("valenceprecision", eps);
    mvalenceprecision = eps;
}


double BVSCalculator::getValencePrecision() const
{
    return mvalenceprecision;
}


double BVSCalculator::getRmaxUsed() const
{
    double rv = min(this->getRmax(),
            this->rmaxFromPrecision(this->getValencePrecision()));
    return rv;
}

// Protected Methods ---------------------------------------------------------

// PairQuantity overloads

void BVSCalculator::resetValue()
{
    // structure data need to be cached for rmaxFromPrecision
    this->cacheStructureData();
    this->resizeValue(this->countSites());
    this->PairQuantity::resetValue();
}


void BVSCalculator::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmax(this->getRmaxUsed());
}


void BVSCalculator::addPairContribution(const BaseBondGenerator& bnds,
        int summationscale)
{
    const string& a0 = mstructure_cache.baresymbols[bnds.site0()];
    const string& a1 = mstructure_cache.baresymbols[bnds.site1()];
    int v0 = mstructure_cache.valences[bnds.site0()];
    int v1 = mstructure_cache.valences[bnds.site1()];
    const BVParametersTable& bvtb = *(this->getBVParamTable());
    const BVParam& bp = bvtb.lookup(a0, v0, a1, v1);
    // do nothing if there are no bond parameters for this pair
    if (&bp == &bvtb.none())    return;
    double valencehalf = bp.bondvalence(bnds.distance()) / 2.0;
    int pm0 = (v0 >= 0) ? 1 : -1;
    int pm1 = (v1 >= 0) ? 1 : -1;
    const double& o0 = mstructure->siteOccupancy(bnds.site0());
    const double& o1 = mstructure->siteOccupancy(bnds.site1());
    mvalue[bnds.site0()] += summationscale * pm0 * valencehalf * o1;
    mvalue[bnds.site1()] += summationscale * pm1 * valencehalf * o0;
}

// Private Methods -----------------------------------------------------------

void BVSCalculator::cacheStructureData()
{
    int cntsites = this->countSites();
    mstructure_cache.baresymbols.resize(cntsites);
    mstructure_cache.valences.resize(cntsites);
    for (int i = 0; i < cntsites; ++i)
    {
        const string& smbl = mstructure->siteAtomType(i);
        mstructure_cache.baresymbols[i] = atomBareSymbol(smbl);
        mstructure_cache.valences[i] = atomValence(smbl);
    }
}


double BVSCalculator::rmaxFromPrecision(double eps) const
{
    const BVParametersTable& bvtb = *(this->getBVParamTable());
    // build a set of unique (symbol, valence) pairs
    typedef boost::unordered_set< pair<string, int> > SymbolValenceSet;
    SymbolValenceSet allsymvals;
    for (int i = 0; i < this->countSites(); ++i)
    {
        allsymvals.insert(make_pair(mstructure_cache.baresymbols[i],
                    mstructure_cache.valences[i]));
    }
    // find all used bond parameters
    BVParametersTable::SetOfBVParam bpused;
    SymbolValenceSet::const_iterator sv0, sv1;
    for (sv0 = allsymvals.begin(); sv0 != allsymvals.end(); ++sv0)
    {
        for (sv1 = sv0; sv1 != allsymvals.end(); ++sv1)
        {
            const BVParam& bp = bvtb.lookup(
                    sv0->first, sv0->second, sv1->first, sv1->second);
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

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::BVSCalculator)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::BVSCalculator)

// End of file
