/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class PQEvaluatorBasic -- robust PairQuantity evaluator, the result
*     is always calculated from scratch.
*
* class PQEvaluatorOptimized -- optimized PairQuantity evaluator with fast
*     quantity updates
*
*****************************************************************************/


#include <stdexcept>
#include <sstream>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/StructureDifference.hpp>

using namespace std;

// Local Constants -----------------------------------------------------------

namespace {
// tolerated load variance for splitting outer loop for parallel evaluation
const double CPU_LOAD_VARIANCE = 0.1;
}

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class PQEvaluatorBasic
//////////////////////////////////////////////////////////////////////////////

PQEvaluatorBasic::PQEvaluatorBasic() :
    musefullsum(false),
    mcpuindex(0), mncpu(1)
{ }


PQEvaluatorType PQEvaluatorBasic::typeint() const
{
    return BASIC;
}


void PQEvaluatorBasic::updateValue(
        PairQuantity& pq, StructureAdapterConstPtr stru)
{
    mtypeused = BASIC;
    pq.setStructure(stru);
    BaseBondGeneratorPtr bnds = pq.mstructure->createBondGenerator();
    pq.configureBondGenerator(*bnds);
    int cntsites = pq.mstructure->countSites();
    // loop counter
    long n = mcpuindex;
    // split outer loop for many atoms.  The CPUs should have similar load.
    bool chop_outer = (mncpu <= ((cntsites - 1) * CPU_LOAD_VARIANCE + 1));
    bool chop_inner = !chop_outer;
    for (int i0 = 0; i0 < cntsites; ++i0)
    {
        if (chop_outer && (n++ % mncpu))    continue;
        bnds->selectAnchorSite(i0);
        int i1hi = musefullsum ? cntsites : (i0 + 1);
        bnds->selectSiteRange(0, i1hi);
        for (bnds->rewind(); !bnds->finished(); bnds->next())
        {
            if (chop_inner && (n++ % mncpu))    continue;
            int i1 = bnds->site1();
            if (!pq.getPairMask(i0, i1))   continue;
            int summationscale = (musefullsum || i0 == i1) ? 1 : 2;
            pq.addPairContribution(*bnds, summationscale);
        }
    }
    mvalue_ticker.click();
}


void PQEvaluatorBasic::useFullSum(bool flag)
{
    musefullsum = flag;
}


void PQEvaluatorBasic::setupParallelRun(int cpuindex, int ncpu)
{
    // make sure ncpu is at least one
    if (ncpu < 1)
    {
        const char* emsg = "Number of CPU ncpu must be at least 1.";
        throw invalid_argument(emsg);
    }
    mcpuindex = cpuindex;
    mncpu = ncpu;
}

//////////////////////////////////////////////////////////////////////////////
// class PQEvaluatorOptimized
//////////////////////////////////////////////////////////////////////////////

PQEvaluatorType PQEvaluatorOptimized::typeint() const
{
    return OPTIMIZED;
}


void PQEvaluatorOptimized::updateValue(
        PairQuantity& pq, StructureAdapterConstPtr stru)
{
    mtypeused = OPTIMIZED;
    // revert to normal calculation if there is no structure or
    // if PairQuantity uses mask
    if (pq.ticker() >= mvalue_ticker || !pq.getStructure() || pq.hasMask())
    {
        this->PQEvaluatorBasic::updateValue(pq, stru);
        return;
    }
    // do not do fast updates if they take more work
    StructureDifference sd = pq.getStructure()->diff(stru);
    if (!sd.allowsfastupdate())
    {
        this->PQEvaluatorBasic::updateValue(pq, stru);
        return;
    }
    // Remove contributions from the extra sites in the old structure
    assert(sd.stru0 == pq.mstructure);
    int cntsites0 = sd.stru0->countSites();
    BaseBondGeneratorPtr bnds0 = sd.stru0->createBondGenerator();
    // loop counter
    long n = mcpuindex;
    // the loop is adjusted to always perform a full sum and to split the
    // inner loop in case of parallel calculation.
    bnds0->selectSiteRange(0, cntsites0);
    std::vector<int>::const_iterator ii0;
    for (ii0 = sd.pop0.begin(); ii0 != sd.pop0.end(); ++ii0)
    {
        const int& i0 = *ii0;
        bnds0->selectAnchorSite(i0);
        for (bnds0->rewind(); !bnds0->finished(); bnds0->next())
        {
            if (n++ % mncpu)    continue;
            int i1 = bnds0->site1();
            assert(pq.getPairMask(i0, i1));
            const int summationscale = (i0 == i1) ? -1 : -2;
            pq.addPairContribution(*bnds0, summationscale);
        }
        bnds0->selectSite(i0, false);
    }
    // Add contributions from the new atoms in the updated structure
    // save current value to override the resetValue call from setStructure
    assert(sd.stru1);
    pq.stashPartialValue();
    pq.setStructure(sd.stru1);
    pq.restorePartialValue();
    int cntsites1 = sd.stru1->countSites();
    BaseBondGeneratorPtr bnds1 = sd.stru1->createBondGenerator();
    bnds1->selectSiteRange(0, cntsites1);
    std::vector<int>::const_iterator ii1;
    for (ii1 = sd.add1.begin(); ii1 != sd.add1.end(); ++ii1)
    {
        bnds1->selectSite(*ii1, false);
    }
    for (ii1 = sd.add1.begin(); ii1 != sd.add1.end(); ++ii1)
    {
        const int& i0 = *ii1;
        bnds1->selectAnchorSite(i0);
        bnds1->selectSite(i0, true);
        for (bnds1->rewind(); !bnds1->finished(); bnds1->next())
        {
            if (n++ % mncpu)    continue;
            int i1 = bnds1->site1();
            assert(pq.getPairMask(i0, i1));
            const int summationscale = (i0 == i1) ? +1 : +2;
            pq.addPairContribution(*bnds1, summationscale);
        }
    }
    mvalue_ticker.click();
}

// Factory for PairQuantity evaluators ---------------------------------------

PQEvaluatorPtr createPQEvaluator(PQEvaluatorType pqtp)
{
    PQEvaluatorPtr rv;
    switch (pqtp)
    {
        case BASIC:
            rv.reset(new PQEvaluatorBasic());
            break;

        case OPTIMIZED:
            rv.reset(new PQEvaluatorOptimized());
            break;

        default:
            ostringstream emsg;
            emsg << "Invalid PQEvaluatorType value " << pqtp;
            throw invalid_argument(emsg.str());
    }
    return rv;
}


}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::PQEvaluatorBasic)
BOOST_CLASS_EXPORT(diffpy::srreal::PQEvaluatorBasic)

// End of file
