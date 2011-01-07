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

* class PQEvaluatorOptimized -- optimized PairQuantity evaluator with fast
*     quantity updates
*
* $Id$
*
*****************************************************************************/


#include <stdexcept>
#include <sstream>
#include <boost/unordered_map.hpp>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>

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

PQEvaluatorType PQEvaluatorBasic::typeint() const
{
    return BASIC;
}


void PQEvaluatorBasic::updateValue(PairQuantity& pq)
{
    this->updateTypeMaskData(pq);
    auto_ptr<BaseBondGenerator> bnds;
    bnds.reset(pq.mstructure->createBondGenerator());
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
        bnds->selectSiteRange(0, i0 + 1);
        for (bnds->rewind(); !bnds->finished(); bnds->next())
        {
            if (chop_inner && (n++ % mncpu))    continue;
            int i1 = bnds->site1();
            bool exclude = !pq.getPairMask(i0, i1) ||
                !(this->getTypeMaskOfSites(i0, i1));
            if (exclude)  continue;
            int summationscale = (i0 == i1) ? 1 : 2;
            pq.addPairContribution(*bnds, summationscale);
        }
    }
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

// Protected Methods ---------------------------------------------------------

bool PQEvaluatorBasic::getTypeMaskOfSites(int i, int j) const
{
    assert(0 <= i && i < int(mtypeindex.size()));
    assert(0 <= j && j < int(mtypeindex.size()));
    int tpi = mtypeindex[i];
    int tpj = mtypeindex[j];
    assert(tpi < mtypemaskmatrix.rows());
    assert(tpj < mtypemaskmatrix.columns());
    bool rv = mtypemaskmatrix(tpi, tpj);
    return rv;
}


void PQEvaluatorBasic::updateTypeMaskData(const PairQuantity& pq)
{
    boost::unordered_map<string, int> tpidx;
    int cntsites = pq.countSites();
    mtypeindex.resize(cntsites);
    for (int i = 0; i < cntsites; ++i)
    {
        const string& smbl = pq.mstructure->siteAtomType(i);
        tpidx.emplace(smbl, tpidx.size());
        mtypeindex[i] = tpidx[smbl];
    }
    int cnttypes = tpidx.size();
    mtypemaskmatrix.resize(cnttypes, cnttypes);
    boost::unordered_map<string, int>::const_iterator smbi, smbj;
    for (smbi = tpidx.begin(); smbi != tpidx.end(); ++smbi)
    {
        for (smbj = smbi; smbj != tpidx.end(); ++smbj)
        {
            bool tpmask = pq.getTypeMask(smbi->first, smbj->first);
            mtypemaskmatrix(smbi->second, smbj->second) = tpmask;
            mtypemaskmatrix(smbj->second, smbi->second) = tpmask;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
// class PQEvaluatorOptimized
//////////////////////////////////////////////////////////////////////////////

PQEvaluatorType PQEvaluatorOptimized::typeint() const
{
    return OPTIMIZED;
}


void PQEvaluatorOptimized::updateValue(PairQuantity& pq)
{
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
