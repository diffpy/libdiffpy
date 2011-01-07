/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2008 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class PairQuantity -- general implementation of pair quantity calculator
*
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/mathutils.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

PairQuantity::PairQuantity()
{
    this->setRmin(0.0);
    this->setRmax(DEFAULT_BONDGENERATOR_RMAX);
    this->setEvaluator(BASIC);
    mcountsites = 0;
    this->maskAllPairs(true);
    // attributes
    this->registerDoubleAttribute("rmin", this,
            &PairQuantity::getRmin, &PairQuantity::setRmin);
    this->registerDoubleAttribute("rmax", this,
            &PairQuantity::getRmax, &PairQuantity::setRmax);
}

// Public Methods ------------------------------------------------------------

void PairQuantity::setStructure(StructureAdapterPtr stru)
{
    mstructure = stru;
    mcountsites = mstructure->countSites();
    mstructure->customPQConfig(this);
    this->resetValue();
}


const QuantityType& PairQuantity::value() const
{
    return mvalue;
}


void PairQuantity::mergeParallelValue(const QuantityType& pvalue)
{
    if (pvalue.size() != mvalue.size())
    {
        throw invalid_argument("Merged value array must have the same size.");
    }
    transform(mvalue.begin(), mvalue.end(), pvalue.begin(),
            mvalue.begin(), plus<double>());
}


void PairQuantity::setRmin(double rmin)
{
    mrmin = rmin;
}


const double& PairQuantity::getRmin() const
{
    return mrmin;
}


void PairQuantity::setRmax(double rmax)
{
    mrmax = rmax;
}


const double& PairQuantity::getRmax() const
{
    return mrmax;
}


void PairQuantity::setEvaluator(PQEvaluatorType evtp)
{
    if (mevaluator.get() && mevaluator->typeint() == evtp)  return;
    mevaluator = createPQEvaluator(evtp);
    this->resetValue();
}


void PairQuantity::setupParallelRun(int cpuindex, int ncpu)
{
    mevaluator->setupParallelRun(cpuindex, ncpu);
}


int PairQuantity::countSites() const
{
    return mcountsites;
}


void PairQuantity::maskAllPairs(bool mask)
{
    minvertpairmask.clear();
    minverttypemask.clear();
    mdefaultpairmask = mask;
}


void PairQuantity::invertMask()
{
    mdefaultpairmask = !mdefaultpairmask;
}


void PairQuantity::setPairMask(int i, int j, bool mask)
{
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    if (mask == mdefaultpairmask)  minvertpairmask.erase(ij);
    else    minvertpairmask.insert(ij);
}


bool PairQuantity::getPairMask(int i, int j) const
{
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    bool rv = minvertpairmask.count(ij) ?
        !mdefaultpairmask : mdefaultpairmask;
    return rv;
}


void PairQuantity::
setTypeMask(const string& smbli, const string& smblj, bool mask)
{
    pair<string,string> smblij = (smbli > smblj) ?
        make_pair(smblj, smbli) : make_pair(smbli, smblj);
    if (mask == mdefaultpairmask)  minverttypemask.erase(smblij);
    else    minverttypemask.insert(smblij);
}


bool PairQuantity::getTypeMask(const string& smbli, const string& smblj) const
{
    pair<string,string> smblij = (smbli > smblj) ?
        make_pair(smblj, smbli) : make_pair(smbli, smblj);
    bool rv = minverttypemask.count(smblij) ?
        !mdefaultpairmask : mdefaultpairmask;
    return rv;
}

// Protected Methods ---------------------------------------------------------

void PairQuantity::resizeValue(size_t sz)
{
    mvalue.resize(sz);
}


void PairQuantity::resetValue()
{
    fill(mvalue.begin(), mvalue.end(), 0.0);
}


void PairQuantity::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmin(this->getRmin());
    bnds.setRmax(this->getRmax());
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT(diffpy::srreal::PairQuantity)

// End of file
