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
#include <diffpy/srreal/BaseBondGenerator.hpp>
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

const QuantityType& PairQuantity::eval(StructureAdapterPtr stru)
{
    mstructure = stru;
    mcountsites = mstructure->countSites();
    mstructure->customPQConfig(this);
    mevaluator->updateValue(*this);
    return this->value();
}


const QuantityType& PairQuantity::value() const
{
    return mvalue;
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
    mevaluator.reset(createPQEvaluator(evtp));
    this->resetValue();
}


int PairQuantity::countSites() const
{
    return mcountsites;
}


void PairQuantity::maskAllPairs(bool mask)
{
    minvertpairmask.clear();
    mdefaultpairmask = mask;
}


void PairQuantity::maskSitePair(int i, int j, bool mask)
{
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    if (mask == mdefaultpairmask)  minvertpairmask.erase(ij);
    else    minvertpairmask.insert(ij);
}


bool PairQuantity::getPairMask(int i, int j)
{
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    bool rv = minvertpairmask.count(ij) ?
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


void PairQuantity::configureBondGenerator(BaseBondGenerator& bnds)
{
    bnds.setRmin(this->getRmin());
    bnds.setRmax(this->getRmax());
}

}   // namespace srreal
}   // namespace diffpy

// End of file
