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

using namespace std;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

BVSCalculator::BVSCalculator()
{
    // default configuration
    // FIXME: this->setBVParamTable();
    this->setRmax(10.0);
}

// Public Methods ------------------------------------------------------------

// results

double BVSCalculator::bvmsd() const
{
    // FIXME
    return 0.0;
}


double BVSCalculator::bvrmsd() const
{
    double rvsq = this->bvmsd();
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


const BVParam& BVSCalculator::bvpar(int idx0, int idx1) const
{
    // FIXME
    static BVParam rv;
    return rv;
}

// Protected Methods ---------------------------------------------------------

void BVSCalculator::resetValue()
{ }


void BVSCalculator::configureBondGenerator(BaseBondGenerator& bnds)
{ }


void BVSCalculator::addPairContribution(const BaseBondGenerator& bnds)
{ }


}   // namespace srreal
}   // namespace diffpy

// End of file
